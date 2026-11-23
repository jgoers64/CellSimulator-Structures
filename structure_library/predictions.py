"""Explicit, bounded AlphaFold retrieval for a selected nuclear protein.

This module does not choose predictions for the entire inventory. API-returned
URLs, versions and file hashes are pinned in ignored cache. A source match and a
high confidence score do not establish experimental accuracy or MD readiness.
"""
from __future__ import annotations

import json
import math
import re
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlparse

from .coordinates import (
    atomic_write, bounded_decompress, category_rows, digest, download,
    fetch_source, read_json, require_gemmi, safe_path, write_json,
)


MAX_METADATA_BYTES = 2_000_000
MAX_COORDINATE_BYTES = 40_000_000
MAX_PAE_BYTES = 32_000_000
MAX_PAE_CELLS = 4_000_000
RESOLVER_VERSION = "alphafold-nuclear-v1"
SAFE_ID = re.compile(r"[a-z0-9][a-z0-9_-]*")
ACCESSION = re.compile(
    r"(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-[0-9]+)?"
)


def _model_id(record):
    return record.get("modelEntityId") or record.get("entryId")


def _positive_int(value, label):
    if value is None:
        return None
    if isinstance(value, bool):
        raise ValueError(f"Invalid {label}")
    number = int(value)
    if number <= 0 or str(number) != str(value):
        raise ValueError(f"Invalid {label}: {value!r}")
    return number


def _af_file_url(url):
    if not isinstance(url, str):
        raise ValueError("AlphaFold API did not supply a coordinate/evidence URL")
    parsed = urlparse(url)
    if (parsed.scheme != "https" or parsed.hostname != "alphafold.ebi.ac.uk"
            or not parsed.path.startswith("/files/") or parsed.username or parsed.password
            or parsed.fragment or parsed.port not in (None, 443)):
        raise ValueError("AlphaFold returned an unsupported source URL; review it before importing")
    return url


def _get_file(url, budget):
    raw = download(_af_file_url(url), max_bytes=budget)
    if urlparse(url).path.endswith(".gz"):
        raw = bounded_decompress(raw, max_bytes=budget)
    return raw


def _reference_length(entry, accession):
    values = {
        item["sequence_length"] for item in entry.get("source_lookups", [])
        if item.get("uniprot_accession") == accession and item.get("sequence_length") is not None
    }
    if len(values) > 1:
        raise ValueError("Registry contains conflicting reference sequence lengths")
    return _positive_int(next(iter(values)), "reference sequence length") if values else None


def _sequence_metadata(record, reference_length, explicit_model):
    sequence = record.get("sequence", record.get("uniprotSequence"))
    if not isinstance(sequence, str) or not re.fullmatch(r"[A-Z]+", sequence):
        raise ValueError("AlphaFold metadata lacks a usable explicit protein sequence")
    start = _positive_int(record.get("sequenceStart", record.get("uniprotStart")), "sequence start")
    end = _positive_int(record.get("sequenceEnd", record.get("uniprotEnd")), "sequence end")
    if (start is None) != (end is None):
        raise ValueError("AlphaFold supplied an incomplete sequence range")
    if start is not None and (end < start or end - start + 1 != len(sequence)):
        raise ValueError("AlphaFold sequence length and alignment range disagree; manual review is required")
    if reference_length is not None and end is not None and end > reference_length:
        raise ValueError("Prediction range exceeds the registry reference length; review sequence-version differences")
    fragment = None
    if start is not None:
        if start > 1 or (reference_length is not None and end < reference_length):
            fragment = True
        elif reference_length is not None and start == 1 and end == reference_length:
            fragment = False
    suffix = re.search(r"-F([0-9]+)$", _model_id(record) or "")
    if suffix and int(suffix.group(1)) > 1:
        fragment = True
    if (fragment is True or start is None) and not explicit_model:
        raise ValueError(f"Prediction is a fragment or has unresolved sequence bounds; select model_id={_model_id(record)!r} explicitly")
    return {
        "api_sequence": sequence,
        "api_sequence_sha256": digest(sequence.encode("ascii")),
        "sequence_start": start, "sequence_end": end,
        "sequence_length": len(sequence), "reference_sequence_length": reference_length,
        "reference_length_source": "pinned UniProt registry snapshot" if reference_length else None,
        "reference_coverage_fraction": (len(sequence) / reference_length
                                        if reference_length and start is not None else None),
        "fragment": fragment,
        "reference_sequence_identity": "not independently sequence-aligned to current UniProt; accession and available length checked",
        "coverage_note": "Sequence coverage does not mean every atom, proteoform, biological assembly or conformation is correct.",
    }


def _coordinate_evidence(raw, record, accession, sequence):
    gemmi = require_gemmi()
    block = gemmi.cif.read_string(raw.decode("utf-8")).sole_block()
    if block.find_value("_entry.id") != _model_id(record):
        raise ValueError("AlphaFold coordinate entry ID does not match its API metadata")
    polymers = category_rows(block, "_entity_poly.")
    if len(polymers) != 1 or polymers[0].get("type") != "polypeptide(L)":
        raise ValueError("Resolver supports a single explicit protein polymer; other assemblies need separate curation")
    cif_sequence = "".join(str(polymers[0].get("pdbx_seq_one_letter_code_can", "")).split())
    if cif_sequence != sequence["api_sequence"]:
        raise ValueError("Coordinate polymer sequence differs from AlphaFold API sequence")
    references = category_rows(block, "_ma_target_ref_db_details.")
    accessions = {r.get("db_accession") for r in references if r.get("db_name") in ("UNP", "UniProt")}
    if not accessions:
        accessions = {r.get("pdbx_db_accession") for r in category_rows(block, "_struct_ref.")
                      if r.get("db_name") in ("UNP", "UniProt")}
    if accessions != {accession}:
        raise ValueError("Coordinate source does not identify exactly the requested UniProt accession")
    for reference in references:
        taxon = reference.get("ncbi_taxonomy_id")
        if taxon is not None and _positive_int(taxon, "coordinate taxonomy") != 9606:
            raise ValueError("Coordinate source taxonomy is not human")
        begin = _positive_int(reference.get("seq_db_align_begin"), "coordinate alignment start")
        end = _positive_int(reference.get("seq_db_align_end"), "coordinate alignment end")
        if begin is not None and sequence["sequence_start"] is not None and begin != sequence["sequence_start"]:
            raise ValueError("Coordinate and API sequence starts disagree")
        if end is not None and sequence["sequence_end"] is not None and end != sequence["sequence_end"]:
            raise ValueError("Coordinate and API sequence ends disagree")
    structure = gemmi.make_structure_from_block(block)
    if len(structure) != 1 or len(structure[0]) != 1:
        raise ValueError("Resolver supports one predicted model and one chain; ensembles/complexes need explicit curation")
    chain = structure[0][0]
    if record.get("chainId") and chain.name != record["chainId"]:
        raise ValueError("Predicted coordinate chain differs from API chain")
    expected_residues = {int(r["num"]): r["mon_id"] for r in category_rows(block, "_entity_poly_seq.")}
    positions, count = set(), 0
    plddt_low, plddt_high = 100.0, 0.0
    unresolved_positions = 0
    for residue in chain:
        if residue.label_seq is None:
            unresolved_positions += 1
        else:
            if not 1 <= residue.label_seq <= sequence["sequence_length"]:
                raise ValueError("Coordinate residue lies outside the API model sequence")
            if expected_residues.get(residue.label_seq) != residue.name:
                raise ValueError("Coordinate residue identity differs from its polymer sequence")
            positions.add(residue.label_seq)
        for atom in residue:
            if not all(math.isfinite(x) for x in (atom.pos.x, atom.pos.y, atom.pos.z, atom.b_iso)):
                raise ValueError("Prediction contains nonfinite coordinates or confidence")
            if not 0 <= atom.b_iso <= 100:
                raise ValueError("Predicted pLDDT field is outside its 0-100 scale")
            count += 1
            plddt_low = min(plddt_low, atom.b_iso)
            plddt_high = max(plddt_high, atom.b_iso)
    if count == 0:
        raise ValueError("Prediction contains no atoms")
    return {
        "api_and_coordinate_sequences_match": True,
        "coordinate_accession_matches": True,
        "coordinate_atom_count": count,
        "coordinate_residues_with_sequence_positions": len(positions),
        "coordinate_residues_without_sequence_positions": unresolved_positions,
        "coordinate_sequence_coverage_fraction": (len(positions) / sequence["sequence_length"]
                                                  if not unresolved_positions else None),
        "plddt_min_per_atom": plddt_low, "plddt_max_per_atom": plddt_high,
    }


def _pae_evidence(record, directory, root, length):
    url = record.get("paeDocUrl")
    if not url:
        return {"status": "not_provided_by_api", "source_url": None, "path": None, "sha256": None}
    _af_file_url(url)
    if length * length > MAX_PAE_CELLS:
        return {
            "status": "not_cached_size_budget", "source_url": url, "path": None, "sha256": None,
            "reason": f"PAE matrix would exceed {MAX_PAE_CELLS:,} cells; relative-domain confidence remains unreviewed.",
        }
    try:
        raw = _get_file(url, MAX_PAE_BYTES)
    except HTTPError as exc:
        if exc.code != 404:
            raise
        return {"status": "not_available_at_source", "source_url": url, "path": None,
                "sha256": None, "reason": "API-provided PAE URL returned HTTP 404 at resolution time"}
    data = json.loads(raw)
    if not isinstance(data, list) or len(data) != 1 or not isinstance(data[0], dict):
        raise ValueError("Unsupported PAE document shape; no confidence values inferred")
    matrix = data[0].get("predicted_aligned_error")
    if not isinstance(matrix, list) or len(matrix) != length:
        raise ValueError("PAE dimensions differ from the selected model sequence")
    for row in matrix:
        if not isinstance(row, list) or len(row) != length:
            raise ValueError("PAE matrix is not square for the selected sequence")
        if any(isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(x) or x < 0 for x in row):
            raise ValueError("PAE contains invalid error estimates")
    maximum = data[0].get("max_predicted_aligned_error")
    if maximum is not None and (isinstance(maximum, bool) or not isinstance(maximum, (int, float)) or not math.isfinite(maximum) or maximum < 0):
        raise ValueError("PAE maximum is invalid")
    path = directory / "pae.json"
    atomic_write(path, raw)
    return {
        "status": "cached", "source_url": url, "path": path.relative_to(root).as_posix(),
        "sha256": digest(raw), "matrix_shape": [length, length], "units": "angstrom",
        "reported_max_predicted_aligned_error": maximum,
        "interpretation": "Predicted relative-position uncertainty; not an experimental coordinate-error measurement.",
    }


def _cached_metadata(path, root, accession, model_id=None):
    if path.stat().st_size > MAX_METADATA_BYTES:
        raise ValueError("Cached prediction metadata exceeds its byte budget")
    metadata = read_json(path)
    source = metadata.get("source", {})
    if (metadata.get("provenance_class") != "predicted_model"
            or source.get("uniprot_accession") != accession
            or (model_id is not None and source.get("model_id") != model_id)):
        raise ValueError("Cached prediction does not match the selected identity")
    api_path = safe_path(root, source["api_response_path"])
    if not api_path.is_file() or api_path.stat().st_size > MAX_METADATA_BYTES or digest(api_path.read_bytes()) != source["api_response_sha256"]:
        raise ValueError("Cached API evidence is missing or changed; review the cached prediction before replacing it")
    fetch_source(metadata, root)  # Uses the original version/URL/hash; never silently refreshes it.
    pae = metadata.get("confidence", {}).get("pae", {})
    if pae.get("status") == "cached":
        pae_path = safe_path(root, pae["path"])
        if not pae_path.exists():
            raw = _get_file(pae["source_url"], MAX_PAE_BYTES)
            if digest(raw) != pae["sha256"]:
                raise ValueError("Pinned PAE document changed at the source")
            atomic_write(pae_path, raw)
        if pae_path.stat().st_size > MAX_PAE_BYTES or digest(pae_path.read_bytes()) != pae["sha256"]:
            raise ValueError("Cached PAE differs from its pinned source")
    return metadata


def resolve_prediction(entry, root, *, accession=None, model_id=None):
    """Resolve one explicit protein identity and pin a source-compatible view asset.

    Multiple accessions/models require selectors. Known fragments require an
    explicit model selector even when the API returned only one model. Removing
    an existing selection pin is a deliberate new-version review operation.
    """
    root = Path(root).resolve()
    if entry.get("kind") != "protein_gene" or not SAFE_ID.fullmatch(entry.get("id", "")):
        raise ValueError("Prediction resolution requires a registered protein_gene identity")
    accessions = entry.get("uniprot_accessions", [])
    if accession is None:
        if len(accessions) != 1:
            raise ValueError("Select one explicit accession from this gene entry: " + ", ".join(accessions or ["none recorded"]))
        accession = accessions[0]
    if not isinstance(accession, str) or not ACCESSION.fullmatch(accession) or accession not in accessions:
        raise ValueError("Selected accession is not a normalized identifier on this registry entry")
    if model_id is not None and (not isinstance(model_id, str) or not re.fullmatch(r"[A-Za-z0-9_-]+", model_id)):
        raise ValueError("Unsafe model selector")
    selection_key = f"{entry['id']}-{accession.lower()}-{model_id.lower() if model_id else 'default'}"
    selection_path = root / "cache" / "predictions" / "selections" / (selection_key + ".json")
    if selection_path.exists():
        if selection_path.stat().st_size > MAX_METADATA_BYTES:
            raise ValueError("Cached prediction selection exceeds its byte budget")
        pin = read_json(selection_path)
        if pin.get("inventory_entry_id") != entry["id"] or pin.get("uniprot_accession") != accession:
            raise ValueError("Prediction selection pin does not match its requested identity")
        cached = _cached_metadata(safe_path(root, pin["metadata_path"]), root, accession, model_id)
        previous_length = cached.get("sequence", {}).get("reference_sequence_length")
        current_length = _reference_length(entry, accession)
        if previous_length is not None and current_length is not None and previous_length != current_length:
            raise ValueError("Registry reference length changed since this prediction was pinned; review the sequence/version before reuse")
        return cached
    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{accession}"
    try:
        api_raw = download(api_url, max_bytes=MAX_METADATA_BYTES)
    except HTTPError as exc:
        if exc.code == 404:
            raise RuntimeError(f"AlphaFold has no API prediction for {accession}; no structure was generated") from exc
        raise RuntimeError(f"AlphaFold metadata request failed for {accession}: HTTP {exc.code}") from exc
    except (OSError, URLError) as exc:
        raise RuntimeError(f"AlphaFold metadata could not be retrieved for {accession}") from exc
    records = json.loads(api_raw)
    if not isinstance(records, list) or not all(isinstance(r, dict) for r in records):
        raise ValueError("Unexpected AlphaFold API response")
    if not records:
        raise RuntimeError(f"AlphaFold returned no predictions for {accession}; no coordinates were generated")
    if model_id is None and len(records) != 1:
        choices = ", ".join(str(_model_id(r)) for r in records)
        raise ValueError(f"Multiple AlphaFold models are available; select model_id explicitly: {choices}")
    matches = [r for r in records if model_id is None or _model_id(r) == model_id]
    if len(matches) != 1:
        raise ValueError("Selected model is absent or ambiguous in the AlphaFold API result")
    record = matches[0]
    selected_id = _model_id(record)
    if not isinstance(selected_id, str) or not re.fullmatch(r"[A-Za-z0-9_-]+", selected_id):
        raise ValueError("AlphaFold supplied an unusable model identifier")
    if record.get("uniprotAccession") != accession:
        raise ValueError("AlphaFold model accession differs from the selected registry accession")
    taxon = record.get("taxId")
    if taxon is not None and _positive_int(taxon, "API taxonomy") != 9606:
        raise ValueError("AlphaFold model taxonomy is not human")
    if record.get("isComplex") is True or record.get("entityType") not in (None, "protein"):
        raise ValueError("This resolver accepts individual protein predictions; complexes require curated assembly metadata")
    version = _positive_int(record.get("latestVersion"), "prediction version")
    if version is None:
        raise ValueError("AlphaFold did not provide a version to pin")
    sequence = _sequence_metadata(record, _reference_length(entry, accession), model_id is not None)
    ident = f"alphafold-{selected_id.lower()}-v{version}"
    directory = root / "cache" / "predictions" / ident
    metadata_path = directory / "metadata.json"
    if metadata_path.exists():
        metadata = _cached_metadata(metadata_path, root, accession, selected_id)
    else:
        coordinate_url = _af_file_url(record.get("cifUrl"))
        raw = _get_file(coordinate_url, MAX_COORDINATE_BYTES)
        evidence = _coordinate_evidence(raw, record, accession, sequence)
        pae = _pae_evidence(record, directory, root, sequence["sequence_length"])
        api_path = directory / "api-response.json"
        atomic_write(api_path, api_raw)
        coordinate_path = root / "cache" / "nuclear" / f"{ident}.cif"
        atomic_write(coordinate_path, raw)
        metadata = {
            "schema_version": 1, "id": ident,
            "name": f"{record.get('uniprotDescription') or entry['name']} (AlphaFold prediction)",
            "provenance_class": "predicted_model",
            "availability": "cached_coordinates",
            "source": {
                "database": "AlphaFold Protein Structure Database",
                "entry_url": f"https://alphafold.ebi.ac.uk/entry/{selected_id}",
                "api_url": api_url,
                "api_response_path": api_path.relative_to(root).as_posix(),
                "api_response_sha256": digest(api_raw),
                "coordinate_url": coordinate_url, "coordinate_sha256": digest(raw),
                "coordinate_path": coordinate_path.relative_to(root).as_posix(),
                "uniprot_accession": accession, "model_id": selected_id,
                "version": version, "model_created_date": record.get("modelCreatedDate"),
                "sequence_version_date": record.get("sequenceVersionDate"),
                "prediction_method": record.get("toolUsed"), "provider_id": record.get("providerId"),
            },
            "organism": {"taxonomy_id": taxon, "name": record.get("organismScientificName"),
                         "api_taxonomy_status": "human_verified" if taxon is not None else "not_provided"},
            "assembly": {"id": None, "description": "Single predicted chain; no physiological oligomer or complex is inferred."},
            "sequence": {**sequence, **evidence},
            "confidence": {
                "plddt_storage": "Original atom B_iso_or_equiv values retained as source_b_factor/plddt_per_atom by the importer",
                "api_global_metric_value": record.get("globalMetricValue"),
                "api_fraction_plddt_very_low": record.get("fractionPlddtVeryLow"),
                "api_fraction_plddt_low": record.get("fractionPlddtLow"),
                "api_fraction_plddt_confident": record.get("fractionPlddtConfident"),
                "api_fraction_plddt_very_high": record.get("fractionPlddtVeryHigh"),
                "pae": pae,
                "interpretation": "Model confidence is not experimental validation, an accuracy percentage, or evidence of a unique cellular conformation.",
            },
            "review": {"resolver": RESOLVER_VERSION, "identity_and_source_checks": "completed",
                       "scientific_structure_review": "not_completed", "scientific_accuracy_percent": None},
            "view_only": True, "simulation_ready": False,
            "limitations": [
                "Prediction of the selected sequence and state; not an experimentally measured structure.",
                "Gene isoforms, post-translational modifications, cofactors and oligomeric partners are not supplied automatically.",
                "Low-confidence regions and uncertain domain arrangements require review; PAE does not establish experimental correctness.",
                "No MD parameters, missing atoms, hydrogens, population counts or nuclear placement are generated.",
            ],
        }
        write_json(metadata_path, metadata)
    write_json(selection_path, {
        "inventory_entry_id": entry["id"], "uniprot_accession": accession,
        "model_id": metadata["source"]["model_id"],
        "metadata_path": metadata_path.relative_to(root).as_posix(),
        "policy": "Reuse pinned source/version; a later provider version requires deliberate new selection review.",
    })
    return metadata
