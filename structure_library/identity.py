"""Join nuclear evidence to existing protein identities without installing copies.

Gene annotations, protein records, coordinate representations and scene copies
are different things. This module only reads metadata and, when available, the
target simulator's existing Human Proteome database. It never creates a database,
downloads a structure, assigns populations, or rewrites an inventory record.
"""
from __future__ import annotations

from contextlib import closing
import json
from pathlib import Path
import re
import sqlite3


_CANONICAL_ACCESSION = (
    r"(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|"
    r"[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})"
)
_CANONICAL = re.compile(_CANONICAL_ACCESSION)
_WITH_ISOFORM = re.compile(_CANONICAL_ACCESSION + r"(?:-[1-9][0-9]*)?")
_VERIFIED = "reviewed_snapshot_match"


def protein_ref(accession: str) -> str:
    """Return the simulator's existing canonical UniProt reference.

    Isoform suffixes are deliberately rejected: the current simulator resolves
    canonical records only. Removing a suffix would change the requested entity.
    """
    if not isinstance(accession, str):
        raise ValueError("A UniProt accession must be a string")
    normalized = accession.strip().upper()
    if not _CANONICAL.fullmatch(normalized):
        raise ValueError(f"Expected a canonical UniProt accession: {accession!r}")
    return "human-proteome/" + normalized


def _text(value) -> str:
    return value.strip() if isinstance(value, str) else ""


def _evidence(entry: dict, accession: str | None) -> dict:
    source = entry.get("source_record") or {}
    return {
        "inventory_id": entry["id"],
        "gene_symbol": _text(entry.get("gene_symbol")),
        "ensembl_gene_id": _text(source.get("ensembl_gene_id")),
        "source_url": _text(source.get("source_url")),
        "source_release": _text(source.get("source_release")),
        "evidence_scope": _text(source.get("evidence_scope")),
        "uniprot_url": (f"https://www.uniprot.org/uniprotkb/{accession}/entry"
                        if accession and _WITH_ISOFORM.fullmatch(accession) else None),
    }


def _unique_dicts(rows: list[dict]) -> list[dict]:
    # Stable ordering makes a repeated reconciliation produce the same report.
    by_value = {json.dumps(row, sort_keys=True, ensure_ascii=False): row for row in rows}
    return [by_value[key] for key in sorted(by_value)]


def _collection_links(entries) -> tuple[list[dict], list[dict], set[str]]:
    grouped: dict[str, dict] = {}
    unresolved: list[dict] = []
    gene_ids: set[str] = set()
    identity_by_id: dict[str, tuple] = {}
    for entry in entries:
        if not isinstance(entry, dict):
            raise ValueError("Inventory entries must be objects")
        if entry.get("kind") != "protein_gene":
            continue
        ident = _text(entry.get("id"))
        if not ident:
            raise ValueError("Protein-gene evidence needs a stable inventory ID")
        gene_ids.add(ident)
        accessions = entry.get("uniprot_accessions", [])
        lookups = entry.get("source_lookups", [])
        if not isinstance(accessions, list) or not isinstance(lookups, list):
            raise ValueError(f"Invalid accession/source lookup list for {ident}")
        if not all(isinstance(value, str) for value in accessions):
            raise ValueError(f"Accession values must be strings for {ident}")
        if not all(isinstance(value, dict) for value in lookups):
            raise ValueError(f"Source lookups must be objects for {ident}")
        normalized = sorted({_text(value).upper() for value in accessions})
        fingerprint = (_text(entry.get("gene_symbol")), tuple(normalized))
        if ident in identity_by_id and identity_by_id[ident] != fingerprint:
            raise ValueError(f"Conflicting identities use the same inventory ID: {ident}")
        identity_by_id[ident] = fingerprint
        if not normalized:
            unresolved.append({**_evidence(entry, None), "uniprot_accession": None,
                               "reason": "no_accession_on_inventory_entry"})
            continue
        for accession in normalized:
            reason = None
            if not _CANONICAL.fullmatch(accession):
                reason = ("isoform_requires_explicit_support"
                          if _WITH_ISOFORM.fullmatch(accession)
                          else "invalid_accession")
            statuses = {_text(row.get("mapping_status")) for row in lookups
                        if _text(row.get("uniprot_accession")).upper() == accession}
            if reason is None and statuses != {_VERIFIED}:
                reason = ("conflicting_source_mapping" if len(statuses) > 1 else
                          next(iter(statuses)) if statuses else "no_source_mapping")
                reason = reason or "no_source_mapping"
            if reason:
                unresolved.append({**_evidence(entry, accession),
                                   "uniprot_accession": accession, "reason": reason})
                continue
            group = grouped.setdefault(accession, {
                "canonical_ref": protein_ref(accession),
                "uniprot_accession": accession,
                "mapping_status": _VERIFIED,
                "inventory_ids": set(), "gene_symbols": set(),
                "ensembl_gene_ids": set(), "evidence_links": [],
            })
            evidence = _evidence(entry, accession)
            group["inventory_ids"].add(ident)
            if evidence["gene_symbol"]:
                group["gene_symbols"].add(evidence["gene_symbol"])
            if evidence["ensembl_gene_id"]:
                group["ensembl_gene_ids"].add(evidence["ensembl_gene_id"])
            group["evidence_links"].append(evidence)
    result = []
    for accession in sorted(grouped):
        group = grouped[accession]
        for key in ("inventory_ids", "gene_symbols", "ensembl_gene_ids"):
            group[key] = sorted(group[key])
        group["evidence_links"] = _unique_dicts(group["evidence_links"])
        result.append(group)
    return result, _unique_dicts(unresolved), gene_ids


def collection_protein_links(entries) -> list[dict]:
    """Group verified gene-to-accession mappings under existing simulator refs.

    Only an exact ``reviewed_snapshot_match`` in the entry's source lookups
    establishes a verified mapping. Unresolved associations are excluded here
    and reported explicitly by :func:`reconcile_inventory`. This verification
    concerns the recorded identifier, not coordinates or a particular proteoform.
    """
    return _collection_links(entries)[0]


def _read_proteome(target: Path) -> tuple[dict, set[str]]:
    path = target / "macromolecules/catalog/human_proteome.sqlite"
    report = {"path": str(path), "status": "missing", "indexed_entries": None}
    if not path.exists():
        return report, set()
    if not path.is_file():
        return {**report, "status": "unreadable", "error": "Database path is not a file"}, set()
    try:
        # mode=ro is essential: sqlite3.connect(path) would create a missing DB.
        with closing(sqlite3.connect(path.as_uri() + "?mode=ro", uri=True)) as conn:
            columns = {row[1] for row in conn.execute("PRAGMA table_info(proteins)")}
            if "accession" not in columns:
                raise ValueError("Human Proteome database has no proteins.accession column")
            accessions = set()
            for (accession,) in conn.execute("SELECT accession FROM proteins"):
                normalized = protein_ref(accession).split("/", 1)[1]
                if normalized != accession or normalized in accessions:
                    raise ValueError("Human Proteome database has noncanonical or repeated accessions")
                accessions.add(normalized)
        return {**report, "status": "indexed", "indexed_entries": len(accessions)}, accessions
    except (OSError, sqlite3.Error, ValueError) as exc:
        return {**report, "status": "unreadable", "error": str(exc)}, set()


def reconcile_inventory(entries, target) -> dict:
    """Report canonical reuse against an existing checkout; no writes/network.

    ``reference_only`` means the target database cannot be inspected. It is not
    evidence that a protein is missing on the user's other computers. A verified
    source accession absent from a readable index is ``missing_from_index``;
    neither status authorizes creating a second molecule definition.
    """
    root = Path(target).expanduser().resolve()
    if not (root / "main.py").is_file() or not (root / "sim/loader.py").is_file():
        raise ValueError("Target must be an existing CellSimulator checkout")
    links, unresolved, gene_ids = _collection_links(entries)
    database, indexed = _read_proteome(root)
    for link in links:
        link["catalog_status"] = (
            "reference_only" if database["status"] != "indexed" else
            "known" if link["uniprot_accession"] in indexed else "missing_from_index"
        )
    counts = {
        "inventory_gene_records": len(gene_ids),
        "unique_verified_protein_refs": len(links),
        "known_protein_refs": sum(link["catalog_status"] == "known" for link in links),
        "reference_only_protein_refs": sum(link["catalog_status"] == "reference_only" for link in links),
        "missing_from_index_protein_refs": sum(link["catalog_status"] == "missing_from_index" for link in links),
        "unique_unresolved_accessions": len({row["uniprot_accession"] for row in unresolved
                                             if row["uniprot_accession"]}),
        "genes_without_accessions": len({row["inventory_id"] for row in unresolved
                                          if row["reason"] == "no_accession_on_inventory_entry"}),
        "unresolved_gene_records": len({row["inventory_id"] for row in unresolved}),
        "collapsed_gene_to_protein_links": sum(len(link["inventory_ids"]) - 1 for link in links),
    }
    return {
        "schema_version": 1, "target": str(root), "proteome_database": database,
        "counts": counts, "protein_links": links, "unresolved_mappings": unresolved,
        "coordinate_availability": "not_inspected",
        "population_assignment": "not_assigned",
        "notes": [
            "The source gene records retain their localization evidence; shared accessions reuse one protein reference.",
            "An accession identifies a protein record, not every isoform, modification or structural state.",
            "A missing or unreadable local database does not establish that the user's proteome is absent elsewhere.",
            "No names, sequences or PDB cross-references are used as substitutes for an exact identity mapping.",
            "No structures, molecules, databases or scene populations were created or changed.",
        ],
    }


def _human_entity_ids(entities: dict) -> set[str]:
    taxa: dict[str, set[str]] = {}
    for category, taxon_key in (
        ("genetically_engineered_source", "pdbx_gene_src_ncbi_taxonomy_id"),
        ("natural_source", "pdbx_ncbi_taxonomy_id"),
    ):
        for source in entities.get(category, []):
            ident = _text(source.get("entity_id"))
            if ident:
                # Use biological origin, never the recombinant expression host.
                taxa.setdefault(ident, set()).add(str(source.get(taxon_key, "")).strip())
    proteins = {_text(row.get("entity_id")) for row in entities.get("polymers", [])
                if row.get("type") == "polypeptide(L)"}
    return {ident for ident, values in taxa.items() if values == {"9606"} and ident in proteins}


def _component_identity(accession, isoform=None) -> dict | None:
    accession = _text(accession).upper()
    if not _WITH_ISOFORM.fullmatch(accession):
        return None
    # Source CIF usually stores text, but a JSON numeric isoform must not be
    # mistaken for no isoform and redirected to the canonical protein.
    annotation = str(isoform).strip().upper() if isoform is not None else ""
    if annotation in {".", "?"}:
        annotation = ""
    status = "source_entity_reference"
    if annotation:
        base = accession.split("-", 1)[0]
        selected = base + "-" + annotation if re.fullmatch(r"[1-9][0-9]*", annotation) else annotation
        if (_WITH_ISOFORM.fullmatch(selected) and "-" in selected
                and selected.split("-", 1)[0] == base
                and ("-" not in accession or accession == selected)):
            accession = selected
        else:
            status = "conflicting_or_unresolved_isoform_annotation"
    canonical = bool(_CANONICAL.fullmatch(accession)) and status == "source_entity_reference"
    if not canonical and status == "source_entity_reference":
        status = "isoform_not_supported_by_simulator"
    return {
        "uniprot_accession": accession,
        "canonical_ref": protein_ref(accession) if canonical else None,
        "isoform_annotation": annotation or None,
        "mapping_status": status,
    }


def asset_component_refs(metadata: dict) -> list[dict]:
    """Link source-identified human components; never equate assembly and protein.

    Returned dictionaries contain a ``canonical_ref`` only when the existing
    simulator can represent the declared accession. Explicit isoforms retain
    their suffix and have a null canonical ref, rather than being collapsed.
    These are provenance links; current-catalog membership, complete sequence
    coverage, modifications and structural equivalence are not established.
    """
    source = metadata.get("source") or {}
    if metadata.get("provenance_class") == "predicted_model":
        organism = metadata.get("organism") or {}
        if str(organism.get("taxonomy_id", "")) != "9606":
            return []
        identity = _component_identity(source.get("uniprot_accession"))
        if identity is None:
            return []
        return [{**identity, "entity_ids": [], "chain_ids": [],
                 "relationship": "predicted_sequence_reference", "taxonomy_id": 9606,
                 "source_model_id": source.get("model_id"),
                 "full_protein_identity": False, "sequence_state_equivalence": False,
                 "catalog_verification": "not_checked"}]
    entities = metadata.get("entities") or {}
    human = _human_entity_ids(entities)
    polymers = {_text(row.get("entity_id")): row for row in entities.get("polymers", [])}
    grouped: dict[tuple, dict] = {}
    for row in entities.get("sequence_references", []):
        entity = _text(row.get("entity_id"))
        if entity not in human or _text(row.get("db_name")).upper() not in {"UNP", "UNIPROT"}:
            continue
        identity = _component_identity(row.get("pdbx_db_accession"), row.get("pdbx_db_isoform"))
        if identity is None:
            continue
        key = (identity["uniprot_accession"], identity["mapping_status"], identity["isoform_annotation"] or "")
        group = grouped.setdefault(key, {
            **identity, "entity_ids": set(), "chain_ids": set(),
            "relationship": "assembly_component_reference", "taxonomy_id": 9606,
            "source_pdb_id": source.get("pdb_id"),
            "assembly_id": (metadata.get("assembly") or {}).get("id"),
            "full_protein_identity": False, "sequence_state_equivalence": False,
            "catalog_verification": "not_checked",
            "reference_scope": "source_entry_entity",
            "assembly_membership": "not_verified_from_metadata",
            "coordinate_presence": "not_checked",
        })
        group["entity_ids"].add(entity)
        group["chain_ids"].update(chain.strip() for chain in
                                  _text(polymers[entity].get("pdbx_strand_id")).split(",") if chain.strip())
    result = []
    for key in sorted(grouped):
        group = grouped[key]
        group["entity_ids"] = sorted(group["entity_ids"])
        group["chain_ids"] = sorted(group["chain_ids"])
        result.append(group)
    return result
