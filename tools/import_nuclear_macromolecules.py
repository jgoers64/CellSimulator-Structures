"""Pin curated nuclear mmCIF sources and archive quality metadata.

Explicit maintenance command, not a startup hook or test runner. Source atom files
stay in ignored cache; reviewable metadata and checksums are committed to Git.
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from structure_library.coordinates import (atomic_write, bounded_decompress, category_rows,
    digest, download, require_gemmi, write_json, assembly_atom_upper_bound)


def candidates(source):
    result = json.loads(Path(source).read_text(encoding="utf-8"))["structural_candidates"]
    return [{
        "pdb_id": "2CV5", "name": "Human nucleosome core particle (existing BioSim reference)",
        "assembly": {"preferred_id": "1", "description": "Histone octamer and experimental DNA duplex"},
        "provenance_class": "experimental_coordinate_model",
        "limitations": ["Experimental DNA construct and resolved histone regions only.",
                        "The existing application already uses this template; registry entry prevents a separate untracked copy.",
                        "This structure does not establish a chromosome conformation or nucleosome population."],
        "sources": [{"kind": "PDB", "url": "https://www.rcsb.org/structure/2CV5"}],
    }] + result


def ingest(candidate, root):
    gemmi = require_gemmi()
    pdb_id = candidate["pdb_id"].upper()
    asset_id = "pdb-" + pdb_id.lower()
    folder = Path(root) / "structures" / "nuclear" / asset_id
    cache = Path(root) / "cache" / "nuclear"
    cache.mkdir(parents=True, exist_ok=True)
    url = f"https://files.rcsb.org/download/{pdb_id}.cif.gz"
    raw_path = cache / (asset_id + ".cif")
    if raw_path.exists():
        raw = raw_path.read_bytes()
    else:
        raw = bounded_decompress(download(url))
        atomic_write(raw_path, raw)
    block = gemmi.cif.read_string(raw.decode()).sole_block()
    if block.find_value("_entry.id").upper() != pdb_id:
        raise ValueError(f"Unexpected coordinate identity for {pdb_id}")
    structure = gemmi.make_structure_from_block(block)
    selected = candidate["assembly"]["preferred_id"]
    assembly = next((a for a in structure.assemblies if a.name == selected), None)
    if assembly is None:
        raise ValueError(f"Missing assembly {selected} in {pdb_id}")
    entry_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    entry_path = cache / (asset_id + "-entry.json")
    if entry_path.exists():
        entry_raw = entry_path.read_bytes()
    else:
        entry_raw = download(entry_url, max_bytes=4_000_000)
        atomic_write(entry_path, entry_raw)
    entry = json.loads(entry_raw)
    if entry["rcsb_id"].upper() != pdb_id:
        raise ValueError("Wrong entry metadata returned")
    report_keys = [key for key in entry if key.startswith("pdbx_vrpt_")]
    full_report = f"https://files.rcsb.org/pub/pdb/validation_reports/{pdb_id[1:3].lower()}/{pdb_id.lower()}/{pdb_id.lower()}_full_validation.pdf"
    coverage = entry.get("rcsb_entry_info", {})
    # Complete missing-region records can be large; keep them as a separate,
    # compressed source-derived artifact instead of copying them into each view.
    missing_residues = category_rows(block, "_pdbx_unobs_or_zero_occ_residues.")
    missing_atoms = category_rows(block, "_pdbx_unobs_or_zero_occ_atoms.")
    coverage_path = folder / "coverage.json.gz"
    write_json(coverage_path, {
        "pdb_id": pdb_id, "coordinate_sha256": digest(raw),
        "unobserved_or_zero_occupancy_residue_records": missing_residues,
        "unobserved_or_zero_occupancy_atom_records": missing_atoms,
        "note": "Archive construct records, not a claim of canonical full-length coverage.",
    }, compact=True)
    version = entry.get("rcsb_accession_info", {})
    provenance = candidate["provenance_class"]
    metadata = {
        "schema_version": 1, "id": asset_id, "name": candidate["name"],
        "provenance_class": provenance, "coordinate_status": "pinned_source_available",
        "source": {
            "database": "wwPDB / RCSB PDB", "pdb_id": pdb_id,
            "entry_url": f"https://www.rcsb.org/structure/{pdb_id}",
            "coordinate_url": url, "coordinate_format": "PDBx/mmCIF",
            "coordinate_sha256": digest(raw), "checksum_scope": "uncompressed original mmCIF bytes",
            "coordinate_bytes": len(raw), "coordinate_units": "angstrom",
            "retrieved_at_utc": datetime.fromtimestamp(raw_path.stat().st_mtime, timezone.utc).isoformat(timespec="seconds"),
            "entry_version": version, "primary_citation": entry.get("rcsb_primary_citation"),
            "metadata_url": entry_url, "metadata_sha256": digest(entry_raw),
        },
        "assembly": {"id": selected, "description": candidate["assembly"]["description"],
                     "author_determined": assembly.author_determined,
                     "software_determined": assembly.software_determined,
                     "oligomeric_details": assembly.oligomeric_details,
                     "atom_count_upper_bound": assembly_atom_upper_bound(structure[0], assembly),
                     "coordinate_scope": "Original asymmetric unit plus explicitly selected archive assembly operators"},
        "entities": {
            "descriptions": category_rows(block, "_entity."),
            "polymers": category_rows(block, "_entity_poly."),
            "genetically_engineered_source": category_rows(block, "_entity_src_gen."),
            "natural_source": category_rows(block, "_entity_src_nat."),
            "synthetic_source": category_rows(block, "_pdbx_entity_src_syn."),
            "sequence_references": category_rows(block, "_struct_ref."),
        },
        "coverage": {
            "deposited_atom_count": coverage.get("deposited_atom_count"),
            "modeled_polymer_residues": coverage.get("deposited_modeled_polymer_monomer_count"),
            "deposited_polymer_residues": coverage.get("deposited_polymer_monomer_count"),
            "unmodeled_polymer_residues": coverage.get("deposited_unmodeled_polymer_monomer_count"),
            "deposited_hydrogen_count": coverage.get("deposited_hydrogen_atom_count"),
            "unobserved_residue_record_count": len(missing_residues),
            "unobserved_atom_record_count": len(missing_atoms),
            "full_records_path": coverage_path.relative_to(root).as_posix(),
            "full_records_sha256": digest(coverage_path.read_bytes()),
            "note": "Counts describe deposited constructs; they do not establish full canonical-sequence coverage or nuclear abundance.",
        },
        "quality": {
            "methods": entry.get("exptl"), "resolution_angstrom": coverage.get("resolution_combined"),
            "refinement": entry.get("refine"),
            "archive_validation_summary": {key: entry[key] for key in report_keys},
            "full_validation_report_url": full_report,
            "full_report_review_status": "not_reviewed",
            "local_density_or_map_review_status": "not_reviewed",
            "scientific_accuracy_percent": None,
        },
        "review": {"status": "identity_and_source_metadata_reviewed_with_limitations",
                   "scope": "Pinned archive identity, assembly presence, source checksum, coverage and summary metrics; no independent experimental re-refinement or local-map review."},
        "nuclear_context": {"copies_per_default_g1_nucleus": None, "placement": "not_assigned",
                            "evidence_note": "A structure entry establishes a model, not its quantity or position in the default cell."},
        "rendering": {"view_only": True, "simulation_ready": False,
                      "automatic_atom_budget": 250000,
                      "large_assembly_requires_chunking": assembly_atom_upper_bound(structure[0], assembly) > 250000},
        "limitations": candidate["limitations"], "research_reference": candidate,
    }
    write_json(folder / "metadata.json", metadata)
    atomic_write(folder / "archive-entry.json", entry_raw)
    return {"id": asset_id, "name": metadata["name"],
            "metadata_path": (folder / "metadata.json").relative_to(root).as_posix(),
            "kind": "macromolecule", "provenance_class": provenance,
            "source_status": "pinned_source_available",
            "atom_count_upper_bound": metadata["assembly"]["atom_count_upper_bound"]}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=ROOT / "collections/human-nucleus/inventory.json")
    parser.add_argument("--output", type=Path, default=ROOT)
    args = parser.parse_args(argv)
    root = args.output.resolve()
    manifest = []
    for candidate in candidates(args.source):
        record = ingest(candidate, root)
        manifest.append(record)
        print(f"Pinned {record['id']}: up to {record['atom_count_upper_bound']:,} assembly atoms", flush=True)
    write_json(root / "collections/human-nucleus/macromolecular-assets.json", {
        "schema_version": 1, "collection_id": "human-nucleus",
        "raw_coordinate_policy": "Original mmCIF fetched lazily into ignored cache; pinned checksum checked before every conversion.",
        "assets": manifest,
    })


if __name__ == "__main__":
    main()
