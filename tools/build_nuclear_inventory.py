"""Build the nuclear identity registry without assigning unsupported coordinates.

The source is the reviewed BioSim research inventory JSON. Rebuilding is offline
unless --fetch-uniprot-crossrefs is explicitly selected. UniProt cross-references
are discovery evidence only: they neither select a full protein/complex model nor
verify the availability, biological state, coverage, or quality of coordinates.

Example (from the structures repository):
    python tools/build_nuclear_inventory.py --source inventory-research.json \
        --output . --uniprot-tsv reviewed-human-crossrefs.tsv
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import re
from collections import Counter
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen


COLLECTION_ID = "human-nucleus"
SOURCE_LIMIT = 50_000_000
UNIPROT_LIMIT = 10_000_000
UNIPROT_FIELDS = "accession,id,gene_names,length,xref_pdb,xref_alphafolddb"
UNIPROT_QUERY = "(reviewed:true) AND (organism_id:9606)"
UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/stream?" + urlencode({
    "query": UNIPROT_QUERY, "format": "tsv", "fields": UNIPROT_FIELDS,
})
UNIPROT_ACCESSION = re.compile(
    r"(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-[0-9]+)?"
)
RNA_CATEGORIES = {
    "Ribosomal RNA", "Spliceosomal small nuclear RNA", "RNA processing guides",
    "Transcript and processing families", "Regulatory noncoding RNA",
}


def digest(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def load_json(path: Path):
    raw = path.read_bytes()
    if len(raw) > SOURCE_LIMIT:
        raise ValueError(f"Source JSON exceeds the {SOURCE_LIMIT}-byte budget")
    return json.loads(raw), raw


def normalized_accessions(values):
    """HPA can put comma-separated accessions inside one exported field."""
    result = []
    for field in values:
        for value in re.split(r"[,;]", field):
            value = value.strip()
            if value and not UNIPROT_ACCESSION.fullmatch(value):
                raise ValueError(f"Unrecognized UniProt accession: {value!r}")
            if value and value not in result:
                result.append(value)
    return result


def slug(name):
    value = re.sub(r"[^a-z0-9]+", "-", name.lower()).strip("-")
    if not value:
        raise ValueError("Cannot create an identity from an empty molecule name")
    return value


def fetch_uniprot():
    request = Request(UNIPROT_URL, headers={
        "User-Agent": "CellSimulator-Structures/1.0 (nuclear cross-reference catalog)",
        "Accept": "text/plain; format=tsv",
    })
    with urlopen(request, timeout=60) as response:
        raw = response.read(UNIPROT_LIMIT + 1)
        provenance = {
            "url": UNIPROT_URL,
            "sha256": digest(raw),
            "bytes": len(raw),
            "release": response.headers.get("x-uniprot-release"),
            "release_date": response.headers.get("x-uniprot-release-date"),
            "retrieved_server_date": response.headers.get("date"),
            "scope": "Reviewed human UniProtKB entries; cross-references only",
        }
    if len(raw) > UNIPROT_LIMIT:
        raise ValueError("UniProt cross-reference source exceeded its byte budget")
    return raw, provenance


def read_uniprot(raw, provenance):
    if len(raw) > UNIPROT_LIMIT:
        raise ValueError("UniProt cross-reference source exceeded its byte budget")
    if provenance.get("sha256") and provenance["sha256"] != digest(raw):
        raise ValueError("UniProt source SHA256 differs from its provenance record")
    reader = csv.DictReader(io.StringIO(raw.decode("utf-8")), delimiter="\t")
    required = {"Entry", "Entry Name", "Gene Names", "Length", "PDB", "AlphaFoldDB"}
    if not required.issubset(reader.fieldnames or []):
        raise ValueError("UniProt TSV is missing expected cross-reference columns")
    rows = {}
    for row in reader:
        accession = row["Entry"]
        if not UNIPROT_ACCESSION.fullmatch(accession) or accession in rows:
            raise ValueError(f"Invalid or repeated UniProt primary accession: {accession!r}")
        rows[accession] = row
    provenance = {
        **provenance, "sha256": digest(raw), "bytes": len(raw),
        "rows": len(rows), "row_selection": UNIPROT_QUERY,
        "attribution": "UniProt Consortium; UniProtKB structural cross-references",
        "license": "CC BY 4.0",
        "license_url": "https://www.uniprot.org/help/license",
        "interpretation": "A recorded reference is not a downloaded, selected, complete or scientifically reviewed structure.",
    }
    return rows, provenance


def lookup(accession, crossrefs, have_snapshot):
    row = crossrefs.get(accession)
    refs = {
        "uniprot_accession": accession,
        "uniprot_url": f"https://www.uniprot.org/uniprotkb/{accession}/entry",
        "alphafold_metadata_url": f"https://alphafold.ebi.ac.uk/api/prediction/{accession}",
        "mapping_status": (
            "reviewed_snapshot_match" if row else
            "not_found_in_reviewed_snapshot" if have_snapshot else "not_queried"
        ),
        "pdb_ids": [], "alphafold_accessions": [],
        "sequence_length": None,
    }
    if row:
        refs.update({
            "pdb_ids": [x.strip() for x in row["PDB"].split(";") if x.strip()],
            "alphafold_accessions": [x.strip() for x in row["AlphaFoldDB"].split(";") if x.strip()],
            "sequence_length": int(row["Length"]),
        })
    return refs


def protein_entry(row, crossrefs, have_snapshot, legacy_ids):
    ensembl = row["ensembl_gene_id"]
    if not re.fullmatch(r"ENSG[0-9]+", ensembl):
        raise ValueError(f"Invalid Ensembl gene identity: {ensembl!r}")
    accessions = normalized_accessions(row["uniprot_accessions"])
    representative_links = []
    if "ubiquitin" in legacy_ids and row["gene_symbol"] in {"UBB", "UBC", "UBA52", "RPS27A"}:
        representative_links.append({
            "catalog_entry_id": "ubiquitin", "pdb_id": "1UBQ",
            "relationship": "representative_mature_ubiquitin_component_only",
            "full_gene_product_model": False,
            "note": "Free mature ubiquitin is one processed component of the gene product; 1UBQ is not the full fusion/polyubiquitin precursor. The HPA gene-level localization does not independently localize each processed component.",
        })
    return {
        "id": "protein-" + ensembl.lower(),
        "name": row["protein_description"],
        "kind": "protein_gene",
        "gene_symbol": row["gene_symbol"],
        "uniprot_accessions": accessions,
        "copies_per_default_g1_nucleus": None,
        "abundance_rank_default_g1": None,
        "source_lookups": [lookup(a, crossrefs, have_snapshot) for a in accessions],
        "existing_representative_links": representative_links,
        "source_record": row,
    }


def molecule_kind(row):
    if row["category"] in RNA_CATEGORIES:
        return "rna_family_or_transcript"
    if row["category"] == "Nuclear envelope lipids":
        return "lipid_class"
    if row["name"] == "Genomic DNA":
        return "genome"
    if row["name"] == "Free proteinogenic amino acids":
        return "chemical_family"
    if "ionic" in row["category"].lower() or row["category"] == "Structural metal":
        return "ionic_pool"
    return "chemical_species_candidate"


def molecule_entry(row):
    return {
        "id": "molecule-" + slug(row["name"]),
        "name": row["name"], "kind": molecule_kind(row),
        "copies_per_default_g1_nucleus": None,
        "abundance_rank_default_g1": None,
        "source_record": row,
    }


def write_json(path, data, record_key=None):
    """One JSON record per line keeps large registries compact and diffable."""
    path.parent.mkdir(parents=True, exist_ok=True)
    if record_key is None:
        text = json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False) + "\n"
    else:
        head = {k: v for k, v in data.items() if k != record_key}
        prefix = json.dumps(head, ensure_ascii=False, indent=2, allow_nan=False).rstrip()
        records = ",\n".join(
            "    " + json.dumps(row, ensure_ascii=False, separators=(",", ":"), allow_nan=False)
            for row in data[record_key]
        )
        text = prefix[:-1].rstrip() + ",\n  " + json.dumps(record_key) + ": [\n" + records + "\n  ]\n}\n"
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(text, encoding="utf-8")
    temp.replace(path)


def build(source, output, uniprot_raw=None, uniprot_provenance=None):
    data, raw = load_json(source)
    source_digest = digest(raw)
    crossrefs, crossref_provenance = {}, None
    if uniprot_raw is not None:
        crossrefs, crossref_provenance = read_uniprot(uniprot_raw, uniprot_provenance or {})
    legacy_path = output / "catalog.json"
    legacy_ids = {r["id"] for r in json.loads(legacy_path.read_text())["entries"]} if legacy_path.exists() else set()
    rows = data["hpa_protein_catalog"]
    entries = [protein_entry(r, crossrefs, uniprot_raw is not None, legacy_ids) for r in rows]
    entries += [molecule_entry(r) for r in data["other_molecules_and_families"]]
    if len({r["id"] for r in entries}) != len(entries):
        raise ValueError("Repeated stable identity; source needs explicit disambiguation")
    for row in rows + data["other_molecules_and_families"]:
        if row.get("copies_per_default_g1_nucleus") is not None or row.get("abundance_rank_default_g1") is not None:
            raise ValueError("This builder cannot transfer population/rank estimates into the default G1 model")
    protein_entries = [r for r in entries if r["kind"] == "protein_gene"]
    requested_accessions = {a for r in protein_entries for a in r["uniprot_accessions"]}
    counts = {
        "identity_entries": len(entries), "protein_genes": len(rows),
        "other_molecules_and_families": len(data["other_molecules_and_families"]),
        "structural_candidates": len(data["structural_candidates"]),
        "hpa_evidence_scopes": dict(Counter(r["evidence_scope"] for r in rows)),
        "unique_normalized_uniprot_accessions": len(requested_accessions),
        "genes_without_uniprot_accessions": sum(not r["uniprot_accessions"] for r in protein_entries),
        "normalized_comma_or_semicolon_fields": sum(
            any(re.search(r"[,;]", a) for a in r["uniprot_accessions"]) for r in rows
        ),
        "genes_with_reviewed_uniprot_mapping": sum(
            any(x["mapping_status"] == "reviewed_snapshot_match" for x in r["source_lookups"])
            for r in protein_entries
        ),
        "genes_with_recorded_pdb_reference": sum(
            any(x["pdb_ids"] for x in r["source_lookups"]) for r in protein_entries
        ),
        "genes_with_recorded_alphafold_reference": sum(
            any(x["alphafold_accessions"] for x in r["source_lookups"]) for r in protein_entries
        ),
        "unmapped_accessions_in_reviewed_snapshot": sorted(requested_accessions - crossrefs.keys()) if uniprot_raw is not None else None,
    }
    collection = {
        "schema_version": 1, "collection_id": COLLECTION_ID,
        "name": "Human nucleus molecular identity registry",
        "scope": "Research inventory across cell contexts, not the composition of a single measured G1 nucleus.",
        "complete_nuclear_inventory": False,
        "default_profile": data["metadata"]["default_biosim_profile"],
        "counts": counts,
        "source_inventory": {
            "filename": source.name, "sha256": source_digest,
            "metadata": data["metadata"],
            "record_preservation": "Original protein/molecule source_record fields are unchanged; normalized identifiers are separate.",
        },
        "uniprot_cross_reference_snapshot": crossref_provenance,
        "rules": [
            "Identity registration and a structure cross-reference do not mean a coordinate asset exists locally.",
            "Current coordinate availability is joined from the separate asset manifests, never inferred from this identity record.",
            "A gene entry can encode several isoforms/proteoforms; never equate it automatically with one full-length structure.",
            "PDB cross-references may describe fragments, mutants, complexes or different states; they require selection and review.",
            "AlphaFold references are predicted models; resolve actual metadata/download URLs and confidence before use.",
            "Missing reviewed-UniProt mapping is not evidence of protein or structure absence.",
            "Unknown counts are null; the separate HeLa reference does not populate the default G1 model.",
            "RNA and lipid family records require concrete sequence/species identities before coordinate selection.",
            "No atoms, missing residues, protonation states, chromosome folds, or population counts are generated here.",
        ],
        "source_url_templates": {
            "pdb_entry": "https://www.rcsb.org/structure/{pdb_id}",
            "alphafold_api_documentation": "https://alphafold.com/api-docs",
            "alphafold_version_guidance": "https://www.ebi.ac.uk/pdbe/news/alphafold-database-release-notes",
            "uniprot_api_guidance": "https://www.uniprot.org/help/api_queries",
        },
        "structural_candidates": data["structural_candidates"],
        "structural_review": data["structural_review"],
        "quantitative_anchors": data["quantitative_anchors"],
        "nolin_source_provenance": data["nolin_source_provenance"],
        "entries": entries,
    }
    directory = output / "collections" / COLLECTION_ID
    write_json(directory / "inventory.json", collection, "entries")
    write_json(directory / "abundance-reference.json", {
        "schema_version": 1, "collection_id": COLLECTION_ID,
        "name": "Separate HeLaM derived nuclear-pool protein-group estimates",
        "not_default_g1_population": True,
        "source_inventory_sha256": source_digest,
        "provenance": data["hela_abundance_provenance"],
        "group_count": len(data["hela_ranked_protein_groups"]),
        "groups": data["hela_ranked_protein_groups"],
    }, "groups")
    if uniprot_raw is not None:
        source_dir = directory / "sources"
        source_dir.mkdir(parents=True, exist_ok=True)
        (source_dir / "uniprot-reviewed-human-crossrefs.tsv").write_bytes(uniprot_raw)
        write_json(source_dir / "uniprot-reviewed-human-crossrefs.provenance.json", crossref_provenance)
    return counts


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True, help="Reviewed research inventory JSON")
    parser.add_argument("--output", type=Path, required=True, help="Structures repository root")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--uniprot-tsv", type=Path, help="Replay a saved reviewed-human cross-reference TSV")
    group.add_argument("--fetch-uniprot-crossrefs", action="store_true", help="Explicitly retrieve current UniProt cross-references")
    args = parser.parse_args(argv)
    raw, provenance = None, None
    if args.fetch_uniprot_crossrefs:
        raw, provenance = fetch_uniprot()
    elif args.uniprot_tsv:
        raw = args.uniprot_tsv.read_bytes()
        metadata_path = args.uniprot_tsv.with_suffix(".provenance.json")
        provenance = json.loads(metadata_path.read_text()) if metadata_path.exists() else {
            "url": UNIPROT_URL,
            "release": None,
            "scope": "User-supplied cross-reference snapshot; no HTTP release header supplied",
        }
    counts = build(args.source, args.output, raw, provenance)
    print(json.dumps(counts, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
