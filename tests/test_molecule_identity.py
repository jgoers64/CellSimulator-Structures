"""Offline identity regressions, prepared for explicit test execution only."""
from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import sqlite3
import tempfile
import unittest

from structure_library.identity import (
    asset_component_refs, collection_protein_links, protein_ref, reconcile_inventory,
)


def gene(ident="protein-ensg00000275714", symbol="H3C1", accession="P68431",
         status="reviewed_snapshot_match"):
    return {
        "id": ident, "kind": "protein_gene", "gene_symbol": symbol,
        "uniprot_accessions": [] if accession is None else [accession],
        "source_lookups": ([] if accession is None else [{
            "uniprot_accession": accession, "mapping_status": status,
        }]),
        "source_record": {
            "ensembl_gene_id": ident.removeprefix("protein-").upper(),
            "source_url": f"https://www.proteinatlas.org/{ident}/subcellular",
            "source_release": "Synthetic test evidence",
            "evidence_scope": "interphase_annotation",
        },
    }


def checkout(folder):
    root = Path(folder) / "CellSimulator"
    (root / "sim").mkdir(parents=True)
    (root / "main.py").write_text("", encoding="utf-8")
    (root / "sim/loader.py").write_text("", encoding="utf-8")
    return root


def database(root, accessions):
    path = root / "macromolecules/catalog/human_proteome.sqlite"
    path.parent.mkdir(parents=True)
    with sqlite3.connect(path) as conn:
        # Older indexes are supported; reconciliation only needs exact IDs.
        conn.execute("CREATE TABLE proteins (accession TEXT PRIMARY KEY)")
        conn.executemany("INSERT INTO proteins VALUES (?)", [(value,) for value in accessions])
    return path


def assembly():
    return {
        "source": {"pdb_id": "9ZZZ"}, "assembly": {"id": "1"},
        "entities": {
            "polymers": [
                {"entity_id": "1", "type": "polypeptide(L)", "pdbx_strand_id": "A,B"},
                {"entity_id": "2", "type": "polypeptide(L)", "pdbx_strand_id": "C"},
            ],
            "genetically_engineered_source": [
                {"entity_id": "1", "pdbx_gene_src_ncbi_taxonomy_id": "9606",
                 "pdbx_host_org_ncbi_taxonomy_id": "562"},
            ],
            "natural_source": [{"entity_id": "2", "pdbx_ncbi_taxonomy_id": "9606"}],
            "sequence_references": [
                {"entity_id": "1", "db_name": "UNP", "pdbx_db_accession": "P68431"},
                {"entity_id": "2", "db_name": "UNP", "pdbx_db_accession": "P68431"},
            ],
        },
    }


class CanonicalProteinTests(unittest.TestCase):
    def test_existing_simulator_reference_and_strict_identifiers(self):
        for accession in ("P68431", "O60814", "Q9Y6K9", "A0A024RBG1"):
            self.assertEqual(protein_ref(" " + accession.lower() + " "),
                             "human-proteome/" + accession)
        for accession in ("P68431-2", "H3C1", "HISTONE", "../P68431", "", None):
            with self.subTest(accession=accession), self.assertRaises(ValueError):
                protein_ref(accession)

    def test_two_genes_share_one_protein_but_preserve_both_evidence_records(self):
        entries = [gene(), gene("protein-ensg00000286522", "H3C2")]
        before = deepcopy(entries)
        links = collection_protein_links(entries + [deepcopy(entries[0])])
        self.assertEqual(len(links), 1)
        self.assertEqual(links[0]["canonical_ref"], "human-proteome/P68431")
        self.assertEqual(links[0]["gene_symbols"], ["H3C1", "H3C2"])
        self.assertEqual(len(links[0]["inventory_ids"]), 2)
        self.assertEqual(len(links[0]["evidence_links"]), 2)
        self.assertEqual(entries, before)
        self.assertEqual(links, collection_protein_links(list(reversed(entries))))

    def test_names_do_not_merge_distinct_accessions(self):
        entries = [gene(), gene("protein-ensg00000000002", "H3C1", "P62805")]
        self.assertEqual(len(collection_protein_links(entries)), 2)

    def test_lookups_must_verify_the_same_accession_without_conflicting_status(self):
        unknown = gene(status="not_found_in_reviewed_snapshot")
        mismatched = gene("protein-ensg00000000002")
        mismatched["source_lookups"][0]["uniprot_accession"] = "P62805"
        conflict = gene("protein-ensg00000000003")
        conflict["source_lookups"].append({"uniprot_accession": "P68431", "mapping_status": "not_queried"})
        self.assertEqual(collection_protein_links([unknown, mismatched, conflict]), [])

    def test_conflicting_inventory_identity_is_not_silently_merged(self):
        with self.assertRaisesRegex(ValueError, "Conflicting identities"):
            collection_protein_links([gene(), gene(accession="P62805")])

    def test_bundled_histone_gene_evidence_reuses_the_existing_protein_namespace(self):
        root = Path(__file__).resolve().parents[1]
        entries = json.loads((root / "collections/human-nucleus/inventory.json").read_text(encoding="utf-8"))["entries"]
        matching = [row for row in collection_protein_links(entries)
                    if row["uniprot_accession"] == "P68431"]
        self.assertEqual(len(matching), 1)
        self.assertTrue({"H3C1", "H3C2"}.issubset(matching[0]["gene_symbols"]))
        self.assertEqual(matching[0]["canonical_ref"], "human-proteome/P68431")


class ReconciliationTests(unittest.TestCase):
    def test_missing_local_database_is_reference_only_and_is_not_created(self):
        with tempfile.TemporaryDirectory() as folder:
            target = checkout(folder)
            before = sorted(path.relative_to(target) for path in target.rglob("*"))
            result = reconcile_inventory([gene()], target)
            self.assertEqual(result["proteome_database"]["status"], "missing")
            self.assertEqual(result["counts"]["reference_only_protein_refs"], 1)
            self.assertEqual(result["counts"]["missing_from_index_protein_refs"], 0)
            self.assertEqual(before, sorted(path.relative_to(target) for path in target.rglob("*")))

    def test_existing_index_distinguishes_known_and_absent_without_mutation(self):
        with tempfile.TemporaryDirectory() as folder:
            target = checkout(folder)
            path = database(target, ["P68431"])
            before = path.read_bytes()
            result = reconcile_inventory([gene(), gene("protein-ensg00000000002", "H4C1", "P62805")], target)
            self.assertEqual(result["proteome_database"]["status"], "indexed")
            self.assertEqual(result["counts"]["known_protein_refs"], 1)
            self.assertEqual(result["counts"]["missing_from_index_protein_refs"], 1)
            self.assertEqual(result["counts"]["reference_only_protein_refs"], 0)
            self.assertEqual(path.read_bytes(), before)

    def test_unreadable_database_does_not_claim_every_protein_is_missing(self):
        with tempfile.TemporaryDirectory() as folder:
            target = checkout(folder)
            path = target / "macromolecules/catalog/human_proteome.sqlite"
            path.parent.mkdir(parents=True)
            path.write_bytes(b"not a sqlite file")
            result = reconcile_inventory([gene()], target)
            self.assertEqual(result["proteome_database"]["status"], "unreadable")
            self.assertEqual(result["counts"]["reference_only_protein_refs"], 1)
            self.assertEqual(result["counts"]["missing_from_index_protein_refs"], 0)
            self.assertEqual(path.read_bytes(), b"not a sqlite file")

    def test_unresolved_gene_and_isoform_rows_are_retained_without_false_refs(self):
        with tempfile.TemporaryDirectory() as folder:
            target = checkout(folder)
            entries = [gene(accession=None),
                       gene("protein-ensg00000000002", accession="P68431-2"),
                       gene("protein-ensg00000000003", accession="Q6ZMK1", status="not_found_in_reviewed_snapshot")]
            result = reconcile_inventory(entries, target)
            self.assertEqual(result["protein_links"], [])
            self.assertEqual(result["counts"]["unresolved_gene_records"], 3)
            self.assertEqual(result["counts"]["unique_unresolved_accessions"], 2)
            self.assertEqual(result["counts"]["genes_without_accessions"], 1)
            self.assertIn("P68431-2", {row["uniprot_accession"] for row in result["unresolved_mappings"]})


class AssemblyReferenceTests(unittest.TestCase):
    def test_human_components_group_and_expression_host_does_not_change_origin(self):
        links = asset_component_refs(assembly())
        self.assertEqual(len(links), 1)
        self.assertEqual(links[0]["canonical_ref"], "human-proteome/P68431")
        self.assertEqual(links[0]["entity_ids"], ["1", "2"])
        self.assertEqual(links[0]["chain_ids"], ["A", "B", "C"])
        self.assertEqual(links[0]["relationship"], "assembly_component_reference")
        self.assertFalse(links[0]["full_protein_identity"])
        self.assertFalse(links[0]["sequence_state_equivalence"])

    def test_nonhuman_unknown_and_mixed_origin_entities_are_not_human_refs(self):
        for taxon in (None, "10090"):
            metadata = assembly()
            metadata["entities"]["natural_source"][0]["pdbx_ncbi_taxonomy_id"] = taxon
            metadata["entities"]["genetically_engineered_source"][0]["pdbx_gene_src_ncbi_taxonomy_id"] = taxon
            self.assertEqual(asset_component_refs(metadata), [])
        metadata = assembly()
        metadata["entities"]["natural_source"].append({"entity_id": "1", "pdbx_ncbi_taxonomy_id": "10090"})
        self.assertEqual(asset_component_refs(metadata)[0]["entity_ids"], ["2"])

    def test_isoform_accession_or_annotation_is_preserved_without_canonical_redirect(self):
        for full_accession, annotation in (("P68431-2", None), ("P68431", "P68431-2"),
                                           ("P68431", "2"), ("P68431", 2)):
            metadata = assembly()
            metadata["entities"]["sequence_references"] = [{
                "entity_id": "1", "db_name": "UNP", "pdbx_db_accession": full_accession,
                "pdbx_db_isoform": annotation,
            }]
            link = asset_component_refs(metadata)[0]
            self.assertEqual(link["uniprot_accession"], "P68431-2")
            self.assertIsNone(link["canonical_ref"])
            self.assertEqual(link["mapping_status"], "isoform_not_supported_by_simulator")

    def test_conflicting_isoform_annotation_blocks_canonical_link(self):
        metadata = assembly()
        metadata["entities"]["sequence_references"] = [{
            "entity_id": "1", "db_name": "UNP", "pdbx_db_accession": "P68431",
            "pdbx_db_isoform": "P62805-2",
        }]
        link = asset_component_refs(metadata)[0]
        self.assertIsNone(link["canonical_ref"])
        self.assertEqual(link["mapping_status"], "conflicting_or_unresolved_isoform_annotation")

    def test_human_prediction_links_record_but_never_claims_sequence_state_equivalence(self):
        metadata = {"provenance_class": "predicted_model", "organism": {"taxonomy_id": 9606},
                    "source": {"uniprot_accession": "P68431", "model_id": "AF-P68431-F1"}}
        link = asset_component_refs(metadata)[0]
        self.assertEqual(link["canonical_ref"], "human-proteome/P68431")
        self.assertEqual(link["relationship"], "predicted_sequence_reference")
        self.assertFalse(link["sequence_state_equivalence"])
        metadata["source"]["uniprot_accession"] = "P68431-2"
        self.assertIsNone(asset_component_refs(metadata)[0]["canonical_ref"])
        metadata["organism"]["taxonomy_id"] = None
        self.assertEqual(asset_component_refs(metadata), [])


if __name__ == "__main__":
    unittest.main()
