"""Offline installation safeguards, prepared for explicit execution only.

Fixtures are synthetic metadata/geometry used to exercise identity decisions.
Their accession, CCD and source labels do not validate a biological structure.
No test downloads data, launches the simulator, or needs optional dependencies.
"""
from __future__ import annotations

from copy import deepcopy
import hashlib
import json
from pathlib import Path
import sqlite3
import subprocess
import tempfile
import unittest
from unittest import mock
from zipfile import ZipFile

from structure_library.installation import (
    TargetIndex, install_payload, plan_payload, representation_digest,
)


def chemical(*, name="Synthetic component", ident="chemical-ccd-tst", ccd="TST"):
    return {
        "id": ident, "name": name, "pdb_ccd_id": ccd,
        "coordinate_type": "idealized_chemical_component", "coordinate_units": "angstrom",
        "view_only": True, "simulation_ready": False,
        "source_sha256": "a" * 64, "formal_charge": 0, "total_charge": 0.0,
        "atoms": [{"element": "C", "formal_charge": 0}, {"element": "O", "formal_charge": 0}],
        "coords": [[0.0, 0.0, 0.0], [1.25, 0.0, 0.0]],
        "bonds": [[0, 1]], "bond_orders": [2],
    }


def pdb_asset(pdb="9ZZZ"):
    return {
        "asset_id": "pdb-" + pdb.lower(), "name": "Synthetic assembly",
        "source_pdb_id": pdb, "source_sha256": "b" * 64,
        "coordinate_type": "experimental_mmcif", "coordinate_units": "angstrom",
        "view_only": True, "simulation_ready": False,
        "total_charge": 0.0,
        "processing": {"assembly_id": "1", "model_selection": "first deposited model",
                       "atom_selection": "synthetic selection A"},
        "structure_provenance": {
            "source": {"pdb_id": pdb}, "assembly": {"id": "1"},
            "provenance_class": "experimental_coordinate_model",
        },
        "atoms": [{"element": "C", "atom_name": "CA", "chain_id": "A",
                   "residue_name": "ALA", "residue_index": 0}],
        "coords": [[1.0, 2.0, 3.0]], "bonds": [], "bond_orders": [],
    }


def prediction(accession="P68431", taxonomy=9606):
    payload = pdb_asset()
    payload["asset_id"] = "alphafold-synthetic-v1"
    payload.pop("source_pdb_id")
    payload["coordinate_type"] = "predicted_alphafold"
    payload["structure_provenance"] = {
        "provenance_class": "predicted_model",
        "source": {"uniprot_accession": accession, "model_id": "AF-SYNTHETIC-F1"},
        "organism": {"taxonomy_id": taxonomy},
    }
    return payload


def water(*, nist=False):
    payload = chemical(name="Water", ident="chemical-ccd-hoh", ccd="HOH")
    payload.update({"atoms": [{"element": "O"}, {"element": "H"}, {"element": "H"}],
                    "coords": [[0.0, 0.0, 0.0], [0.7, 0.6, 0.0], [-0.7, 0.6, 0.0]],
                    "bonds": [[0, 1], [0, 2]], "bond_orders": [1, 1]})
    if nist:
        payload.pop("id")
        payload.pop("pdb_ccd_id")
        payload.pop("formal_charge")
        payload["provenance"] = {
            "database": "NIST CCCBDB",
            "source_url": "https://cccbdb.nist.gov/expgeom2x.asp?casno=7732185",
        }
    return payload


class InstallationTests(unittest.TestCase):
    def setUp(self):
        self.folder = tempfile.TemporaryDirectory()
        self.addCleanup(self.folder.cleanup)
        self.target = Path(self.folder.name) / "CellSimulator"
        (self.target / "sim").mkdir(parents=True)
        (self.target / "main.py").write_text("", encoding="utf-8")
        (self.target / "sim/loader.py").write_text("", encoding="utf-8")

    def write_json(self, relative, payload):
        path = self.target / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(payload), encoding="utf-8")
        return path

    def create_proteome(self, accessions=("P68431",)):
        path = self.target / "macromolecules/catalog/human_proteome.sqlite"
        path.parent.mkdir(parents=True, exist_ok=True)
        with sqlite3.connect(path) as conn:
            conn.execute("CREATE TABLE proteins (accession TEXT PRIMARY KEY)")
            conn.executemany("INSERT INTO proteins VALUES (?)", [(a,) for a in accessions])
        return path

    def assert_no_new_nuclear_json(self):
        self.assertEqual(list((self.target / "molecules/nucleus").glob("*.json")), [])

    def test_canonical_prediction_reuses_existing_human_proteome_identity(self):
        self.create_proteome()
        payload = prediction()
        plan = plan_payload(payload, self.target)
        self.assertEqual(plan["action"], "reuse_protein_identity")
        self.assertEqual(plan["canonical_ref"], "human-proteome/P68431")
        self.assertFalse(plan["exact_representation_match"])
        self.assertEqual(install_payload(payload, self.target), "human-proteome/P68431")
        self.assert_no_new_nuclear_json()

    def test_normalized_prediction_accession_matches_existing_index(self):
        self.create_proteome()
        plan = plan_payload(prediction(" p68431 "), self.target)
        self.assertEqual(plan["action"], "reuse_protein_identity")
        self.assertEqual(plan["canonical_ref"], "human-proteome/P68431")

    def test_missing_proteome_is_not_created_and_does_not_install_duplicate(self):
        payload = prediction()
        plan = plan_payload(payload, self.target)
        self.assertEqual(plan["action"], "needs_proteome_index")
        self.assertEqual(plan["catalog_status"], "database_missing")
        with self.assertRaises(ValueError):
            install_payload(payload, self.target)
        self.assertFalse((self.target / "macromolecules/catalog/human_proteome.sqlite").exists())
        self.assert_no_new_nuclear_json()

    def test_isoform_is_not_stripped_to_a_canonical_record(self):
        self.create_proteome()
        payload = prediction("P68431-2")
        plan = plan_payload(payload, self.target)
        self.assertEqual(plan["action"], "needs_review")
        self.assertNotEqual(plan.get("canonical_ref"), "human-proteome/P68431")
        with self.assertRaises(ValueError):
            install_payload(payload, self.target)
        self.assert_no_new_nuclear_json()

    def test_prediction_requires_explicit_human_taxonomy(self):
        self.create_proteome()
        for taxonomy in (None, 10090):
            with self.subTest(taxonomy=taxonomy):
                payload = prediction(taxonomy=taxonomy)
                self.assertEqual(plan_payload(payload, self.target)["action"], "needs_review")
                with self.assertRaises(ValueError):
                    install_payload(payload, self.target)
        self.assert_no_new_nuclear_json()

    def test_existing_nist_water_is_reused_without_rewriting_geometry(self):
        path = self.write_json("molecules/simple_molecules/triatomic/water.json", water(nist=True))
        original = path.read_bytes()
        plan = plan_payload(water(), self.target)
        self.assertEqual(plan["action"], "reuse_chemical_identity")
        self.assertFalse(plan["exact_representation_match"])
        self.assertEqual(install_payload(water(), self.target), path)
        self.assertEqual(path.read_bytes(), original)
        self.assert_no_new_nuclear_json()

    def test_name_alone_does_not_authorize_merging_or_adding_a_copy(self):
        existing = chemical()
        existing.pop("id")
        existing.pop("pdb_ccd_id")
        path = self.write_json("molecules/user/example.json", existing)
        original = path.read_bytes()
        plan = plan_payload(chemical(), self.target)
        self.assertEqual(plan["action"], "needs_review")
        with self.assertRaises(ValueError):
            install_payload(chemical(), self.target)
        self.assertEqual(path.read_bytes(), original)
        self.assert_no_new_nuclear_json()

    def test_same_ccd_with_different_charge_requires_review(self):
        existing = chemical(ident="older-local-id")
        path = self.write_json("molecules/user/state.json", existing)
        changed = chemical()
        changed["formal_charge"] = 1
        changed["total_charge"] = 1.0
        changed["atoms"][0]["formal_charge"] = 1
        self.assertEqual(plan_payload(changed, self.target)["action"], "needs_review")
        with self.assertRaises(ValueError):
            install_payload(changed, self.target)
        self.assertEqual(json.loads(path.read_text())["formal_charge"], 0)

    def test_charged_incoming_hoh_is_not_aliased_to_neutral_nist_water(self):
        self.write_json("molecules/simple_molecules/triatomic/water.json", water(nist=True))
        changed = water()
        changed["formal_charge"] = 1
        changed["total_charge"] = 1.0
        changed["atoms"][0]["formal_charge"] = 1
        self.assertEqual(plan_payload(changed, self.target)["action"], "needs_review")
        with self.assertRaises(ValueError):
            install_payload(changed, self.target)

    def test_same_pdb_with_different_selection_is_not_exact_reuse(self):
        old = pdb_asset()
        old["asset_id"] = "older-9zzz-representation"
        path = self.write_json("molecules/user/assembly.json", old)
        changed = pdb_asset()
        changed["processing"]["atom_selection"] = "synthetic selection B"
        changed["processing"]["assembly_id"] = "2"
        changed["structure_provenance"]["assembly"]["id"] = "2"
        plan = plan_payload(changed, self.target)
        self.assertEqual(plan["action"], "needs_review")
        self.assertFalse(plan["exact_representation_match"])
        with self.assertRaises(ValueError):
            install_payload(changed, self.target)
        self.assertEqual(json.loads(path.read_text())["processing"]["assembly_id"], "1")

    def builtin(self, *, checksum=True, changed=False):
        raw = b"synthetic shared-geometry fixture; no biological model claimed"
        asset = self.target / "macromolecules/assets/nucleosome-2cv5.atoms.json.gz"
        asset.parent.mkdir(parents=True, exist_ok=True)
        asset.write_bytes(raw + (b"changed" if changed else b""))
        metadata = {"id": "nucleosome-2cv5", "pdb_id": "2CV5",
                    "coordinate_type": "experimental_atomic_template",
                    "asset": "assets/nucleosome-2cv5.atoms.json.gz"}
        if checksum:
            metadata["asset_sha256"] = hashlib.sha256(raw).hexdigest()
        return self.write_json("macromolecules/builtins/nucleosome-2cv5.json", metadata)

    def test_verified_builtin_nucleosome_is_reused_as_different_representation(self):
        metadata_path = self.builtin()
        plan = plan_payload(pdb_asset("2CV5"), self.target)
        self.assertEqual(plan["action"], "reuse_builtin_representation")
        self.assertFalse(plan["exact_representation_match"])
        self.assertEqual(install_payload(pdb_asset("2CV5"), self.target), metadata_path)
        self.assert_no_new_nuclear_json()

    def test_builtin_missing_or_changed_hash_is_not_claimed_verified(self):
        for checksum, changed in ((False, False), (True, True)):
            with self.subTest(checksum=checksum, changed=changed):
                self.builtin(checksum=checksum, changed=changed)
                self.assertEqual(plan_payload(pdb_asset("2CV5"), self.target)["action"], "needs_review")
                with self.assertRaises(ValueError):
                    install_payload(pdb_asset("2CV5"), self.target)
        self.assert_no_new_nuclear_json()

    def test_repeat_install_finds_identical_asset_after_its_file_is_renamed(self):
        payload = chemical()
        original_path = install_payload(payload, self.target)
        renamed = self.target / "molecules/user/renamed-source.json"
        renamed.parent.mkdir(parents=True)
        original_path.rename(renamed)
        snapshot = renamed.read_bytes()
        index = TargetIndex(self.target)
        plan = plan_payload(payload, self.target, index=index)
        self.assertEqual(plan["action"], "reuse_exact_asset")
        self.assertTrue(plan["exact_representation_match"])
        self.assertEqual(install_payload(payload, self.target, index=index), renamed)
        self.assertEqual(renamed.read_bytes(), snapshot)
        self.assertFalse(original_path.exists())
        self.assertEqual(len(list((self.target / "molecules").rglob("*.json"))), 1)

    def test_same_asset_id_with_changed_content_never_overwrites(self):
        payload = chemical()
        path = install_payload(payload, self.target)
        original = path.read_bytes()
        changed = deepcopy(payload)
        changed["coords"][1][0] = 2.0
        self.assertEqual(plan_payload(changed, self.target)["action"], "needs_review")
        with self.assertRaises(ValueError):
            install_payload(changed, self.target)
        self.assertEqual(path.read_bytes(), original)

    def test_stale_plan_rechecks_a_file_added_by_another_import(self):
        payload = chemical()
        index = TargetIndex(self.target)
        added = self.write_json("molecules/user/added-after-plan.json", payload)
        self.assertEqual(install_payload(payload, self.target, index=index), added)
        self.assert_no_new_nuclear_json()

    def test_existing_import_lock_prevents_concurrent_creation(self):
        lock = self.target / "build/locks/nuclear-library.lock"
        lock.parent.mkdir(parents=True)
        lock.write_text("another importer\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "Another nuclear import"):
            install_payload(chemical(), self.target)
        self.assertTrue(lock.exists())
        self.assert_no_new_nuclear_json()

    def test_index_from_another_checkout_cannot_redirect_installation(self):
        index = TargetIndex(self.target)
        with self.assertRaisesRegex(ValueError, "different CellSimulator"):
            install_payload(chemical(), self.target / "different", index=index)
        self.assert_no_new_nuclear_json()

    def test_partial_checkout_does_not_establish_missing_molecules(self):
        (self.target / ".git").mkdir()
        listing = subprocess.CompletedProcess([], 0, stdout=b"molecules/missing.json\0", stderr=b"")
        with mock.patch("structure_library.installation.subprocess.run", return_value=listing):
            index = TargetIndex(self.target)
        self.assertTrue(any("Incomplete checkout" in message for message in index.errors))
        self.assertEqual(plan_payload(chemical(), self.target, index=index)["action"], "needs_review")
        with self.assertRaises(ValueError):
            install_payload(chemical(), self.target, index=index)
        self.assert_no_new_nuclear_json()

    def test_reference_without_atoms_cannot_be_installed_as_empty_molecule(self):
        reference = {"asset_id": "pdb-9zzz", "source_pdb_id": "9ZZZ", "name": "Reference only"}
        with self.assertRaises(ValueError):
            install_payload(reference, self.target)
        self.assert_no_new_nuclear_json()

    def test_mismatched_coordinate_count_is_not_installed(self):
        payload = chemical()
        payload["coords"] = payload["coords"][:1]
        with self.assertRaises(ValueError):
            install_payload(payload, self.target)
        self.assert_no_new_nuclear_json()

    def test_charge_is_part_of_exact_representation_digest(self):
        original = chemical()
        changed = deepcopy(original)
        changed["total_charge"] = 1.0
        self.assertNotEqual(representation_digest(original), representation_digest(changed))
        changed = deepcopy(original)
        changed["formal_charge"] = 1
        self.assertNotEqual(representation_digest(original), representation_digest(changed))

    def write_pack(self, index, *, member_refs=()):
        path = self.target / "molecules/chebi_conformers.zip"
        path.parent.mkdir(parents=True, exist_ok=True)
        with ZipFile(path, "w") as archive:
            archive.writestr("index.json", json.dumps(index))
            for ref in member_refs:
                archive.writestr(ref + ".json", json.dumps(chemical()))
        return path

    def test_unsupported_or_malformed_packed_metadata_requires_review(self):
        valid = {"ref": "synthetic/chebi-1", "file": "synthetic/chebi-1.json",
                 "chebi_id": "CHEBI:1", "display_name": "Synthetic packed component"}
        malformed = [
            {"schema_version": 2, "entries": []},
            [],
            {"schema_version": 1, "entries": [None]},
            {"schema_version": 1, "entries": [valid, dict(valid)]},
            {"schema_version": 1, "entries": [{**valid, "ref": "../outside"}]},
        ]
        for number, index_data in enumerate(malformed):
            with self.subTest(case=number):
                self.write_pack(index_data, member_refs=("synthetic/chebi-1",))
                index = TargetIndex(self.target)
                self.assertTrue(index.errors)
                self.assertEqual(plan_payload(chemical(), self.target, index=index)["action"], "needs_review")
                with self.assertRaises(ValueError):
                    install_payload(chemical(), self.target, index=index)
        self.assert_no_new_nuclear_json()

    def test_packed_index_cannot_claim_missing_geometry_member_exists(self):
        self.write_pack({"schema_version": 1, "entries": [
            {"ref": "synthetic/chebi-1", "file": "synthetic/chebi-1.json",
             "chebi_id": "CHEBI:1", "display_name": "Synthetic packed component"},
        ]})
        index = TargetIndex(self.target)
        self.assertTrue(index.errors)
        self.assertEqual(plan_payload(chemical(), self.target, index=index)["action"], "needs_review")

    def test_corrupt_pack_is_reported_as_inspection_error(self):
        path = self.target / "molecules/chebi_conformers.zip"
        path.parent.mkdir(parents=True)
        path.write_bytes(b"invalid ZIP fixture")
        index = TargetIndex(self.target)
        self.assertTrue(index.errors)
        self.assertEqual(plan_payload(chemical(), self.target, index=index)["action"], "needs_review")
        self.assertEqual(path.read_bytes(), b"invalid ZIP fixture")


if __name__ == "__main__":
    unittest.main()
