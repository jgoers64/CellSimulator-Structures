"""Offline CCD identity and chemical-graph regressions; run only on request."""
from __future__ import annotations

from collections import Counter
import importlib.util
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "nuclear_chemical_importer", ROOT / "tools" / "import_nuclear_chemicals.py")
chemicals = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(chemicals)

# Synthetic, minimal CCD water fixture with the dictionary's ideal-coordinate
# convention. It contains a distinct model-coordinate field so the test also
# guards against accidentally substituting the bound source pose.
WATER = b"""data_HOH
_chem_comp.id HOH
_chem_comp.name WATER
_chem_comp.type NON-POLYMER
_chem_comp.formula 'H2 O'
_chem_comp.pdbx_formal_charge 0
_chem_comp.pdbx_release_status REL
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.charge
_chem_comp_atom.pdbx_stereo_config
_chem_comp_atom.pdbx_aromatic_flag
_chem_comp_atom.pdbx_leaving_atom_flag
_chem_comp_atom.model_Cartn_x
_chem_comp_atom.model_Cartn_y
_chem_comp_atom.model_Cartn_z
_chem_comp_atom.pdbx_model_Cartn_x_ideal
_chem_comp_atom.pdbx_model_Cartn_y_ideal
_chem_comp_atom.pdbx_model_Cartn_z_ideal
HOH O O 0 N N N 99 98 97 -0.064 0.000 0.000
HOH H1 H 0 N N N 99 98 97 0.512 0.000 -0.776
HOH H2 H 0 N N N 99 98 97 0.512 0.000 0.776
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
_chem_comp_bond.pdbx_stereo_config
HOH O H1 SING N N
HOH O H2 SING N N
#
"""


def water_asset(raw=WATER):
    return chemicals.build_asset(
        raw, ("Water", "HOH", "H2 O", 0, "Water"),
        "fixture/source.cif", "2000-01-01T00:00:00Z")


class NuclearChemicalTests(unittest.TestCase):
    def test_biochemical_names_do_not_become_unrelated_ccd_ids(self):
        mapping = {row[0]: row[1] for row in chemicals.SPECS}
        # CCD GMP is guanosine, CMP is cyclic AMP, MAL is maltose, and SUC
        # is sucrose. These explicit mappings protect chemically distinct names.
        for name, ccd_id in {"GMP": "5GP", "CMP": "C5P", "Malate": "LMR",
                             "Succinate": "SIN", "Propionyl-CoA": "1VU"}.items():
            with self.subTest(name=name):
                self.assertEqual(mapping[name], ccd_id)

    def test_water_uses_ideal_pose_and_preserves_authoritative_graph(self):
        metadata, payload = water_asset()
        self.assertEqual(Counter(atom["element"] for atom in payload["atoms"]), {"H": 2, "O": 1})
        self.assertEqual(payload["coords"][0], [-0.064, 0.0, 0.0])
        self.assertNotIn([99.0, 98.0, 97.0], payload["coords"])
        names = [atom["atom_name"] for atom in payload["atoms"]]
        graph = {frozenset((names[a], names[b])) for a, b in payload["bonds"]}
        self.assertEqual(graph, {frozenset(("O", "H1")), frozenset(("O", "H2"))})
        self.assertEqual(payload["bond_orders"], [1, 1])
        self.assertEqual(payload["formal_charge"], 0)
        self.assertEqual(payload["coordinate_type"], "idealized_chemical_component")
        self.assertIsNone(metadata["copies_per_default_g1_nucleus"])
        self.assertFalse(payload["simulation_ready"])

    def test_source_formula_and_formal_charge_must_match_reviewed_state(self):
        variants = [
            WATER.replace(b"_chem_comp.formula 'H2 O'", b"_chem_comp.formula 'H3 O'"),
            WATER.replace(b"_chem_comp.pdbx_formal_charge 0", b"_chem_comp.pdbx_formal_charge 1"),
        ]
        for raw in variants:
            with self.subTest(raw=raw[:120]), self.assertRaisesRegex(ValueError, "Chemical state differs"):
                water_asset(raw)

    def test_atom_inventory_cannot_disagree_with_source_formula(self):
        changed = WATER.replace(b"HOH H2 H 0 N N N", b"HOH H2 O 0 N N N")
        with self.assertRaisesRegex(ValueError, "atom inventory"):
            water_asset(changed)

    def test_atom_charges_cannot_disagree_with_component_formal_charge(self):
        changed = WATER.replace(b"HOH O O 0 N N N", b"HOH O O 1 N N N")
        with self.assertRaisesRegex(ValueError, "charges do not sum"):
            water_asset(changed)

    def test_disconnected_hydrogen_rejects_incomplete_chemical_graph(self):
        changed = WATER.replace(b"HOH O H2 SING N N\n", b"")
        with self.assertRaisesRegex(ValueError, "Disconnected chemical graph"):
            water_asset(changed)

    def test_monatomic_magnesium_is_valid_without_covalent_bonds(self):
        raw = (ROOT / "structures/chemicals/ccd-mg/source.cif").read_bytes()
        spec = next(row for row in chemicals.SPECS if row[1] == "MG")
        metadata, payload = chemicals.build_asset(
            raw, spec, "structures/chemicals/ccd-mg/source.cif", "2000-01-01T00:00:00Z")
        self.assertEqual(len(payload["atoms"]), 1)
        self.assertEqual(payload["atoms"][0]["element"], "Mg")
        self.assertEqual(payload["formal_charge"], 2)
        self.assertEqual(payload["bonds"], [])
        self.assertEqual(metadata["bond_count"], 0)


if __name__ == "__main__":
    unittest.main()
