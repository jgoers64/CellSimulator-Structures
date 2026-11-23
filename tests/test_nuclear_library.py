"""Focused offline source-integrity and coordinate-conversion regressions.

Prepared for explicit test execution only. No download or application startup is
needed. Tests requiring the optional Gemmi importer skip when it is unavailable.
The synthetic mmCIF describes a short chain with an unresolved middle residue,
one zinc ion, one water, and a two-copy assembly; it is not biological evidence.
"""
from __future__ import annotations

import gzip
import importlib.util
import io
from pathlib import Path
import tempfile
import unittest
from unittest import mock

from structure_library import coordinates


HAS_GEMMI = importlib.util.find_spec("gemmi") is not None

MMCIF = b"""data_9ZZZ
_entry.id 9ZZZ
loop_
_entity.id
_entity.type
1 polymer
2 non-polymer
3 water
_entity_poly.entity_id 1
_entity_poly.type 'polypeptide(L)'
_entity_poly.pdbx_strand_id A
loop_
_entity_poly_seq.entity_id
_entity_poly_seq.num
_entity_poly_seq.mon_id
_entity_poly_seq.hetero
1 1 ALA n
1 2 GLY n
1 3 GLY n
loop_
_struct_asym.id
_struct_asym.entity_id
A 1
B 2
C 3
loop_
_chem_comp.id
_chem_comp.type
ALA 'L-peptide linking'
GLY 'peptide linking'
ZN non-polymer
HOH non-polymer
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
ALA N CA SING
ALA CA C SING
ALA C O DOUB
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1 N N . ALA A 1 1 ? 0.125 -0.250 0.500 1.0 20.0 1 ALA A N 1
ATOM 2 C CA . ALA A 1 1 ? 1.500 -0.250 0.500 1.0 20.0 1 ALA A CA 1
ATOM 3 C C . ALA A 1 1 ? 2.800 -0.250 0.500 1.0 20.0 1 ALA A C 1
ATOM 4 O O . ALA A 1 1 ? 3.400 0.750 0.500 1.0 20.0 1 ALA A O 1
ATOM 5 N N . GLY A 1 3 ? 4.100 -0.250 0.500 1.0 20.0 3 GLY A N 1
HETATM 6 ZN ZN . ZN B 2 . ? 30.000 1.000 2.000 1.0 20.0 10 ZN B ZN 1
HETATM 7 O O . HOH C 3 . ? 40.000 2.000 3.000 1.0 20.0 1 HOH C O 1
_pdbx_struct_assembly.id 1
_pdbx_struct_assembly.details 'synthetic test assembly'
_pdbx_struct_assembly.oligomeric_details dimeric
_pdbx_struct_assembly.oligomeric_count 2
_pdbx_struct_assembly_gen.assembly_id 1
_pdbx_struct_assembly_gen.oper_expression '(1,2)'
_pdbx_struct_assembly_gen.asym_id_list A
loop_
_pdbx_struct_oper_list.id
_pdbx_struct_oper_list.type
_pdbx_struct_oper_list.matrix[1][1]
_pdbx_struct_oper_list.matrix[1][2]
_pdbx_struct_oper_list.matrix[1][3]
_pdbx_struct_oper_list.matrix[2][1]
_pdbx_struct_oper_list.matrix[2][2]
_pdbx_struct_oper_list.matrix[2][3]
_pdbx_struct_oper_list.matrix[3][1]
_pdbx_struct_oper_list.matrix[3][2]
_pdbx_struct_oper_list.matrix[3][3]
_pdbx_struct_oper_list.vector[1]
_pdbx_struct_oper_list.vector[2]
_pdbx_struct_oper_list.vector[3]
1 'identity operation' 1 0 0 0 1 0 0 0 1 0 0 0
2 'test translation' 1 0 0 0 1 0 0 0 1 10 0 0
#
"""


def source_metadata(raw=MMCIF, *, pdb_id="9ZZZ", assembly_id=None):
    return {
        "id": "pdb-9zzz",
        "name": "Synthetic incomplete-chain fixture",
        "provenance_class": "experimental_coordinate_model",
        "assembly": {"id": assembly_id},
        "source": {
            "database": "test fixture",
            "pdb_id": pdb_id,
            "coordinate_url": "https://files.rcsb.org/download/9ZZZ.cif",
            "entry_url": "https://www.rcsb.org/structure/9ZZZ",
            "coordinate_sha256": coordinates.digest(raw),
        },
    }


class SourceIntegrityTests(unittest.TestCase):
    def test_repository_paths_reject_escape_and_absolute_paths(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder) / "repo"
            root.mkdir()
            self.assertEqual(coordinates.safe_path(root, "structures/a.json"),
                             root / "structures" / "a.json")
            for relative in ("../outside.json", str(Path(folder) / "outside.json")):
                with self.subTest(path=relative), self.assertRaises(ValueError):
                    coordinates.safe_path(root, relative)

    def test_repository_path_rejects_symlink_escape(self):
        with tempfile.TemporaryDirectory() as folder:
            root, outside = Path(folder) / "repo", Path(folder) / "outside"
            root.mkdir()
            outside.mkdir()
            try:
                (root / "linked").symlink_to(outside, target_is_directory=True)
            except (OSError, NotImplementedError):
                self.skipTest("Creating symlinks is unavailable on this platform")
            with self.assertRaises(ValueError):
                coordinates.safe_path(root, "linked/source.cif")

    def test_decompression_bounds_cover_concatenated_members(self):
        self.assertEqual(coordinates.bounded_decompress(gzip.compress(b"a" * 32), 32),
                         b"a" * 32)
        too_large = (gzip.compress(b"a" * 33),
                     gzip.compress(b"a" * 16) + gzip.compress(b"b" * 17))
        for raw in too_large:
            with self.subTest(compressed_bytes=len(raw)), self.assertRaises(ValueError):
                coordinates.bounded_decompress(raw, 32)

    def test_changed_download_is_not_cached(self):
        with tempfile.TemporaryDirectory() as folder:
            with mock.patch.object(coordinates, "download", return_value=b"wrong model"):
                with self.assertRaisesRegex(ValueError, "Archive content changed"):
                    coordinates.fetch_source(source_metadata(), folder)
            self.assertFalse((Path(folder) / "cache/nuclear/pdb-9zzz.cif").exists())

    def test_changed_cache_is_not_silently_replaced(self):
        with tempfile.TemporaryDirectory() as folder:
            cached = Path(folder) / "cache/nuclear/pdb-9zzz.cif"
            cached.parent.mkdir(parents=True)
            cached.write_bytes(b"user changed this source")
            with mock.patch.object(coordinates, "download") as fetch:
                with self.assertRaisesRegex(ValueError, "Cached source differs"):
                    coordinates.fetch_source(source_metadata(), folder)
            fetch.assert_not_called()
            self.assertEqual(cached.read_bytes(), b"user changed this source")

    def test_download_limit_applies_without_content_length(self):
        class Response(io.BytesIO):
            def geturl(self):
                return "https://files.rcsb.org/download/9ZZZ.cif"

        with mock.patch.object(coordinates, "urlopen", return_value=Response(b"x" * 33)):
            with self.assertRaisesRegex(ValueError, "byte limit"):
                coordinates.download("https://files.rcsb.org/download/9ZZZ.cif", max_bytes=32)


@unittest.skipUnless(HAS_GEMMI, "Optional Gemmi importer is not installed")
class CoordinateConversionTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.path = Path(self.directory.name) / "fixture.cif"
        self.path.write_bytes(MMCIF)

    def test_correct_hash_does_not_excuse_wrong_entry_identity(self):
        with self.assertRaisesRegex(ValueError, "entry ID"):
            coordinates.read_structure(self.path, source_metadata(pdb_id="1VYM"))

    def test_assembly_budget_rejects_before_expansion(self):
        import gemmi
        # Five chain-A atoms expanded twice exceed nine, even though the raw
        # asymmetric unit itself has only seven atoms including ion and water.
        with mock.patch.object(gemmi, "make_assembly") as expand:
            with self.assertRaisesRegex(ValueError, "Assembly requires"):
                coordinates.convert_mmcif(source_metadata(assembly_id="1"), self.path, max_atoms=9)
        expand.assert_not_called()

    def test_source_geometry_gap_and_metal_are_preserved(self):
        payload = coordinates.convert_mmcif(source_metadata(), self.path)
        self.assertEqual(payload["coordinate_units"], "angstrom")
        self.assertEqual(payload["coords"][0], [0.125, -0.25, 0.5])
        self.assertEqual(payload["coords"][1], [1.5, -0.25, 0.5])
        self.assertEqual(len(payload["atoms"]), 6)
        self.assertEqual(payload["atoms"][-1]["element"], "Zn")
        self.assertEqual(payload["coords"][-1], [30.0, 1.0, 2.0])
        observed = {atom["label_seq_id"] for atom in payload["atoms"]
                    if atom["label_asym_id"] == "A"}
        self.assertEqual(observed, {1, 3})
        self.assertNotIn([2, 4], payload["bonds"])
        self.assertEqual(payload["processing"]["missing_atoms_added"], 0)
        self.assertEqual(payload["processing"]["hydrogens_added"], 0)
        self.assertEqual(payload["processing"]["omitted_water_atoms"], 1)
        self.assertTrue(payload["view_only"])
        self.assertFalse(payload["simulation_ready"])


if __name__ == "__main__":
    unittest.main()
