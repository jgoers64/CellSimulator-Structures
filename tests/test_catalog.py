"""Offline catalog and static protein loader tests: no RCSB network required."""
import hashlib
import importlib.util
import json
import tempfile
import unittest
from unittest import mock
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("cellcatalog", ROOT / "catalog.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def atom_line(serial, atom, res, chain, idx, x, y, z, el, rec="ATOM"):
    return f"{rec:<6}{serial:5d} {atom:4s} {res:>3s} {chain}{idx:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{20.0:6.2f}          {el:>2s}  "


class CatalogTests(unittest.TestCase):
    def setUp(self):
        self.meta = mod.metadata(mod.catalog()[0])
        self.pdbid = self.meta["source"]["pdb_id"]
        header = "HEADER    TEST PROTEIN".ljust(62) + self.pdbid
        self.lines = [
            header,
            "MODEL        1",
            atom_line(1,"N","CYS","A",1,0.0,0.0,0.0,"N"),
            atom_line(2,"CA","CYS","A",1,1.45,0.0,0.0,"C"),
            atom_line(3,"C","CYS","A",1,2.75,0.0,0.0,"C"),
            atom_line(4,"O","CYS","A",1,3.5,0.6,0.0,"O"),
            atom_line(5,"SG","CYS","A",1,1.45,2.0,0.0,"S"),
            atom_line(6,"N","GLY","A",2,3.9,0.0,0.0,"N"),
            atom_line(7,"CA","GLY","A",2,5.2,0.0,0.0,"C"),
            atom_line(8,"C","GLY","A",2,6.5,0.0,0.0,"C"),
            atom_line(9,"O","GLY","A",2,7.1,0.6,0.0,"O"),
            atom_line(10,"SG","CYS","B",1,1.45,4.05,0.0,"S"),
            *[atom_line(i,f"C{i-10}","ALA","B",2, i*1.4,7.0,0.0,"C") for i in range(11,24)],
            atom_line(24,"FE","HEM","A",10,0.0,5.0,0.0,"FE","HETATM"),
            atom_line(25,"O","HOH","A",20,0.0,4.5,0.0,"O","HETATM"),
            "ENDMDL",
            "MODEL        2",
            atom_line(26,"N","GLY","A",1,900,900,900,"N"),
            "ENDMDL",
        ]
        self.raw = ("\n".join(self.lines) + "\n").encode()

    def test_catalog(self):
        self.assertEqual(mod.validate(), 5)
        self.assertEqual(len(mod.select_entries("all")), 5)
        with self.assertRaises(ValueError):
            mod.select_entries("../../foo")

    def test_parser_first_model_and_disulfide(self):
        parsed = mod.parse_pdb(self.raw, self.meta)
        self.assertEqual(len(parsed["atoms"]), 23)
        self.assertNotIn([900.0,900.0,900.0],parsed["coords"])
        self.assertNotIn("HOH", [a["residue_name"] for a in parsed["atoms"]])
        self.assertNotIn("HEM", [a["residue_name"] for a in parsed["atoms"]])
        sg = [i for i,a in enumerate(parsed["atoms"]) if a["atom_name"] == "SG"]
        self.assertEqual(len(sg),2)
        self.assertIn(sg, parsed["bonds"])
        self.assertIn([2,5],parsed["bonds"])

    def test_opt_in_heme(self):
        meta = json.loads(json.dumps(self.meta))
        meta["structure"]["include_hetero_residues"] = ["HEM"]
        parsed = mod.parse_pdb(self.raw, meta)
        self.assertIn("HEM",[a["residue_name"] for a in parsed["atoms"]])
        self.assertEqual(parsed["atoms"][-1]["element"],"FE")

    def test_build_structure_view_only_and_source_hash(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "test.pdb"
            path.write_bytes(self.raw)
            data = mod.build_structure(self.meta, path)
            self.assertTrue(data["view_only"])
            self.assertFalse(data["simulation_ready"])
            self.assertEqual(data["source_sha256"], hashlib.sha256(self.raw).hexdigest())
            self.assertEqual(len(data["atoms"]),len(data["coords"]))
            self.assertEqual(data["coordinate_type"],"experimental_pdb")

    def test_end_to_end_install_without_network(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            checkout = root / "CellSimulator"
            (checkout / "sim").mkdir(parents=True)
            (checkout / "main.py").write_text("# main",encoding="utf-8")
            (checkout / "sim" / "loader.py").write_text("# loader",encoding="utf-8")
            fake_pdb = root / "2hiu.pdb"
            fake_pdb.write_bytes(self.raw)
            with mock.patch.object(mod, "fetch_pdb", return_value=fake_pdb):
                output = mod.install(mod.catalog()[0], checkout)
            self.assertEqual(output.name, "insulin-2hiu.json")
            self.assertTrue(output.is_file())
            payload = json.loads(output.read_text(encoding="utf-8"))
            self.assertTrue(payload["view_only"])
            self.assertEqual(len(payload["atoms"]), 23)
            self.assertEqual(payload["source"], "RCSB PDB 2HIU; DOI 10.2210/pdb2HIU/pdb")


if __name__ == "__main__":
    unittest.main()
