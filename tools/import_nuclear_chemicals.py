"""Import explicitly identified CCD chemical states, without inventing nuclear poses.

The wwPDB Chemical Component Dictionary supplies chemical graphs and computed
ideal coordinates. Those coordinates are NOT experimentally measured positions
inside a nucleus. No molecular dynamics parameters or partial charges are made.

Rebuild from the previously reviewed raw files (no network):
    python tools/import_nuclear_chemicals.py --offline

Initial import can use --source-cache DIR and --inventory FILE. Downloads are
limited to fixed files.rcsb.org CCD endpoints, at most 2 MB per component.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]
MAX_SOURCE_BYTES = 2_000_000
MAX_ATOMS = 500
MANIFEST = "collections/human-nucleus/chemical-assets.json"
CCD_DOCUMENTATION = "https://www.wwpdb.org/data/ccd"

# These IDs were checked against CCD names, formulae and stereochemical
# descriptors. Familiar biochemical abbreviations are NOT necessarily CCD IDs:
# GMP is guanosine, CMP is cyclic AMP, MAL is maltose, SUC is sucrose, etc.
# Formula and formal charge specify the representative state actually imported.
SPECS = [
    ("Water", "HOH", "H2 O", 0, "Water"),
    ("Potassium (K+)", "K", "K", 1, "Potassium ion"),
    ("Chloride (Cl-)", "CL", "Cl", -1, "Chloride ion"),
    ("Sodium (Na+)", "NA", "Na", 1, "Sodium ion"),
    ("Magnesium (Mg2+)", "MG", "Mg", 2, "Magnesium ion"),
    ("Calcium (Ca2+)", "CA", "Ca", 2, "Calcium ion"),
    ("ATP", "ATP", "C10 H16 N5 O13 P3", 0, "ATP, CCD neutral acid state"),
    ("ADP", "ADP", "C10 H15 N5 O10 P2", 0, "ADP, CCD neutral acid state"),
    ("AMP", "AMP", "C10 H14 N5 O7 P", 0, "AMP, CCD neutral acid state"),
    ("GTP", "GTP", "C10 H16 N5 O14 P3", 0, "GTP, CCD neutral acid state"),
    ("GDP", "GDP", "C10 H15 N5 O11 P2", 0, "GDP, CCD neutral acid state"),
    ("GMP", "5GP", "C10 H14 N5 O8 P", 0, "GMP, CCD neutral acid state"),
    ("CTP", "CTP", "C9 H16 N3 O14 P3", 0, "CTP, CCD neutral acid state"),
    ("CDP", "CDP", "C9 H15 N3 O11 P2", 0, "CDP, CCD neutral acid state"),
    ("CMP", "C5P", "C9 H14 N3 O8 P", 0, "CMP, CCD neutral acid state"),
    ("UTP", "UTP", "C9 H15 N2 O15 P3", 0, "UTP, CCD neutral acid state"),
    ("UDP", "UDP", "C9 H14 N2 O12 P2", 0, "UDP, CCD neutral acid state"),
    ("UMP", "U5P", "C9 H13 N2 O9 P", 0, "UMP, CCD neutral acid state"),
    ("dATP", "DTP", "C10 H16 N5 O12 P3", 0, "dATP, CCD neutral acid state"),
    ("dGTP", "DGT", "C10 H16 N5 O13 P3", 0, "dGTP, CCD neutral acid state"),
    ("dCTP", "DCP", "C9 H16 N3 O13 P3", 0, "dCTP, CCD neutral acid state"),
    ("dTTP", "TTP", "C10 H17 N2 O14 P3", 0, "dTTP, CCD neutral acid state"),
    ("Free NAD+", "NAD", "C21 H27 N7 O14 P2", 0, "NAD+, CCD net-neutral protonation state"),
    ("NADH", "NAI", "C21 H29 N7 O14 P2", 0, "NADH, CCD neutral protonation state"),
    ("NADP+", "NAP", "C21 H28 N7 O17 P3", 0, "NADP+, CCD net-neutral protonation state"),
    ("NADPH", "NDP", "C21 H30 N7 O17 P3", 0, "NADPH, CCD neutral protonation state"),
    ("Nicotinamide", "NCA", "C6 H6 N2 O", 0, "Nicotinamide"),
    ("S-adenosylmethionine", "SAM", "C15 H22 N6 O5 S", 0, "S-adenosylmethionine, CCD net-neutral state"),
    ("S-adenosylhomocysteine", "SAH", "C14 H20 N6 O5 S", 0, "S-adenosyl-L-homocysteine"),
    ("Coenzyme A", "COA", "C21 H36 N7 O16 P3 S", 0, "Coenzyme A, CCD neutral acid state"),
    ("Acetyl-CoA", "ACO", "C23 H38 N7 O17 P3 S", 0, "Acetyl-CoA, CCD neutral acid state"),
    ("Propionyl-CoA", "1VU", "C24 H40 N7 O17 P3 S", 0, "Propionyl-CoA, CCD neutral acid state"),
    ("Reduced glutathione (GSH)", "GSH", "C10 H17 N3 O6 S", 0, "Reduced glutathione, CCD neutral acid state"),
    ("Oxidized glutathione (GSSG)", "GDS", "C20 H32 N6 O12 S2", 0, "Oxidized glutathione, CCD neutral acid state"),
    ("Spermidine", "SPD", "C7 H19 N3", 0, "Spermidine, CCD unprotonated state"),
    ("Spermine", "SPM", "C10 H26 N4", 0, "Spermine, CCD unprotonated state"),
    ("Putrescine", "PUT", "C4 H12 N2", 0, "Putrescine, CCD unprotonated state"),
    ("Pyruvate", "PYR", "C3 H4 O3", 0, "Pyruvic acid, neutral representative of pyruvate family"),
    ("Citrate", "CIT", "C6 H8 O7", 0, "Citric acid, neutral representative of citrate family"),
    ("Isocitrate", "ICT", "C6 H8 O7", 0, "(1R,2S)-Isocitric acid, neutral representative"),
    ("2-Oxoglutarate (alpha-ketoglutarate)", "AKG", "C5 H6 O5", 0, "2-Oxoglutaric acid, neutral representative"),
    ("Succinate", "SIN", "C4 H6 O4", 0, "Succinic acid, neutral representative of succinate family"),
    ("Fumarate", "FUM", "C4 H4 O4", 0, "Fumaric acid, neutral representative of fumarate family"),
    ("Malate", "LMR", "C4 H6 O5", 0, "L-Malic acid, neutral representative of malate family"),
    ("Inorganic phosphate", "PO4", "O4 P", -3, "Phosphate, CCD PO4(3-) state"),
    ("Pyrophosphate", "PPV", "H4 O7 P2", 0, "Pyrophosphoric acid, neutral representative"),
    ("Bicarbonate", "BCT", "C H O3", -1, "Bicarbonate ion"),
    ("Dissolved oxygen", "OXY", "O2", 0, "Dioxygen, CCD ideal bond geometry"),
    ("Dissolved carbon dioxide", "CO2", "C O2", 0, "Carbon dioxide"),
    ("Cholesterol", "CLR", "C27 H46 O", 0, "Cholesterol"),
]

AMINO_ACIDS = [
    ("ALA", "C3 H7 N O2", 0, "L-Alanine"),
    ("ARG", "C6 H15 N4 O2", 1, "L-Arginine"),
    ("ASN", "C4 H8 N2 O3", 0, "L-Asparagine"),
    ("ASP", "C4 H7 N O4", 0, "L-Aspartic acid"),
    ("CYS", "C3 H7 N O2 S", 0, "L-Cysteine"),
    ("GLN", "C5 H10 N2 O3", 0, "L-Glutamine"),
    ("GLU", "C5 H9 N O4", 0, "L-Glutamic acid"),
    ("GLY", "C2 H5 N O2", 0, "Glycine"),
    ("HIS", "C6 H10 N3 O2", 1, "L-Histidine"),
    ("ILE", "C6 H13 N O2", 0, "L-Isoleucine"),
    ("LEU", "C6 H13 N O2", 0, "L-Leucine"),
    ("LYS", "C6 H15 N2 O2", 1, "L-Lysine"),
    ("MET", "C5 H11 N O2 S", 0, "L-Methionine"),
    ("PHE", "C9 H11 N O2", 0, "L-Phenylalanine"),
    ("PRO", "C5 H9 N O2", 0, "L-Proline"),
    ("SER", "C3 H7 N O3", 0, "L-Serine"),
    ("THR", "C4 H9 N O3", 0, "L-Threonine"),
    ("TRP", "C11 H12 N2 O2", 0, "L-Tryptophan"),
    ("TYR", "C9 H11 N O3", 0, "L-Tyrosine"),
    ("VAL", "C5 H11 N O2", 0, "L-Valine"),
]


def cif_tokens(text):
    """Tokenize this bounded CCD mmCIF subset, preserving quoted atom names."""
    i, length = 0, len(text)
    while i < length:
        if text[i].isspace():
            i += 1
            continue
        if text[i] == "#":
            i = text.find("\n", i)
            if i < 0:
                return
            continue
        if text[i] == ";" and (i == 0 or text[i - 1] == "\n"):
            end = text.find("\n;", i + 1)
            if end < 0:
                raise ValueError("Unterminated mmCIF text field")
            yield text[i + 1:end].strip("\n")
            i = end + 2
        elif text[i] in "\"'":
            quote, start = text[i], i + 1
            i = start
            while i < length:
                if text[i] == quote and (i + 1 == length or text[i + 1].isspace()):
                    yield text[start:i]
                    i += 1
                    break
                i += 1
            else:
                raise ValueError("Unterminated quoted mmCIF field")
        else:
            start = i
            while i < length and not text[i].isspace():
                i += 1
            yield text[start:i]


def parse_ccd(raw):
    if len(raw) > MAX_SOURCE_BYTES:
        raise ValueError("CCD input exceeds bounded import size")
    tokens = list(cif_tokens(raw.decode("utf-8")))
    scalars, loops, i, blocks = {}, {}, 0, 0
    while i < len(tokens):
        token = tokens[i]
        if token.startswith("data_"):
            blocks += 1
            if blocks != 1:
                raise ValueError("Expected one CCD data block")
            i += 1
        elif token == "loop_":
            i += 1
            tags = []
            while i < len(tokens) and tokens[i].startswith("_"):
                tags.append(tokens[i]); i += 1
            if not tags:
                raise ValueError("Empty mmCIF loop")
            values = []
            while i < len(tokens) and not (tokens[i].startswith(("_", "data_", "save_")) or tokens[i] in ("loop_", "stop_")):
                values.append(tokens[i]); i += 1
            if len(values) % len(tags):
                raise ValueError("Incomplete mmCIF loop row")
            category = tags[0].split(".")[0]
            loops[category] = [dict(zip(tags, values[j:j + len(tags)])) for j in range(0, len(values), len(tags))]
        elif token.startswith("_"):
            if i + 1 >= len(tokens):
                raise ValueError("Missing mmCIF scalar value")
            scalars[token] = tokens[i + 1]; i += 2
        else:
            raise ValueError(f"Unsupported mmCIF control token: {token}")
    for category in ("_chem_comp_atom", "_chem_comp_bond", "_pdbx_chem_comp_descriptor"):
        one = {key: value for key, value in scalars.items() if key.startswith(category + ".")}
        if one:
            if category in loops:
                raise ValueError("Mixed scalar and loop category")
            loops[category] = [one]
    return scalars, loops


def formula_counts(formula):
    result = Counter()
    for token in formula.split():
        match = re.fullmatch(r"([A-Z][a-z]?)(\d*)", token)
        if not match:
            raise ValueError(f"Unsupported formula token: {token}")
        result[match[1]] += int(match[2] or 1)
    return result


def build_asset(raw, spec, raw_path, source_date):
    inventory_name, ccd_id, expected_formula, expected_charge, title = spec
    scalars, loops = parse_ccd(raw)
    value = lambda key: scalars.get("_chem_comp." + key)
    if value("id") != ccd_id or value("pdbx_release_status") != "REL":
        raise ValueError(f"Wrong or obsolete CCD record for {ccd_id}")
    if value("formula") != expected_formula or int(value("pdbx_formal_charge")) != expected_charge:
        raise ValueError(f"Chemical state differs from reviewed identity for {ccd_id}: {value('formula')}, {value('pdbx_formal_charge')}")
    rows = loops.get("_chem_comp_atom", [])
    if not 1 <= len(rows) <= MAX_ATOMS:
        raise ValueError(f"Invalid or excessive CCD atom count for {ccd_id}")
    atoms, coords, names = [], [], {}
    for row in rows:
        get = lambda key: row.get("_chem_comp_atom." + key)
        name, element = get("atom_id"), get("type_symbol")
        if not name or name in names or not re.fullmatch(r"[A-Z][a-z]?", element.title()):
            raise ValueError(f"Invalid or duplicate CCD atom in {ccd_id}")
        xyz = [float(get("pdbx_model_Cartn_" + axis + "_ideal")) for axis in "xyz"]
        if not all(math.isfinite(number) for number in xyz):
            raise ValueError(f"Missing or nonfinite ideal coordinates for {ccd_id}")
        names[name] = len(atoms)
        atoms.append({"element": element.title(), "atom_name": name, "formal_charge": int(get("charge")),
                      "stereochemistry": get("pdbx_stereo_config"), "aromatic": get("pdbx_aromatic_flag") == "Y",
                      "ccd_leaving_atom": get("pdbx_leaving_atom_flag") == "Y"})
        coords.append(xyz)
    if Counter(atom["element"] for atom in atoms) != formula_counts(expected_formula):
        raise ValueError(f"CCD atom inventory does not match formula for {ccd_id}")
    if sum(atom["formal_charge"] for atom in atoms) != expected_charge:
        raise ValueError(f"CCD atom charges do not sum to formal charge for {ccd_id}")
    bonds, orders, bond_details, unique = [], [], [], set()
    order_map = {"SING": 1, "DOUB": 2, "TRIP": 3, "AROM": 1.5}
    for row in loops.get("_chem_comp_bond", []):
        get = lambda key: row.get("_chem_comp_bond." + key)
        pair = [names[get("atom_id_1")], names[get("atom_id_2")]]
        key = tuple(sorted(pair))
        if key in unique or pair[0] == pair[1] or get("value_order") not in order_map:
            raise ValueError(f"Invalid or unsupported CCD bond in {ccd_id}")
        unique.add(key); bonds.append(pair); orders.append(order_map[get("value_order")])
        bond_details.append({"source_order": get("value_order"), "aromatic": get("pdbx_aromatic_flag") == "Y", "stereochemistry": get("pdbx_stereo_config")})
    if len(atoms) > 1 and not bonds:
        raise ValueError(f"Missing chemical graph in {ccd_id}")
    # Every representative in this curated set is one covalently connected
    # molecule or one monatomic ion. An isolated retained hydrogen is therefore
    # an incomplete CCD graph, not a permissible second component. A future
    # salt/coordination-complex importer must model its components explicitly.
    adjacency = [set() for _ in atoms]
    for first, second in bonds:
        adjacency[first].add(second)
        adjacency[second].add(first)
    visited, pending = set(), [0]
    while pending:
        current = pending.pop()
        if current in visited:
            continue
        visited.add(current)
        pending.extend(adjacency[current] - visited)
    if len(visited) != len(atoms):
        raise ValueError(f"Disconnected chemical graph in {ccd_id}; expected one component")
    descriptors = [{key.split(".", 1)[1]: val for key, val in row.items()}
                   for row in loops.get("_pdbx_chem_comp_descriptor", [])]
    state_note = ("This exact CCD formula, formal charge and stereochemistry define a representative chemical state. "
                  "It is not a selected physiological protonation state or a measurement of nuclear abundance. "
                  "Bound ions, Mg-nucleotide complexes, hydration shells and binding-induced conformations are not supplied.")
    if inventory_name == "Free proteinogenic amino acids":
        state_note += (" The full CCD component includes its terminal/leaving atoms; this asset is a free-component representative, "
                       "not a peptide residue. This expands the family to the 20 canonical amino acids only.")
    if ccd_id in ("NAD", "NAP"):
        state_note += " The biochemical '+' identifies the oxidized pyridinium form; the CCD total formal charge is separately recorded."
    if ccd_id == "OXY":
        state_note += " This bond graph and static geometry do not model oxygen's electronic spin state."
    source = {"database": "wwPDB Chemical Component Dictionary", "ccd_id": ccd_id,
              "url": f"https://files.rcsb.org/ligands/download/{ccd_id}.cif",
              "entry_url": f"https://www.rcsb.org/ligand/{ccd_id}", "documentation_url": CCD_DOCUMENTATION,
              "source_path": raw_path, "sha256": hashlib.sha256(raw).hexdigest(), "retrieved_at_utc": source_date,
              "ccd_modified_date": value("pdbx_modified_date"), "ideal_coordinate_generator": value("pdbx_ideal_coordinates_details")}
    asset_id = "chemical-ccd-" + ccd_id.lower()
    metadata = {"schema_version": 1, "id": asset_id, "name": title, "inventory_name": inventory_name,
                "source_kind": "wwpdb_ccd", "provenance_class": "idealized_chemical_component", "source": source,
                "chemical_identity": {"ccd_name": value("name"), "formula": expected_formula, "formal_charge": expected_charge,
                                      "component_type": value("type"), "descriptors": descriptors},
                "coordinate_selection": "all CCD ideal-coordinate atoms, including dictionary hydrogens and leaving atoms",
                "coordinate_units": "angstrom", "view_only": True, "simulation_ready": False,
                "atom_count": len(atoms), "bond_count": len(bonds), "chemical_state_note": state_note,
                "review_status": "dictionary_identity_and_import_fidelity_checked",
                "scientific_validation": "No experimental nuclear pose or solution-state conformer validation; ideal coordinates are computational.",
                "copies_per_default_g1_nucleus": None, "abundance_rank_default_g1": None}
    coordinate = {"schema_version": 1, "id": asset_id, "name": title, "aliases": [inventory_name],
                  "coordinate_type": "idealized_chemical_component", "coordinate_units": "angstrom", "units": "angstrom",
                  "view_only": True, "simulation_ready": False, "source": f"wwPDB CCD {ccd_id}, computed ideal coordinates",
                  "source_url": source["entry_url"], "source_sha256": source["sha256"], "formal_charge": expected_charge,
                  "total_charge": float(expected_charge), "charge_note": "Exact CCD formal charge; no partial charges or force-field parameters assigned.",
                  "formula": expected_formula, "notes": state_note, "atoms": atoms, "coords": coords, "bonds": bonds,
                  "bond_orders": orders, "bond_details": bond_details}
    return metadata, coordinate


def write_json(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_bytes((json.dumps(payload, ensure_ascii=False, indent=2) + "\n").encode("utf-8"))
    temporary.replace(path)


def load_source(root, ccd_id, cache, offline):
    folder = root / "structures" / "chemicals" / ("ccd-" + ccd_id.lower())
    path = folder / "source.cif"
    meta_path = folder / "metadata.json"
    expected = None
    source_date = datetime.now(timezone.utc).isoformat()
    if meta_path.exists():
        old = json.loads(meta_path.read_text(encoding="utf-8"))
        expected = old["source"]["sha256"]
        source_date = old["source"]["retrieved_at_utc"]
    if path.exists():
        raw = path.read_bytes()
    elif cache and (cache / (ccd_id + ".cif")).exists():
        raw = (cache / (ccd_id + ".cif")).read_bytes()
    elif offline:
        raise FileNotFoundError(f"Reviewed source unavailable offline: {ccd_id}")
    else:
        request = Request(f"https://files.rcsb.org/ligands/download/{ccd_id}.cif",
                          headers={"User-Agent": "CellSimulator-Structures/2.0 (educational chemical catalog)"})
        with urlopen(request, timeout=30) as response:
            raw = response.read(MAX_SOURCE_BYTES + 1)
    if len(raw) > MAX_SOURCE_BYTES or (expected and hashlib.sha256(raw).hexdigest() != expected):
        raise ValueError(f"Source size or pinned hash mismatch for {ccd_id}; review the new source before replacing it")
    return raw, path, source_date


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--source-cache", type=Path)
    parser.add_argument("--inventory", type=Path, help="Original research inventory or repository inventory JSON")
    parser.add_argument("--offline", action="store_true", help="Only use reviewed files already in the repository")
    args = parser.parse_args(argv)
    root = args.root.resolve()
    specs = SPECS + [("Free proteinogenic amino acids", ccd, formula, charge, name + ", CCD free-component state")
                     for ccd, formula, charge, name in AMINO_ACIDS]
    assets = []
    for spec in specs:
        raw, path, source_date = load_source(root, spec[1], args.source_cache, args.offline)
        metadata, coordinate = build_asset(raw, spec, path.relative_to(root).as_posix(), source_date)
        folder = path.parent
        folder.mkdir(parents=True, exist_ok=True)
        path.write_bytes(raw)
        write_json(folder / "metadata.json", metadata)
        write_json(folder / "structure.json", coordinate)
        assets.append({"id": metadata["id"], "name": metadata["name"], "inventory_name": metadata["inventory_name"],
                       "metadata_path": (folder / "metadata.json").relative_to(root).as_posix(),
                       "coordinate_path": (folder / "structure.json").relative_to(root).as_posix(),
                       "source_kind": metadata["source_kind"], "provenance_class": metadata["provenance_class"],
                       "source_sha256": metadata["source"]["sha256"], "atom_count": metadata["atom_count"],
                       "coordinate_sha256": hashlib.sha256((folder / "structure.json").read_bytes()).hexdigest(),
                       "formal_charge": metadata["chemical_identity"]["formal_charge"],
                       "availability": "bundled_coordinates", "population_assignment": "not_assigned"})
        print(f"Imported {spec[1]}: {metadata['atom_count']} atoms, {metadata['name']}", flush=True)
    manifest_path = root / MANIFEST
    inventory_path = args.inventory or root / "collections" / "human-nucleus" / "inventory.json"
    if inventory_path.exists():
        original = json.loads(inventory_path.read_text(encoding="utf-8"))
        if "other_molecules_and_families" in original:
            inventory_rows = original["other_molecules_and_families"]
        else:
            inventory_rows = [entry["source_record"] for entry in original["entries"] if entry["kind"] != "protein_gene"]
        remaining = []
        mapped = {s[0] for s in specs}
        for row in inventory_rows:
            if row["name"] in mapped:
                continue
            reason = ("Protein-bound zinc requires the source protein's coordination geometry; no isolated-ion asset substitutes for that site."
                      if row["name"] == "Protein-bound zinc" else
                      "This family, sequence/state or lipid class needs a specific identity and suitable structural source; no generic atomic geometry is assigned.")
            remaining.append({"inventory_name": row["name"], "status": "requires_specific_structure", "reason": reason})
    elif manifest_path.exists():
        remaining = json.loads(manifest_path.read_text(encoding="utf-8"))["unresolved"]
    else:
        raise ValueError("First import needs --inventory to preserve unresolved inventory entries")
    manifest = {"schema_version": 1, "title": "Human nucleus chemical-state representatives",
                "coordinate_units": "angstrom", "assets": assets, "unresolved": remaining,
                "inventory_rows_with_representatives": len({a["inventory_name"] for a in assets}),
                "asset_count": len(assets), "all_inventory_molecules_structurally_resolved": False,
                "notes": ["All coordinates are computed CCD ideal coordinates, not experimentally measured nuclear conformations.",
                          "Representative protonation states are explicitly named. Do not use their charges or populations as a nucleus recipe.",
                          "All dictionary hydrogens and terminal/leaving atoms are retained; no coordinate geometry or bonds are invented.",
                          "RNA/DNA families and unspecified lipid classes remain unresolved here; macromolecular assemblies are cataloged separately.",
                          "Successful identity/hash/graph checks establish import fidelity, not 100% scientific accuracy."]}
    write_json(manifest_path, manifest)
    print(json.dumps({"assets": len(assets), "inventory_rows": manifest["inventory_rows_with_representatives"], "unresolved": len(remaining)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
