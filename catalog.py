"""Curated experimental protein catalog and CellSimulator view-only installer.

Standard library only. No runtime simulation parameters are created here.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parent
MAX_PDB_BYTES = 16_000_000
RADII = {"C": 0.76, "N": 0.71, "O": 0.66, "S": 1.05, "P": 1.07, "FE": 1.26}


def catalog():
    data = json.loads((ROOT / "catalog.json").read_text(encoding="utf-8"))
    return data["entries"]


def metadata(entry):
    rel = Path(entry["metadata"])
    if rel.is_absolute() or ".." in rel.parts:
        raise ValueError(f"Unsafe metadata path: {rel}")
    return json.loads((ROOT / rel).read_text(encoding="utf-8"))


def select_entries(which):
    entries = catalog()
    if which == "all":
        return entries
    hits = [e for e in entries if e["id"] == which]
    if not hits:
        raise ValueError(f"Unknown entry {which!r}. Try 'python catalog.py list'.")
    return hits


def validate():
    entries = catalog()
    if len({e["id"] for e in entries}) != len(entries):
        raise ValueError("Duplicate catalog entry ID")
    if len({e["pdb_id"] for e in entries}) != len(entries):
        raise ValueError("Duplicate PDB entry ID")
    if not entries:
        raise ValueError("Empty catalog")
    for entry in entries:
        m = metadata(entry)
        if m["id"] != entry["id"] or m["source"]["pdb_id"] != entry["pdb_id"]:
            raise ValueError(f"Inconsistent entry: {entry['id']}")
        if not m["rendering"]["view_only"] or m["rendering"]["simulation_ready"]:
            raise ValueError(f"Unvalidated simulation model in {entry['id']}")
        if m["source"]["coordinate_url"] != f"https://files.rcsb.org/download/{entry['pdb_id']}.pdb":
            raise ValueError(f"Unexpected external coordinate URL in {entry['id']}")
    return len(entries)


def _valid_raw_pdb(raw: bytes, pdb_id: str):
    head = raw[:8000].decode("ascii", errors="replace")
    if not raw.startswith(b"HEADER") or not any(line[:6].strip() == "ATOM" for line in raw.decode("ascii", errors="replace").splitlines()):
        raise ValueError(f"PDB {pdb_id} did not contain a valid HEADER and ATOM records")
    header = head.splitlines()[0]
    deposited_id = header[62:66].strip().upper()
    if deposited_id and deposited_id != pdb_id.upper():
        raise ValueError(f"Wrong PDB downloaded: expected {pdb_id}, got {deposited_id}")


def fetch_pdb(m, cache_dir=None):
    pdb_id = m["source"]["pdb_id"].upper()
    cache_root = Path(cache_dir) if cache_dir else ROOT / "cache" / "pdb"
    cache_root.mkdir(parents=True, exist_ok=True)
    path = cache_root / f"{pdb_id.lower()}.pdb"
    if path.exists():
        raw = path.read_bytes()
        _valid_raw_pdb(raw, pdb_id)
        return path
    url = m["source"]["coordinate_url"]
    if url != f"https://files.rcsb.org/download/{pdb_id}.pdb":
        raise ValueError("Coordinate source URL is not the pinned RCSB download")
    request = Request(url, headers={"User-Agent": "CellSimulator-Structures/1.0 (research educational catalog)"})
    try:
        with urlopen(request, timeout=30) as response:
            raw = response.read(MAX_PDB_BYTES + 1)
    except OSError as exc:
        raise RuntimeError(f"Could not download {pdb_id} from RCSB. Check internet connectivity.") from exc
    if len(raw) > MAX_PDB_BYTES:
        raise ValueError(f"PDB {pdb_id} larger than configured size limit")
    _valid_raw_pdb(raw, pdb_id)
    tmp = path.with_suffix(".tmp")
    tmp.write_bytes(raw)
    tmp.replace(path)
    return path


def _secondary_ranges(lines):
    spans = []
    for line in lines:
        rec = line[:6].strip()
        try:
            if rec == "HELIX":
                start_chain, end_chain = line[19:20].strip(), line[31:32].strip()
                start, end = int(line[21:25]), int(line[33:37])
                label = "helix"
            elif rec == "SHEET":
                start_chain, end_chain = line[21:22].strip(), line[32:33].strip()
                start, end = int(line[22:26]), int(line[33:37])
                label = "sheet"
            else:
                continue
        except (ValueError, IndexError):
            continue
        if start_chain == end_chain and start_chain and start <= end:
            spans.append((start_chain, start, end, label))
    return spans


def _dist2(a, b):
    return (a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2


def _visual_bonds(atoms, coords):
    residues = {}
    for i, atom in enumerate(atoms):
        key = (atom["chain_id"], atom["residue_index"], atom["insertion_code"])
        residues.setdefault(key, []).append(i)
    bonds = set()
    for indices in residues.values():
        for offset, i in enumerate(indices):
            for j in indices[offset + 1:]:
                limit = RADII.get(atoms[i]["element"].upper(), 0.75) + RADII.get(atoms[j]["element"].upper(), 0.75) + 0.32
                d2 = _dist2(coords[i], coords[j])
                if 0.64 <= d2 <= limit*limit:
                    bonds.add((i, j))
    for (chain, res, ins), indices in residues.items():
        if ins:
            continue
        neighbors = residues.get((chain, res + 1, ""))
        if not neighbors:
            continue
        c = next((i for i in indices if atoms[i]["atom_name"] == "C"), None)
        n = next((i for i in neighbors if atoms[i]["atom_name"] == "N"), None)
        if c is not None and n is not None and 1.0 <= _dist2(coords[c], coords[n]) <= 1.8**2:
            bonds.add((min(c, n), max(c, n)))
    cysteine_sg = [i for i,a in enumerate(atoms) if a["residue_name"]=="CYS" and a["atom_name"]=="SG"]
    for offset,i in enumerate(cysteine_sg):
        for j in cysteine_sg[offset+1:]:
            if 1.8**2 <= _dist2(coords[i], coords[j]) <= 2.3**2:
                bonds.add((min(i,j),max(i,j)))
    return [list(pair) for pair in sorted(bonds)]


def parse_pdb(raw: bytes, m):
    text = raw.decode("ascii", errors="replace")
    lines = text.splitlines()
    spans = _secondary_ranges(lines)
    allowed_hetero = set(m["structure"]["include_hetero_residues"])
    atoms, coords, seen = [], [], set()
    for line in lines:
        rec = line[:6].strip()
        if rec == "ENDMDL":
            break
        if rec != "ATOM" and not (rec == "HETATM" and line[17:20].strip() in allowed_hetero):
            continue
        if len(line) < 54:
            raise ValueError("Truncated PDB coordinate line")
        if line[16:17] not in (" ", "A"):
            continue
        name, resname = line[12:16].strip(), line[17:20].strip()
        chain, insertion = line[21:22].strip(), line[26:27].strip()
        try:
            resnum = int(line[22:26])
            xyz = [float(line[i:i+8]) for i in (30,38,46)]
        except ValueError as exc:
            raise ValueError("Invalid PDB residue number or atom position") from exc
        if not all(math.isfinite(v) for v in xyz):
            raise ValueError("Nonfinite PDB atomic position")
        element = line[76:78].strip().upper() if len(line) >= 78 else ""
        if not element:
            element = name.lstrip("0123456789")[0:1].upper()
        if element in ("H", "D"):
            continue
        if not element or not name or not resname:
            raise ValueError("Invalid PDB atom metadata")
        key = (rec, chain, resnum, insertion, name)
        if key in seen:
            continue
        seen.add(key)
        sec = next((kind for c,start,end,kind in spans if c==chain and start<=resnum<=end),"loop")
        atoms.append({"element":element,"atom_name":name,"residue_name":resname,
            "residue_index":resnum-1,"chain_id":chain,"insertion_code":insertion,
            "secondary_structure":sec})
        coords.append(xyz)
    if len(atoms) < 20:
        raise ValueError(f"PDB contained only {len(atoms)} usable atoms; expected protein coordinates")
    return {"atoms":atoms,"coords":coords,"bonds":_visual_bonds(atoms,coords)}


def build_structure(m, pdb_path):
    raw = Path(pdb_path).read_bytes()
    _valid_raw_pdb(raw, m["source"]["pdb_id"])
    model = parse_pdb(raw, m)
    return {
        "name": m["name"],
        "aliases": m["aliases"],
        "coordinate_type": "experimental_pdb",
        "source": f"RCSB PDB {m['source']['pdb_id']}; DOI {m['source']['doi']}",
        "source_url": m["source"]["entry_url"],
        "source_sha256": hashlib.sha256(raw).hexdigest(),
        "view_only": True,
        "simulation_ready": False,
        "total_charge": 0.0,
        "charge_note": "0.0 is a loader placeholder, NOT a calculated molecular charge; no physics parameterization has been performed.",
        "representation": "First deposited coordinate model; protein heavy atoms and selected cofactors; source geometry unrelaxed.",
        "notes": m["notes"] + " Bonds are for visualization only, not molecular dynamics.",
        **model,
    }


def install(entry, target):
    m = metadata(entry)
    target = Path(target).expanduser().resolve()
    if not (target / "main.py").is_file() or not (target / "sim" / "loader.py").is_file():
        raise ValueError(f"Target is not an existing CellSimulator checkout: {target}")
    pdb_path = fetch_pdb(m)
    payload = build_structure(m, pdb_path)
    folder = target / "molecules" / "proteins"
    folder.mkdir(parents=True, exist_ok=True)
    filename = f"{entry['id']}-{entry['pdb_id'].lower()}.json"
    dest = folder / filename
    contents = json.dumps(payload,ensure_ascii=False,separators=(",", ":"))
    temp = dest.with_suffix(".json.tmp")
    temp.write_text(contents, encoding="utf-8")
    temp.replace(dest)
    print("Installed {} ({}) - {} heavy atoms, {} display bonds: {}".format(
        entry["name"], entry["pdb_id"], len(payload["atoms"]), len(payload["bonds"]), dest))
    return dest


def main(argv=None):
    parser = argparse.ArgumentParser(description="Curated RCSB molecular structure catalog for CellSimulator")
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("list", help="Show curated structures")
    sub.add_parser("validate", help="Validate catalog and metadata without network access")
    a = sub.add_parser("fetch", help="Download original RCSB coordinates to ignored local cache")
    a.add_argument("which", help="Entry ID or 'all'")
    b = sub.add_parser("install", help="Download and install view-only JSON into an existing CellSimulator checkout")
    b.add_argument("which", help="Entry ID or 'all'")
    b.add_argument("--target", required=True, help="Path to CellSimulator source directory")
    args = parser.parse_args(argv)
    try:
        if args.command == "list":
            for e in catalog():
                print(f"{e['priority']:2}. {e['name']}  [PDB {e['pdb_id']}]  ({e['id']})")
        elif args.command == "validate":
            print(f"Validated {validate()} curated entries (no network request).")
        elif args.command == "fetch":
            validate()
            for e in select_entries(args.which):
                print(f"Cached PDB {e['pdb_id']} at {fetch_pdb(metadata(e))}")
        elif args.command == "install":
            validate()
            for e in select_entries(args.which):
                install(e, args.target)
        return 0
    except (OSError, ValueError, RuntimeError, KeyError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
