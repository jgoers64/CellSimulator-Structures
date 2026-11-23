"""Bounded source downloads and explicit, view-only mmCIF conversion.

No missing coordinates, force-field parameters, or atomic charges are generated.
An archive model is scientific evidence, not an exact picture of a living nucleus.
"""
from __future__ import annotations

import gzip
import hashlib
import io
import json
import math
import os
import tempfile
from pathlib import Path
from urllib.parse import urlparse
from urllib.request import Request, urlopen

MAX_DOWNLOAD_BYTES = 80_000_000
MAX_SOURCE_BYTES = 256_000_000
DEFAULT_MAX_ATOMS = 250_000
ALLOWED_HOSTS = {"files.rcsb.org", "data.rcsb.org", "alphafold.ebi.ac.uk"}
IMPORTER_VERSION = "nuclear-mmcif-v1"


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def safe_path(root, relative):
    root = Path(root).resolve()
    candidate = (root / relative).resolve()
    if Path(relative).is_absolute() or not candidate.is_relative_to(root):
        raise ValueError(f"Path is outside the structure repository: {relative}")
    return candidate


def atomic_write(path, raw):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    handle, name = tempfile.mkstemp(prefix=".structure-", dir=path.parent)
    try:
        with os.fdopen(handle, "wb") as stream:
            stream.write(raw)
        Path(name).replace(path)
    finally:
        Path(name).unlink(missing_ok=True)


def write_json(path, payload, *, compact=False):
    raw = (json.dumps(payload, ensure_ascii=False, allow_nan=False,
                      indent=None if compact else 2,
                      separators=(",", ":") if compact else None) + "\n").encode()
    if str(path).endswith(".gz"):
        raw = gzip.compress(raw, mtime=0)
    atomic_write(path, raw)


def read_json(path):
    path = Path(path)
    raw = path.read_bytes()
    if path.suffix == ".gz":
        raw = bounded_decompress(raw)
    return json.loads(raw)


def bounded_decompress(raw, max_bytes=MAX_SOURCE_BYTES):
    with gzip.GzipFile(fileobj=io.BytesIO(raw)) as stream:
        result = stream.read(max_bytes + 1)
    if len(result) > max_bytes:
        raise ValueError("Uncompressed source exceeds the configured memory limit")
    return result


def download(url, *, max_bytes=MAX_DOWNLOAD_BYTES):
    parsed = urlparse(url)
    if parsed.scheme != "https" or parsed.hostname not in ALLOWED_HOSTS:
        raise ValueError(f"Unsupported coordinate/metadata host: {url}")
    request = Request(url, headers={"User-Agent": "CellSimulator-Structures/2.0"})
    with urlopen(request, timeout=30) as response:
        final = urlparse(response.geturl())
        if final.scheme != "https" or final.hostname not in ALLOWED_HOSTS:
            raise ValueError("Unexpected source redirect")
        raw = response.read(max_bytes + 1)
    if len(raw) > max_bytes:
        raise ValueError("Download exceeds configured byte limit")
    return raw


def fetch_source(metadata, root):
    """Fetch exactly the pinned uncompressed archive contents into ignored cache."""
    source = metadata["source"]
    ident = metadata["id"]
    if not ident or any(c not in "abcdefghijklmnopqrstuvwxyz0123456789-_" for c in ident):
        raise ValueError("Unsafe asset ID")
    path = Path(root) / "cache" / "nuclear" / f"{ident}.cif"
    expected = source.get("coordinate_sha256")
    if not expected or len(expected) != 64:
        raise ValueError("Asset has no pinned coordinate checksum; resolve it first")
    if path.exists():
        if path.stat().st_size > MAX_SOURCE_BYTES or digest(path.read_bytes()) != expected:
            raise ValueError(f"Cached source differs from pinned model: {ident}")
        return path
    raw = download(source["coordinate_url"])
    if source["coordinate_url"].endswith(".gz"):
        raw = bounded_decompress(raw)
    if digest(raw) != expected:
        raise ValueError(f"Archive content changed for {ident}; review a new version before updating the pin")
    atomic_write(path, raw)
    return path


def require_gemmi():
    try:
        import gemmi
    except ImportError as exc:
        raise RuntimeError("Macromolecule conversion requires: python -m pip install -r requirements-structures.txt") from exc
    return gemmi


def category_rows(block, prefix):
    category = block.get_mmcif_category(prefix)
    if not category:
        return []
    keys = list(category)
    return [dict(zip(keys, row)) for row in zip(*(category[key] for key in keys))]


def read_structure(path, metadata):
    gemmi = require_gemmi()
    raw = Path(path).read_bytes()
    if len(raw) > MAX_SOURCE_BYTES or digest(raw) != metadata["source"]["coordinate_sha256"]:
        raise ValueError("Coordinate input does not match the reviewed source")
    block = gemmi.cif.read_string(raw.decode("utf-8")).sole_block()
    pdb_id = metadata["source"].get("pdb_id")
    if pdb_id:
        entry = block.find_value("_entry.id")
        if not entry or entry.upper() != pdb_id.upper():
            raise ValueError("Coordinate entry ID does not match the asset")
    structure = gemmi.make_structure_from_block(block)
    if not len(structure):
        raise ValueError("Source has no coordinate models")
    return block, structure


def assembly_atom_upper_bound(model, assembly):
    """Count before expansion, so a large pore cannot silently allocate millions of atoms."""
    total = 0
    for generator in assembly.generators:
        subtotal = sum(len(residue) for chain in model for residue in chain
                       if chain.name in generator.chains or residue.subchain in generator.subchains)
        total += subtotal * len(generator.operators)
    return total


def assembly_description(assembly):
    return [{
        "chains": list(g.chains), "label_asym_ids": list(g.subchains),
        "operators": [{"id": op.name, "rotation": [list(row) for row in op.transform.mat.tolist()],
                       "translation_angstrom": list(op.transform.vec.tolist())}
                      for op in g.operators]
    } for g in assembly.generators]


def selected_atoms(residue):
    """One internally consistent residue altloc; common atoms plus the best occupancy label."""
    occupancy = {}
    for atom in residue:
        if atom.altloc not in ("", "\x00", " "):
            occupancy[atom.altloc] = occupancy.get(atom.altloc, 0.0) + atom.occ
    chosen = min(occupancy, key=lambda label: (-occupancy[label], label)) if occupancy else None
    result = []
    seen = set()
    # Common atoms win over alternatives with the same name.
    for atom in sorted(residue, key=lambda a: a.altloc not in ("", "\x00", " ")):
        if atom.altloc not in ("", "\x00", " ", chosen) or atom.name in seen:
            continue
        seen.add(atom.name)
        result.append(atom)
    return result


def convert_mmcif(metadata, path, *, max_atoms=DEFAULT_MAX_ATOMS):
    """Return source atoms and disclosed display connectivity, without MD preparation."""
    from structure_library.identity import asset_component_refs
    if max_atoms < 1:
        raise ValueError("max_atoms must be positive")
    gemmi = require_gemmi()
    block, structure = read_structure(path, metadata)
    model = structure[0]
    assembly_id = metadata.get("assembly", {}).get("id")
    transforms = []
    if assembly_id:
        assembly = next((a for a in structure.assemblies if a.name == assembly_id), None)
        if assembly is None:
            raise ValueError(f"Assembly {assembly_id} is absent from the pinned source")
        upper_bound = assembly_atom_upper_bound(model, assembly)
        if upper_bound > max_atoms:
            raise ValueError(f"Assembly requires up to {upper_bound:,} atoms; budget is {max_atoms:,}. Use chunked rendering or explicitly increase --max-atoms.")
        transforms = assembly_description(assembly)
        model = gemmi.make_assembly(assembly, model, gemmi.HowToNameCopiedChain.AddNumber)
    elif model.count_atom_sites() > max_atoms:
        raise ValueError("Source model exceeds configured atom budget")

    templates = {}
    for row in category_rows(block, "_chem_comp_bond."):
        templates.setdefault(row["comp_id"], []).append(row)
    atoms, coords, bonds, orders = [], [], [], []
    bond_seen = set()
    omitted_waters = omitted_altloc = omitted_polymer_links = 0
    polymer_links = 0
    bfactor_values = []
    def add_bond(i, j, order):
        pair = tuple(sorted((i, j)))
        if i != j and pair not in bond_seen:
            bond_seen.add(pair)
            bonds.append(list(pair))
            orders.append(order)

    for chain in model:
        previous = None
        for residue in chain:
            if residue.is_water():
                omitted_waters += len(residue)
                continue
            chosen = selected_atoms(residue)
            omitted_altloc += len(residue) - len(chosen)
            indices = {}
            for atom in chosen:
                xyz = [atom.pos.x, atom.pos.y, atom.pos.z]
                if not all(math.isfinite(v) for v in xyz) or not math.isfinite(atom.occ):
                    raise ValueError("Source contains nonfinite coordinates or occupancy")
                element = atom.element.name
                if element == "X":
                    raise ValueError("Source contains an unresolved element identity")
                indices[atom.name] = len(atoms)
                record = {
                    "element": element, "atom_name": atom.name,
                    "residue_name": residue.name, "chain_id": chain.name,
                    "residue_index": residue.seqid.num - 1,
                    "insertion_code": residue.seqid.icode.strip(),
                    "label_seq_id": residue.label_seq, "label_asym_id": residue.subchain,
                    "entity_id": residue.entity_id, "source_atom_serial": atom.serial,
                    "occupancy": atom.occ, "source_b_factor": atom.b_iso,
                    "alternate_location": atom.altloc.strip("\x00 "),
                }
                atoms.append(record)
                coords.append(xyz)
                bfactor_values.append(atom.b_iso)
            for bond in templates.get(residue.name, []):
                a, b = indices.get(bond["atom_id_1"]), indices.get(bond["atom_id_2"])
                if a is not None and b is not None:
                    order = {"SING": 1, "DOUB": 2, "TRIP": 3, "AROM": 1.5}.get(bond.get("value_order"))
                    add_bond(a, b, order)
            if previous and residue.label_seq is not None:
                prev_res, prev_indices = previous
                if (prev_res.label_seq is not None and residue.label_seq == prev_res.label_seq + 1
                        and residue.subchain == prev_res.subchain):
                    # Standard polymer links are an explicit chemistry-based display inference.
                    for left, right, max_distance in (("C", "N", 1.8), ("O3'", "P", 2.0)):
                        a, b = prev_indices.get(left), indices.get(right)
                        if a is not None and b is not None:
                            distance = math.dist(coords[a], coords[b])
                            if 1.0 <= distance <= max_distance:
                                add_bond(a, b, 1)
                                polymer_links += 1
                            else:
                                omitted_polymer_links += 1
            previous = (residue, indices) if residue.label_seq is not None else None
    if not atoms or len(atoms) > max_atoms:
        raise ValueError("Empty model or atom budget exceeded")
    predicted = metadata["provenance_class"] == "predicted_model"
    payload = {
        "name": metadata["name"], "asset_id": metadata["id"],
        "coordinate_type": metadata["provenance_class"], "coordinate_units": "angstrom",
        "source_url": metadata["source"]["entry_url"],
        "source": metadata["source"]["database"],
        "source_sha256": metadata["source"]["coordinate_sha256"],
        "view_only": True, "simulation_ready": False,
        "total_charge": 0.0,
        "charge_note": "Loader placeholder only; total molecular charge and partial charges are not calculated.",
        "structure_provenance": metadata,
        "component_refs": asset_component_refs(metadata),
        "processing": {
            "importer": IMPORTER_VERSION, "gemmi_version": gemmi.__version__,
            "model_selection": "first deposited model", "assembly_id": assembly_id,
            "assembly_operators": transforms,
            "coordinate_transform_after_assembly": "identity; no centering, rescaling, relaxation or rounding",
            "omitted_water_atoms": omitted_waters, "omitted_alternate_atoms": omitted_altloc,
            "altloc_policy": "common atoms plus residue-wide highest summed occupancy conformer; alphabetical ties",
            "missing_atoms_added": 0, "hydrogens_added": 0,
            "other_heteroatoms": "retained, including cofactors, modified residues and source-specific ligands",
            "bond_source": "archive chem_comp_bond templates within residues; explicit standard-polymer display inference between consecutive label_seq_id residues",
            "inferred_polymer_display_links": polymer_links,
            "omitted_implausible_polymer_display_links": omitted_polymer_links,
            "connectivity_limit": "Crosslinks and metal coordination are not reconstructed; bond graph is incomplete and not suitable for MD.",
            "scientific_accuracy_percent": None,
        },
        "atoms": atoms, "coords": coords, "bonds": bonds, "bond_orders": orders,
    }
    if predicted:
        payload["prediction_confidence"] = {
            "plddt_per_atom": bfactor_values,
            "note": "AlphaFold local confidence, not experimental temperature factors or an accuracy percentage. Relative-domain uncertainty requires PAE.",
        }
    return payload
