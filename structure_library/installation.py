"""Offline duplicate checks before adding nuclear assets to CellSimulator.

Protein records, assembly representations, and placed molecules are different
things. This module never creates a second Human Proteome protein record, merges
scientific states by name, or removes an existing user file.
"""
from __future__ import annotations

from collections import Counter
from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path
import re
import sqlite3
import subprocess
import tempfile
from urllib.parse import urlparse
from zipfile import BadZipFile, ZipFile

from structure_library.coordinates import read_json, safe_path

INSTALLER = "nuclear-library-v2"
MAX_INDEX_BYTES = 32_000_000


def validate_target(target):
    target = Path(target).expanduser().resolve()
    if not (target / "main.py").is_file() or not (target / "sim/loader.py").is_file():
        raise ValueError("--target must be an existing CellSimulator checkout")
    return target


def _labels(data):
    values = [data.get("name"), data.get("display_name"), data.get("inventory_name")]
    aliases = data.get("aliases", [])
    values.extend([aliases] if isinstance(aliases, str) else aliases)
    return {re.sub(r"\s+", " ", value.strip().casefold())
            for value in values if isinstance(value, str) and value.strip()}


def _descriptors(data):
    identity = data.get("chemical_identity", {})
    found = {}
    for key in ("inchi", "inchi_key", "inchikey"):
        if data.get(key):
            found["inchikey" if key != "inchi" else key] = str(data[key]).strip()
    for item in identity.get("descriptors", []):
        kind = str(item.get("type", "")).casefold()
        if kind in {"inchi", "inchikey"} and item.get("descriptor"):
            found[kind] = str(item["descriptor"]).strip()
    return found


def _pdb(data):
    source = data.get("structure_provenance", {}).get("source", {})
    value = data.get("pdb_id") or data.get("source_pdb_id") or source.get("pdb_id")
    if not value:
        match = re.fullmatch(r"pdb-([a-z0-9]{4})", str(data.get("asset_id", "")))
        value = match[1] if match else ""
    if not value:
        url = urlparse(str(data.get("source_url", "")))
        match = re.fullmatch(r"/download/([A-Za-z0-9]{4})\.(?:pdb|cif)(?:\.gz)?", url.path)
        if url.hostname == "files.rcsb.org" and match:
            value = match[1]
    return str(value).upper()


def _ccd(data):
    value = data.get("pdb_ccd_id") or data.get("chemical_identity", {}).get("ccd_id")
    if not value:
        match = re.fullmatch(r"chemical-ccd-([a-z0-9]+)", str(data.get("id", "")))
        value = match[1] if match else ""
    return str(value).upper()


def _summary(data):
    """Keep one small metadata record per existing item, not its atom arrays."""
    return {
        "id": data.get("asset_id") or data.get("id"),
        "labels": _labels(data), "pdb_id": _pdb(data), "ccd_id": _ccd(data),
        "descriptors": _descriptors(data),
        "source_sha256": data.get("source_sha256"),
        "processing": data.get("processing"),
        "representation_sha256": data.get("representation_sha256"),
        "formal_charge": data.get("formal_charge"),
        "formula": data.get("formula"),
        "installer": data.get("structure_library_installer"),
        "geometry_present": bool(data.get("atoms") and isinstance(data.get("coords"), list)
                                 and len(data["atoms"]) == len(data["coords"])),
    }


class TargetIndex:
    """Scan existing local entries and the packed catalog once per operation.

    The proteome SQLite index is inspected only for the requested accession.
    Large per-protein coordinate caches are never scanned en masse.
    """

    def __init__(self, target):
        self.target = validate_target(target)
        self.records = []
        self.errors = []
        # A sparse/incomplete checkout cannot establish that tracked molecules
        # are absent from the actual application. Inspect only its local index.
        if (self.target / ".git").exists():
            try:
                listing = subprocess.run(["git", "-C", str(self.target), "ls-files", "-z", "--",
                                          "molecules", "macromolecules/builtins"],
                                         capture_output=True, check=True, timeout=5)
                tracked = listing.stdout.decode("utf-8").split("\0")
                missing = [name for name in tracked if name and not (self.target / name).is_file()]
                if missing:
                    self.errors.append(f"Incomplete checkout: {len(missing)} tracked molecule/asset files are missing; "
                                       f"inspect a complete target before installation (for example {missing[0]})")
            except (OSError, UnicodeError, subprocess.SubprocessError) as exc:
                self.errors.append(f"Cannot establish target checkout completeness: {exc}")
        root = self.target / "molecules"
        for path in sorted(root.rglob("*.json")):
            try:
                data = read_json(path)
                if not isinstance(data, dict):
                    raise ValueError("Molecule file must contain an object")
                record = _summary(data)
                record.update(path=path, ref=path.relative_to(root).with_suffix("").as_posix(),
                              location="molecule_file")
                self.records.append(record)
            except (OSError, ValueError, TypeError, AttributeError) as exc:
                self.errors.append(f"{path.relative_to(self.target)}: {exc}")
        # Built-in macromolecules are source metadata, not ordinary molecule JSON.
        for path in sorted((self.target / "macromolecules/builtins").glob("*.json")):
            try:
                data = read_json(path)
                if not isinstance(data, dict) or not _pdb(data):
                    continue
                record = _summary(data)
                record.update(path=path, ref=None, location="builtin_metadata")
                asset = data.get("asset")
                if asset:
                    asset_path = safe_path(self.target / "macromolecules", asset)
                    record["geometry_present"] = False
                    if asset_path.is_file() and data.get("asset_sha256"):
                        record["geometry_present"] = (
                            hashlib.sha256(asset_path.read_bytes()).hexdigest() == data["asset_sha256"])
                record["builtin_id"] = data.get("id")
                self.records.append(record)
            except (OSError, ValueError, TypeError, AttributeError) as exc:
                self.errors.append(f"{path.relative_to(self.target)}: {exc}")
        # Current proteome downloads have tiny provenance sidecars. Read these
        # instead of loading thousands of coordinate arrays for an offline plan.
        cache = self.target / "build/cache/human_proteome"
        for path in sorted(cache.glob("*.meta.json")):
            try:
                data = read_json(path)
                from structure_library.identity import protein_ref
                ref = protein_ref(data.get("uniprot_accession"))
                coordinate = path.with_name(path.name.removesuffix(".meta.json") + ".json")
                record = _summary(data)
                record.update(path=coordinate, ref=ref, location="proteome_cache",
                              geometry_present=coordinate.is_file())
                self.records.append(record)
            except (OSError, ValueError, TypeError, AttributeError) as exc:
                self.errors.append(f"{path.relative_to(self.target)}: {exc}")
        self.legacy_cache_count = sum(not p.with_name(p.stem + ".meta.json").exists()
                                      for p in cache.glob("*.json")
                                      if not p.name.endswith(".meta.json"))
        pack = root / "chebi_conformers.zip"
        if pack.exists():
            try:
                with ZipFile(pack) as archive:
                    info = archive.getinfo("index.json")
                    if info.file_size > MAX_INDEX_BYTES:
                        raise ValueError("Packed catalog index exceeds the metadata budget")
                    index = json.loads(archive.read(info))
                    members = set(archive.namelist())
                if (not isinstance(index, dict) or index.get("schema_version") != 1
                        or not isinstance(index.get("entries"), list)):
                    raise ValueError("Unsupported packed catalog metadata")
                seen_refs = set()
                for entry in index["entries"]:
                    if (not isinstance(entry, dict) or not isinstance(entry.get("ref"), str)
                            or f"{entry['ref']}.json" not in members):
                        raise ValueError("Packed catalog entry is invalid or has no coordinate member")
                    ref = entry["ref"]
                    if (not ref or "\\" in ref or Path(ref).is_absolute()
                            or ".." in Path(ref).parts or ref in seen_refs):
                        raise ValueError("Packed catalog has duplicate or unsafe references")
                    seen_refs.add(ref)
                    record = _summary(entry)
                    record.update(path=pack, ref=entry["ref"], location="packed_molecule",
                                  geometry_present=True)
                    self.records.append(record)
            except (OSError, ValueError, KeyError, TypeError, AttributeError, BadZipFile) as exc:
                self.errors.append(f"molecules/{pack.name}: {exc}")

    def protein_status(self, accession):
        path = self.target / "macromolecules/catalog/human_proteome.sqlite"
        if not path.is_file():
            return "database_missing"
        conn = None
        try:
            conn = sqlite3.connect(path.as_uri() + "?mode=ro", uri=True)
            found = conn.execute("SELECT accession FROM proteins WHERE accession = ?", (accession,)).fetchone()
            return "indexed" if found else "missing_from_index"
        except sqlite3.Error as exc:
            raise ValueError(f"Cannot inspect target Human Proteome database: {exc}") from exc
        finally:
            if conn is not None:
                conn.close()


def representation_digest(payload):
    """Exact content key; unrelated sources or atom selections stay distinct."""
    selected = {key: payload.get(key) for key in (
        "source_sha256", "coordinate_type", "coordinate_units", "processing",
        "atoms", "coords", "bonds", "bond_orders", "bond_details",
        "formal_charge", "total_charge", "prediction_confidence", "chemical_identity",
        "pdb_ccd_id", "pdb_id", "source_pdb_id", "atom_aromatic", "bond_aromatic", "bond_stereo",
    )}
    metadata = payload.get("structure_provenance", {})
    selected["source_selection"] = {key: metadata.get(key) for key in ("source", "assembly", "sequence")}
    return hashlib.sha256(json.dumps(selected, sort_keys=True, separators=(",", ":"),
                                     allow_nan=False).encode()).hexdigest()


def _verified_water(record):
    """Compatibility for the existing NIST water file without changing it."""
    if record["location"] != "molecule_file" or not record["geometry_present"]:
        return False
    data = read_json(record["path"])
    provenance = data.get("provenance", {})
    if (provenance.get("database") != "NIST CCCBDB"
            or provenance.get("source_url") != "https://cccbdb.nist.gov/expgeom2x.asp?casno=7732185"
            or data.get("total_charge") != 0):
        return False
    atoms, bonds = data.get("atoms", []), data.get("bonds", [])
    if Counter(atom.get("element") for atom in atoms) != {"H": 2, "O": 1}:
        return False
    expected = {tuple(sorted((i, j))) for i, atom in enumerate(atoms)
                for j, other in enumerate(atoms) if atom.get("element") == "O" and other.get("element") == "H"}
    try:
        return len(bonds) == 2 and {tuple(sorted(bond)) for bond in bonds} == expected
    except (TypeError, ValueError):
        return False


def _match_result(action, reason, matches, *, exact=False):
    return {
        "action": action, "reason": reason, "exact_representation_match": exact,
        "matches": [{"path": str(row["path"]), "ref": row["ref"],
                     "location": row["location"]} for row in matches],
    }


def plan_payload(payload, target, *, index=None):
    """Conservatively choose reuse, review, or a new assembly/chemical asset."""
    index = index or TargetIndex(target)
    if Path(target).expanduser().resolve() != index.target:
        raise ValueError("TargetIndex belongs to a different CellSimulator checkout")
    target = index.target
    ident = payload.get("asset_id") or payload.get("id")
    if not isinstance(ident, str) or not re.fullmatch(r"[a-z0-9][a-z0-9_-]*", ident):
        raise ValueError("Unsafe or missing asset identifier")
    base = {"asset_id": ident, "target": str(target),
            "component_refs": payload.get("component_refs", [])}
    if index.errors:
        return {**base, "action": "needs_review", "reason": "Existing catalog could not be fully inspected",
                "inspection_errors": index.errors, "matches": [], "exact_representation_match": False}
    metadata = payload.get("structure_provenance", {})
    if metadata.get("provenance_class") == "predicted_model":
        from structure_library.identity import protein_ref
        accession = metadata.get("source", {}).get("uniprot_accession", "")
        if str(metadata.get("organism", {}).get("taxonomy_id")) != "9606":
            return {**base, "action": "needs_review", "reason": "Prediction has no verified human biological source",
                    "matches": [], "exact_representation_match": False}
        try:
            ref = protein_ref(accession)
        except ValueError:
            return {**base, "action": "needs_review", "reason": "Isoform or unresolved protein: the current Human Proteome loader needs a canonical accession",
                    "matches": [], "exact_representation_match": False}
        status = index.protein_status(ref.split("/", 1)[1])
        return {**base, "action": "reuse_protein_identity" if status == "indexed" else "needs_proteome_index",
                "canonical_ref": ref, "catalog_status": status, "matches": [],
                "exact_representation_match": False,
                "reason": "Use the existing Human Proteome entry; no standalone duplicate protein JSON is installed. A catalog entry does not guarantee downloaded geometry."}

    incoming = _summary(payload)
    has_geometry = bool(payload.get("atoms") and payload.get("coords"))
    if has_geometry and len(payload["atoms"]) != len(payload["coords"]):
        return {**base, "action": "needs_review", "reason": "Atom and coordinate counts do not agree",
                "matches": [], "exact_representation_match": False}
    fingerprint = representation_digest(payload) if has_geometry else None
    same_asset = [r for r in index.records if r["id"] == ident]
    for record in same_asset:
        if record["location"] == "molecule_file" and record["geometry_present"] and has_geometry:
            if representation_digest(read_json(record["path"])) == fingerprint:
                return {**base, **_match_result("reuse_exact_asset", "Exact source, selection and geometry already installed", [record], exact=True)}
    if same_asset:
        if not has_geometry and all(r["installer"] == INSTALLER and r["geometry_present"]
                                    and r["source_sha256"] == incoming["source_sha256"]
                                    and r["representation_sha256"] == representation_digest(read_json(r["path"]))
                                    for r in same_asset):
            return {**base, **_match_result("reuse_registered_asset", "Verified registered asset already installed; use export to inspect an alternative representation", same_asset)}
        return {**base, **_match_result("needs_review", "The asset ID already exists with different or unverified contents; no file is overwritten", same_asset)}

    if incoming["pdb_id"]:
        # PDB IDs identify sources, not protein molecules or exact atom selections.
        related = [r for r in index.records if r["pdb_id"] == incoming["pdb_id"]]
        if related:
            unverified = [r for r in related if not r["geometry_present"]]
            if unverified:
                return {**base, **_match_result("needs_review", "Existing source representation has missing or unverified geometry; it is not evidence that the source is absent", unverified)}
            builtin = [r for r in related if r.get("builtin_id") == "nucleosome-2cv5"]
            if builtin and incoming["pdb_id"] == "2CV5":
                return {**base, **_match_result("reuse_builtin_representation", "The nucleus already uses the verified 2CV5 polymer template. The fuller archive conversion is a different representation; export it separately if needed", builtin)}
            return {**base, **_match_result("needs_review", "Same PDB source exists, but assembly/model/atom selection may differ; do not merge or duplicate automatically", related)}
        if index.legacy_cache_count:
            return {**base, "action": "needs_review", "matches": [], "exact_representation_match": False,
                    "reason": f"{index.legacy_cache_count} legacy proteome caches lack source sidecars; their possible PDB overlap must be inspected first"}

    if incoming["ccd_id"]:
        same_component = [r for r in index.records if r["ccd_id"] == incoming["ccd_id"]]
        if same_component:
            conflicts = [r for r in same_component if (
                r["formal_charge"] is not None and incoming["formal_charge"] is not None
                and r["formal_charge"] != incoming["formal_charge"])
                or (r["formula"] and incoming["formula"] and
                    re.sub(r"\s", "", r["formula"]) != re.sub(r"\s", "", incoming["formula"]))
                or any(key in r["descriptors"] and r["descriptors"][key] != value
                       for key, value in incoming["descriptors"].items())]
            if conflicts:
                return {**base, **_match_result("needs_review", "CCD identifier has conflicting chemical state metadata", conflicts)}
            if (incoming["formal_charge"] is None or not incoming["formula"]
                    or any(r["formal_charge"] is None or not r["formula"] for r in same_component)):
                return {**base, **_match_result("needs_review", "CCD label alone does not establish a chemical state; explicit formula and formal charge are required", same_component)}
            if not all(r["geometry_present"] for r in same_component):
                return {**base, **_match_result("needs_review", "Existing CCD entry has no usable geometry", same_component)}
            return {**base, **_match_result("reuse_chemical_identity", "Exact CCD component is already catalogued; retain its existing representation", same_component)}
        if incoming["ccd_id"] == "HOH" and incoming["formal_charge"] == 0:
            water = [r for r in index.records if r["ref"] == "simple_molecules/triatomic/water" and _verified_water(r)]
            if water:
                return {**base, **_match_result("reuse_chemical_identity", "Verified existing NIST water identity; retain its source geometry", water)}
        candidates = [r for r in index.records if (
            incoming["labels"] & r["labels"] or any(value == r["descriptors"].get(key)
            for key, value in incoming["descriptors"].items()))]
        if candidates:
            return {**base, **_match_result("needs_review", "Possible existing chemical: names and standard InChI are candidate links, not proof of identical protonation, tautomer or stereochemical state", candidates)}
    return {**base, "action": "install_new_asset", "matches": [], "exact_representation_match": False,
            "reason": "No existing match found in inspected catalogs; assembly components remain references to protein identities",
            "destination": str(safe_path(target, f"molecules/nucleus/{ident}.json"))}


@contextmanager
def _installation_lock(target):
    """Serialize this importer's writes; never guess that an existing lock is stale."""
    lock = safe_path(target, "build/locks/nuclear-library.lock")
    lock.parent.mkdir(parents=True, exist_ok=True)
    try:
        handle = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except FileExistsError as exc:
        raise ValueError(f"Another nuclear import may be running: {lock}. Check that process before removing a stale lock") from exc
    try:
        with os.fdopen(handle, "w") as stream:
            stream.write(str(os.getpid()) + "\n")
        yield
    finally:
        lock.unlink(missing_ok=True)


def _apply_plan(payload, index, plan):
    action = plan["action"]
    if action.startswith("reuse_"):
        if action == "reuse_protein_identity":
            return plan["canonical_ref"]
        return Path(plan["matches"][0]["path"])
    if action != "install_new_asset":
        raise ValueError(f"{action}: {plan['reason']}")
    if not payload.get("atoms") or len(payload["atoms"]) != len(payload.get("coords", [])):
        raise ValueError("A molecule file must contain coordinates; references are never installed as empty molecules")
    if not (_ccd(payload) or _pdb(payload)):
        raise ValueError("Asset lacks a supported stable chemical or structural source identity")
    dest = Path(plan["destination"])
    if dest.exists():
        raise ValueError(f"Destination appeared after inspection; refusing to overwrite {dest}")
    result = dict(payload, structure_library_installer=INSTALLER,
                  representation_sha256=representation_digest(payload))
    # Atomic creation, never replacement: another installer cannot be overwritten
    # between the existence check and the final write. The temporary file stays
    # on the destination filesystem; hard links are supported on Windows NTFS.
    dest.parent.mkdir(parents=True, exist_ok=True)
    handle, temporary = tempfile.mkstemp(prefix=".structure-", dir=dest.parent)
    try:
        with os.fdopen(handle, "wb") as stream:
            stream.write((json.dumps(result, ensure_ascii=False, allow_nan=False,
                                     separators=(",", ":")) + "\n").encode())
        os.link(temporary, dest)
    finally:
        Path(temporary).unlink(missing_ok=True)
    record = _summary(result)
    record.update(path=dest, ref=dest.relative_to(index.target / "molecules").with_suffix("").as_posix(),
                  location="molecule_file")
    index.records.append(record)
    return dest


def install_payload(payload, target, *, index=None):
    """Install a new supported asset or return its existing location/reference.

    Reinspect while holding the importer lock before any creation. A plan is a
    snapshot, and another import may have added a differently named equivalent.
    """
    index = index or TargetIndex(target)
    plan = plan_payload(payload, target, index=index)
    if plan["action"] != "install_new_asset":
        return _apply_plan(payload, index, plan)
    with _installation_lock(index.target):
        fresh = TargetIndex(index.target)
        current = plan_payload(payload, index.target, index=fresh)
        result = _apply_plan(payload, fresh, current)
        index.records = fresh.records
        index.errors = fresh.errors
        index.legacy_cache_count = fresh.legacy_cache_count
        return result
