"""Search the nuclear collection and explicitly fetch/install individual assets.

Listing never downloads coordinates. The inventory is broader than the set of
reviewed structures and is not a recipe for the population of one G1 nucleus.
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
import sys
from pathlib import Path

from structure_library.coordinates import (DEFAULT_MAX_ATOMS, convert_mmcif,
    digest, fetch_source, read_json, safe_path, write_json)
from structure_library.identity import (asset_component_refs, collection_protein_links,
    reconcile_inventory)
from structure_library.installation import TargetIndex, install_payload, plan_payload

ROOT = Path(__file__).resolve().parent
COLLECTION = ROOT / "collections" / "human-nucleus"


def inventory(root=ROOT):
    return read_json(Path(root) / "collections/human-nucleus/inventory.json")


def assets(root=ROOT):
    result = []
    for filename in ("chemical-assets.json", "macromolecular-assets.json"):
        path = Path(root) / "collections" / "human-nucleus" / filename
        if path.exists():
            result.extend(read_json(path)["assets"])
    # Only predictions explicitly selected by the user are discoverable here.
    # Scanning local pins performs no network calls or bulk model resolution.
    for path in sorted((Path(root) / "cache/predictions").glob("*/metadata.json")):
        meta = read_json(path)
        if meta.get("provenance_class") == "predicted_model":
            result.append({"id": meta["id"], "name": meta["name"],
                           "kind": "prediction", "provenance_class": "predicted_model",
                           "metadata_path": path.relative_to(root).as_posix(),
                           "uniprot_accession": meta["source"]["uniprot_accession"],
                           "source_status": "explicitly_selected_cached_prediction"})
    ids = [entry["id"] for entry in result]
    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate structure asset IDs")
    return result


def matches(entry, query):
    q = query.casefold()
    return any(q in str(entry.get(key, "")).casefold()
               for key in ("id", "name", "gene_symbol", "uniprot_accessions"))


def find_entry(query, root=ROOT):
    entries = inventory(root)["entries"]
    exact = [entry for entry in entries if query.casefold() in {
        entry["id"].casefold(), entry.get("gene_symbol", "").casefold(),
        *[accession.casefold() for accession in entry.get("uniprot_accessions", [])],
    }]
    hits = exact or [entry for entry in entries if matches(entry, query)]
    if len(hits) != 1:
        raise ValueError(f"Expected one inventory entry, found {len(hits)}. Use an exact inventory ID from 'list --query'.")
    return hits[0]


def linked_assets(entry, available):
    """References are never promoted to claims of complete protein coverage."""
    if entry["kind"] != "protein_gene":
        return [asset["id"] for asset in available if asset.get("inventory_name") == entry["name"]]
    candidate_pdbs = {pdb.lower() for source in entry.get("source_lookups", []) for pdb in source.get("pdb_ids", [])}
    return [asset["id"] for asset in available
            if (asset["id"].startswith("pdb-") and asset["id"][4:] in candidate_pdbs)
            or (asset.get("kind") == "prediction" and
                asset.get("uniprot_accession") in entry.get("uniprot_accessions", []))]


def read_chemical(asset, root=ROOT):
    meta = read_json(safe_path(root, asset["metadata_path"]))
    source = safe_path(root, meta["source"]["source_path"]).read_bytes()
    raw = safe_path(root, asset["coordinate_path"]).read_bytes()
    if digest(source) != asset["source_sha256"] or digest(raw) != asset["coordinate_sha256"]:
        raise ValueError(f"Bundled chemical asset checksum mismatch: {asset['id']}")
    data = json.loads(raw)
    if data.get("id") != asset["id"] or data.get("simulation_ready") is not False:
        raise ValueError("Unexpected chemical identity or simulation claim")
    # Preserve the pinned source bytes; enrich the exported/installed view with
    # the identity already recorded in the CCD metadata, never a name-based guess.
    return dict(data, pdb_ccd_id=meta["source"]["ccd_id"],
                chemical_identity={"ccd_id": meta["source"]["ccd_id"],
                                   **meta["chemical_identity"]})


def get_asset(query, root=ROOT):
    found = [asset for asset in assets(root) if asset["id"] == query]
    if len(found) != 1:
        raise ValueError(f"No selected coordinate asset {query!r}. Use 'assets' for usable models or 'resolve' for inventory candidates.")
    return found[0]


def fetch_asset(asset, root=ROOT):
    if asset.get("source_kind") == "wwpdb_ccd":
        read_chemical(asset, root)
        return safe_path(root, asset["coordinate_path"])
    metadata = read_json(safe_path(root, asset["metadata_path"]))
    if asset.get("kind") == "prediction":
        from structure_library.predictions import _cached_metadata
        metadata = _cached_metadata(safe_path(root, asset["metadata_path"]), Path(root),
                                    asset["uniprot_accession"], metadata["source"]["model_id"])
    return fetch_source(metadata, root)


def payload_for_asset(asset, root=ROOT, *, max_atoms=DEFAULT_MAX_ATOMS):
    if asset.get("source_kind") == "wwpdb_ccd":
        return read_chemical(asset, root)
    metadata = read_json(safe_path(root, asset["metadata_path"]))
    path = fetch_asset(asset, root)
    return convert_mmcif(metadata, path, max_atoms=max_atoms)


def plan_asset(asset, target, root=ROOT, *, index=None, max_atoms=DEFAULT_MAX_ATOMS):
    """Make the installation decision using local metadata before downloading."""
    if asset.get("source_kind") == "wwpdb_ccd":
        payload = read_chemical(asset, root)
    else:
        metadata = read_json(safe_path(root, asset["metadata_path"]))
        payload = {"asset_id": asset["id"], "name": asset["name"],
                   "source_sha256": metadata["source"]["coordinate_sha256"],
                   "structure_provenance": metadata,
                   "component_refs": asset_component_refs(metadata)}
    plan = plan_payload(payload, target, index=index)
    upper_bound = asset.get("atom_count_upper_bound", 0)
    if plan["action"] == "install_new_asset" and upper_bound > max_atoms:
        plan.update(action="needs_chunking", reason=(
            f"Selected assembly upper bound {upper_bound:,} exceeds the {max_atoms:,}-atom budget; no source is downloaded"))
    return plan


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("summary", help="Show counts and evidence limitations; no network")
    listing = sub.add_parser("list", help="Search all inventory identities; no coordinate downloads")
    listing.add_argument("--query", default="")
    listing.add_argument("--limit", type=int, default=30)
    listing.add_argument("--all", action="store_true", help="Print every matching identity")
    listing.add_argument("--json", action="store_true")
    sub.add_parser("assets", help="List bundled chemical geometry and pinned macromolecular sources")
    reconcile = sub.add_parser("reconcile", help="Join nuclear gene evidence to the target Human Proteome; read-only, offline")
    reconcile.add_argument("--target", type=Path, required=True)
    reconcile.add_argument("--json", action="store_true", help="Include all canonical references and unresolved mappings")
    plan = sub.add_parser("plan", help="Check reuse and possible duplicates before downloading or installing; offline")
    plan.add_argument("which", help="Asset ID or 'all'")
    plan.add_argument("--target", type=Path, required=True)
    plan.add_argument("--max-atoms", type=int, default=DEFAULT_MAX_ATOMS)
    plan.add_argument("--json", action="store_true")
    resolve = sub.add_parser("resolve", help="Show exact source candidates; optional explicit AlphaFold resolution")
    resolve.add_argument("which")
    resolve.add_argument("--predicted", action="store_true", help="Resolve and cache a labeled AlphaFold model")
    resolve.add_argument("--accession")
    resolve.add_argument("--model-id")
    fetch = sub.add_parser("fetch", help="Cache pinned sources, without modifying CellSimulator")
    fetch.add_argument("which", help="Asset ID or 'all' (curated assets only; never every gene prediction)")
    export = sub.add_parser("export", help="Write a view-only asset JSON, without changing the application")
    export.add_argument("which")
    export.add_argument("--output", type=Path, required=True)
    export.add_argument("--max-atoms", type=int, default=DEFAULT_MAX_ATOMS)
    install = sub.add_parser("install", help="Explicitly install view-only models into CellSimulator")
    install.add_argument("which", help="Asset ID or 'all'; unresolved identities have no fabricated model")
    install.add_argument("--target", type=Path, required=True)
    install.add_argument("--max-atoms", type=int, default=DEFAULT_MAX_ATOMS)
    args = parser.parse_args(argv)
    try:
        available = assets()
        if args.command == "summary":
            data = inventory()
            chemical = [asset for asset in available if asset.get("source_kind") == "wwpdb_ccd"]
            macro = [asset for asset in available if asset.get("kind") == "macromolecule"]
            print(json.dumps({
                "inventory_counts": data["counts"],
                "unique_verified_protein_references": len(collection_protein_links(data["entries"])),
                "protein_identity_policy": "Reuse human-proteome/<UniProt accession>; source gene rows are evidence, not additional protein definitions",
                "bundled_chemical_coordinate_assets": len(chemical),
                "pinned_macromolecular_sources": len(macro),
                "explicitly_selected_predictions_in_local_cache": sum(asset.get("kind") == "prediction" for asset in available),
                "complete_inventory": False,
                "all_identities_have_reviewed_coordinates": False,
                "default_g1_populations_assigned": False,
                "notes": ["Source cross-references are candidates, not selected full-length structures.",
                          "HeLa protein-group quantities remain separate from default G1 populations.",
                          "Chemical states, structures, local quality and nuclear placement require separate evidence."],
            }, indent=2))
        elif args.command == "list":
            entries = [entry for entry in inventory()["entries"] if matches(entry, args.query)]
            shown = entries if args.all else entries[:max(1, args.limit)]
            rows = [{"id": entry["id"], "name": entry["name"], "kind": entry["kind"],
                     "gene_symbol": entry.get("gene_symbol"),
                     "canonical_protein_refs": [link["canonical_ref"] for link in collection_protein_links([entry])],
                     "representative_or_component_assets": linked_assets(entry, available),
                     "default_g1_copies": None} for entry in shown]
            if args.json:
                print(json.dumps({"total_matches": len(entries), "entries": rows}, indent=2))
            else:
                for row in rows:
                    linked = ",".join(row["representative_or_component_assets"]) or "no selected asset"
                    print(f"{row['id']} | {row['gene_symbol'] or ''} {row['name']} | {linked}")
                print(f"Showing {len(shown)} of {len(entries)} matches. A linked assembly may cover only a component or fragment.")
        elif args.command == "assets":
            for asset in available:
                print(f"{asset['id']} | {asset['name']} | {asset['provenance_class']}")
        elif args.command == "reconcile":
            report = reconcile_inventory(inventory()["entries"], args.target)
            if not args.json:
                report = {key: value for key, value in report.items()
                          if key not in {"protein_links", "unresolved_mappings"}}
            print(json.dumps(report, indent=2))
            if report["proteome_database"]["status"] == "unreadable":
                return 1
        elif args.command == "plan":
            index = TargetIndex(args.target)
            selected = available if args.which == "all" else [get_asset(args.which)]
            plans = [plan_asset(asset, args.target, index=index, max_atoms=args.max_atoms)
                     for asset in selected]
            if args.json:
                print(json.dumps({"actions": dict(Counter(p["action"] for p in plans)),
                                  "assets": plans, "network_requests": 0, "files_written": 0}, indent=2))
            else:
                for result in plans:
                    print(f"{result['asset_id']} | {result['action']} | {result['reason']}")
                print("Read-only plan. No downloads, installed files, protein records or scene populations changed.")
        elif args.command == "resolve":
            entry = find_entry(args.which)
            if args.predicted:
                from structure_library.predictions import resolve_prediction
                result = resolve_prediction(entry, ROOT, accession=args.accession, model_id=args.model_id)
            else:
                result = {"inventory_entry": entry, "representative_or_component_assets": linked_assets(entry, available),
                          "note": "Candidate links are not automatically chosen as complete molecular models."}
            print(json.dumps(result, indent=2))
        elif args.command == "fetch":
            selected = available if args.which == "all" else [get_asset(args.which)]
            failed = []
            for asset in selected:
                try:
                    print(f"{asset['id']}: {fetch_asset(asset)}", flush=True)
                except (OSError, ValueError, RuntimeError) as exc:
                    failed.append(asset["id"])
                    print(f"UNAVAILABLE {asset['id']}: {exc}", file=sys.stderr, flush=True)
            if failed:
                print(f"Failed assets: {', '.join(failed)}", file=sys.stderr)
                return 1
        elif args.command == "export":
            payload = payload_for_asset(get_asset(args.which), max_atoms=args.max_atoms)
            write_json(args.output, payload, compact=True)
            print(f"Exported {len(payload['atoms']):,} source/model atoms to {args.output}")
        elif args.command == "install":
            selected = available if args.which == "all" else [get_asset(args.which)]
            failed = []
            index = TargetIndex(args.target)
            for asset in selected:
                try:
                    decision = plan_asset(asset, args.target, index=index, max_atoms=args.max_atoms)
                    if decision["action"].startswith("reuse_"):
                        destination = decision.get("canonical_ref") or decision["matches"][0]["path"]
                        print(f"Reused {asset['id']}: {destination}. {decision['reason']}", flush=True)
                        continue
                    if decision["action"] != "install_new_asset":
                        raise ValueError(f"{decision['action']}: {decision['reason']}")
                    payload = payload_for_asset(asset, max_atoms=args.max_atoms)
                    print(f"Available {asset['id']}: {install_payload(payload, args.target, index=index)}", flush=True)
                except (OSError, ValueError, RuntimeError) as exc:
                    failed.append(asset["id"])
                    print(f"NOT INSTALLED {asset['id']}: {exc}", file=sys.stderr, flush=True)
            if failed:
                print(f"Not installed: {', '.join(failed)}. Large assemblies may require chunked viewing.", file=sys.stderr)
                return 1
        return 0
    except (OSError, ValueError, RuntimeError, KeyError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
