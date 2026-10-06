"""Join the shard outputs of scripts/run_wbm_iem_cross.py --shard K/N into one matrix file.

Usage:
  python scripts/merge_cross_shards.py results/wbm_iem/Harvey_1_03d_iem_cross_ranking_v1_shard*of4.json \
         [--out results/wbm_iem/Harvey_1_03d_iem_cross_ranking_v1.json] \
         [--check-against results/wbm_iem/Harvey_1_03d_iem_results_v0.3.json]

Checks before joining (any failure stops the merge):
  - shards 1..N of one N are all present, once each;
  - the provenance is identical apart from shard, readouts_sha256 and n_readouts;
  - the IEM set-ups agree across shards: same IEMs in the same order, statuses, reactions, healthy pin and v_max
    (within 1e-9 relative);
  - every shard is complete, the shards' readouts are disjoint, and interleaved they rebuild the panel whose
    sha256 the shards recorded.
The merged file has the shape of a single run's output, readouts in panel order, and records the shard files.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SHARD_KEYS = {"shard", "readouts_sha256", "n_readouts"}

spec = importlib.util.spec_from_file_location("run_wbm_iem_cross", os.path.join(ROOT, "scripts", "run_wbm_iem_cross.py"))
X = importlib.util.module_from_spec(spec)
spec.loader.exec_module(X)


def sha256(path):
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def merge(shards):
    """shards: list of (path, records). Returns merged records; raises ValueError on any inconsistency."""
    by_k = {}
    for path, recs in shards:
        if not recs:
            raise ValueError(f"{path}: empty")
        prov = recs[0]["provenance"]
        if "shard" not in prov:
            raise ValueError(f"{path}: not a shard output")
        k, n = (int(x) for x in prov["shard"].split("/"))
        if k in by_k:
            raise ValueError(f"shard {k} given twice")
        by_k[k] = (path, recs, n)
    ns = {n for _, _, n in by_k.values()}
    if len(ns) != 1:
        raise ValueError(f"shards of different N: {sorted(ns)}")
    n = ns.pop()
    if sorted(by_k) != list(range(1, n + 1)):
        raise ValueError(f"need shards 1..{n}; have {sorted(by_k)}")
    base_path, base, _ = by_k[1]
    base_prov = {k: v for k, v in base[0]["provenance"].items() if k not in SHARD_KEYS}
    iems = [(r["iem"], r["call_index"]) for r in base]
    for k, (path, recs, _) in sorted(by_k.items()):
        for r in recs:
            prov = {kk: vv for kk, vv in r["provenance"].items() if kk not in SHARD_KEYS}
            if prov != base_prov:
                diff = sorted(kk for kk in set(prov) | set(base_prov) if prov.get(kk) != base_prov.get(kk))
                raise ValueError(f"{path}: provenance differs from shard 1 in {diff}")
        if [(r["iem"], r["call_index"]) for r in recs] != iems:
            raise ValueError(f"{path}: IEMs differ from shard 1")
        for r, b in zip(recs, base):
            for field in ("iem_reactions", "healthy_pin", "wb_objective_disease_feasible"):
                if r.get(field) != b.get(field):
                    raise ValueError(f"{path}: {r['iem']} {field} differs from shard 1")
            v, w = r.get("vmax_healthy"), b.get("vmax_healthy")
            if (v is None) != (w is None) or (v is not None and abs(v - w) > 1e-9 * max(1.0, abs(w))):
                raise ValueError(f"{path}: {r['iem']} v_max {v} differs from shard 1 ({w})")
            if r["status"] == "readouts_incomplete":
                raise ValueError(f"{path}: {r['iem']} has not finished its readouts")
    # Rebuild the panel by interleaving the shards' readout lists (shard k holds panel[k-1::n]).
    lists = {k: [e["reaction"] for e in recs[0]["readouts"]] for k, (_, recs, _) in by_k.items()}
    total = sum(len(v) for v in lists.values())
    panel = [lists[(i % n) + 1][i // n] for i in range(total)]
    if len(set(panel)) != len(panel):
        raise ValueError("shards overlap")
    if hashlib.sha256(json.dumps(panel).encode()).hexdigest() != base_prov.get("panel_sha256"):
        raise ValueError("interleaved shards do not rebuild the recorded panel")
    merged = []
    for i, b in enumerate(base):
        entries = {}
        for k, (_, recs, _) in by_k.items():
            r = recs[i]
            if [e["reaction"] for e in r["readouts"]] != lists[k] and r["readouts"]:
                raise ValueError(f"shard {k}: {r['iem']} readouts differ from the shard's first IEM")
            for e in r["readouts"]:
                entries[e["reaction"]] = e
        rec = {kk: vv for kk, vv in b.items() if kk not in ("readouts", "check_against_reference", "provenance", "run_fingerprint")}
        rec["readouts"] = [entries[x] for x in panel if x in entries]
        if rec["readouts"]:
            unavailable = sum(1 for e in rec["readouts"] if e["status_healthy"] != "absent" and e["predicted"] == "NA")
            rec["status"] = "partial" if unavailable else "complete"
        merged.append(rec)
    prov = dict(base_prov, readouts_sha256=base_prov["panel_sha256"], n_readouts=len(panel), merged_shards=n,
                merged_from=[{"file": os.path.relpath(os.path.abspath(p), ROOT),
                              "sha256": sha256(p) if os.path.exists(p) else None}
                             for _, (p, _, _) in sorted(by_k.items())])
    fingerprint = hashlib.sha256(json.dumps(prov, sort_keys=True).encode()).hexdigest()
    for rec in merged:
        rec["provenance"], rec["run_fingerprint"] = prov, fingerprint
    return merged, panel


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("shards", nargs="+")
    ap.add_argument("--out")
    ap.add_argument("--check-against")
    args = ap.parse_args()
    shards = []
    for path in args.shards:
        with open(path) as fh:
            shards.append((path, json.load(fh)))
    try:
        merged, panel = merge(shards)
    except ValueError as exc:
        sys.exit(f"not merged: {exc}")
    if args.check_against:
        with open(args.check_against) as fh:
            ref = {(r["iem"], r["call_index"]): r for r in json.load(fh)}
        for rec in merged:
            cmp = X.compare_with_reference(rec, ref)
            if cmp is not None:
                rec["check_against_reference"] = {"file": os.path.relpath(os.path.abspath(args.check_against), ROOT), **cmp}
    out = args.out or args.shards[0].split("_shard")[0] + ".json"
    X.write_json(out, merged)
    cmps = [r["check_against_reference"] for r in merged if "check_against_reference" in r]
    print(f"merged {len(args.shards)} shards: {len(merged)} IEMs x {len(panel)} readouts -> {os.path.relpath(out, ROOT)}; "
          f"statuses {sorted({r['status'] for r in merged})}"
          + (f"; own biomarkers vs reference: {sum(c['n_same_call'] for c in cmps)}/{sum(c['n'] for c in cmps)} same call"
             if cmps else ""))


if __name__ == "__main__":
    main()
