"""Paired comparisons between transfer-study arms on the conditions where both wild types grow.

Two gene sets are reported (docs/studies/transfer-method-v1.md, section 5.1):
  union  (primary)   every gene of either arm's model; a gene absent from one arm's model is predicted to have no
                     effect in that arm (knockout growth = that arm's wild-type growth), as the model implies
  common (secondary) genes present in both models (gembench.comparison.aligned_pairs)
with gene-bootstrap 95% intervals and the changed predictions (common genes) and their agreement with the operational
fitness threshold.

Usage: python scripts/transfer_compare.py --role development --orgs Btheta,MR1 --pairs B0:U,U:UN [--out file.json]
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from gembench import metrics as MET  # noqa: E402
from gembench.comparison import aligned_pairs, compare_runs, load_run  # noqa: E402

RESULTS = os.path.join(ROOT, "results", "transfer_v1")


def changes(a, b):
    sa, sb, fit, genes, grows, _ = aligned_pairs(a, b)
    st, ft = a["params"]["growth_threshold"], a["params"]["fitness_threshold"]
    ok = np.isfinite(fit)
    imp = ok & (fit < ft)
    to_ng = ok & (sa >= st) & (sb < st)
    to_g = ok & (sa < st) & (sb >= st)
    return {"to_no_growth": int(to_ng.sum()), "to_no_growth_agree": int((to_ng & imp).sum()),
            "to_growth": int(to_g.sum()), "to_growth_agree": int((to_g & ~imp & ok).sum()),
            "genes_changed": int((to_ng | to_g).any(axis=1).sum())}


def union_pairs(a, b):
    """Secondary analysis: the union of both arms' genes. A gene absent from one arm's model is predicted to have no
    effect there (knockout growth = that arm's wild-type growth). Conditions where both wild types grow."""
    st = a["params"]["growth_threshold"]
    genes = sorted(set(a["browser_genes"]) | set(b["browser_genes"]))
    conds = sorted(set(a["conditions"]) & set(b["conditions"]))
    gi = [{g: i for i, g in enumerate(r["browser_genes"])} for r in (a, b)]
    ci = [{c: i for i, c in enumerate(r["conditions"])} for r in (a, b)]
    grows = [c for c in conds if all(np.isfinite(r["wt_growth"][ix[c]]) and r["wt_growth"][ix[c]] >= st
                                     for r, ix in zip((a, b), ci))]
    sims, fit = [], np.full((len(genes), len(grows)), np.nan)
    for r, g_ix, c_ix in zip((a, b), gi, ci):
        sim = np.zeros((len(genes), len(grows)))
        for j, c in enumerate(grows):
            sim[:, j] = r["wt_growth"][c_ix[c]]
        for i, g in enumerate(genes):
            if g in g_ix:
                row = [c_ix[c] for c in grows]
                sim[i, :] = r["sim_growth"][g_ix[g], row]
                f = r["fitness"][g_ix[g], row]
                known = np.isfinite(fit[i, :])
                if np.any(known & np.isfinite(f) & (fit[i, :] != f)):
                    raise ValueError(f"fitness differs between arms for {g}")
                fit[i, :] = np.where(np.isfinite(f), f, fit[i, :])
        sims.append(sim)
    return sims[0], sims[1], fit, genes, grows


def gene_counts(sim, fit, st, ft):
    """Per-gene confusion counts [tp, tn, fp, fn] over finite observations ('positive' = no growth defect, as in
    gembench.metrics.mcc: observed fitness >= threshold, predicted knockout growth >= threshold)."""
    ok = np.isfinite(sim) & np.isfinite(fit)
    obs = np.where(ok, fit >= ft, False)
    pred = np.where(ok, sim >= st, False)
    return np.stack([(ok & obs & pred).sum(1), (ok & ~obs & ~pred).sum(1),
                     (ok & ~obs & pred).sum(1), (ok & obs & ~pred).sum(1)], axis=1).astype(float)


def mcc_from_counts(c):
    tp, tn, fp, fn = c
    if tp + tn + fp + fn == 0:
        return float("nan")
    den = (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)
    return float((tp * tn - fp * fn) / np.sqrt(den)) if den > 0 else 0.0


def paired_bootstrap(sa, sb, fit, st, ft, n_boot=1000, seed=0):
    """Gene-resampling bootstrap of MCC(B) - MCC(A), drawing indices exactly as gembench.comparison.compare_runs does;
    MCC depends only on confusion totals, so each resample sums per-gene counts (fast, identical values)."""
    ca, cb = gene_counts(sa, fit, st, ft), gene_counts(sb, fit, st, ft)
    n = ca.shape[0]
    rng = np.random.default_rng(seed)
    d = []
    for _ in range(n_boot if n else 0):
        w = np.bincount(rng.integers(0, n, n), minlength=n)
        d.append(mcc_from_counts(w @ cb) - mcc_from_counts(w @ ca))
    d = np.asarray(d)
    d = d[np.isfinite(d)]
    point = mcc_from_counts(cb.sum(0)) - mcc_from_counts(ca.sum(0))
    return mcc_from_counts(ca.sum(0)), mcc_from_counts(cb.sum(0)), point, \
        (np.percentile(d, [2.5, 97.5]).tolist() if len(d) else [float("nan")] * 2)


def union_delta(a, b, n_boot=1000, seed=0):
    sa, sb, fit, genes, grows = union_pairs(a, b)
    st, ft = a["params"]["growth_threshold"], a["params"]["fitness_threshold"]
    ma, mb, d, ci = paired_bootstrap(sa, sb, fit, st, ft, n_boot, seed)
    return {"union_genes": len(genes), "union_mcc_A": ma, "union_mcc_B": mb, "union_delta": d, "union_delta_ci95": ci}


def common_delta(a, b, n_boot=1000, seed=0):
    """Common-gene comparison: point values from gembench.comparison.compare_runs (no bootstrap there), interval from
    the same index stream summed per gene."""
    c = compare_runs(a, b, n_boot=0, seed=seed)
    sa, sb, fit, genes, _, _ = aligned_pairs(a, b)
    st, ft = a["params"]["growth_threshold"], a["params"]["fitness_threshold"]
    ma, mb, d, ci = paired_bootstrap(sa, sb, fit, st, ft, n_boot, seed)
    if not (np.isclose(ma, c["mcc"]["A"], equal_nan=True) and np.isclose(mb, c["mcc"]["B"], equal_nan=True)):
        raise AssertionError(f"count-based MCC {ma},{mb} differs from gembench.metrics {c['mcc']}")
    c["paired_mcc_difference_B_minus_A"]["ci95"] = ci
    c["paired_mcc_difference_B_minus_A"]["finite_bootstrap_samples"] = n_boot
    c["paired_mcc_difference_B_minus_A"]["requested_bootstrap_samples"] = n_boot
    return c


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--role", default="development")
    ap.add_argument("--orgs", required=True)
    ap.add_argument("--pairs", required=True)
    ap.add_argument("--root", default=RESULTS)
    ap.add_argument("--n-boot", type=int, default=1000)
    ap.add_argument("--out")
    a = ap.parse_args()
    rows = []
    for org in a.orgs.split(","):
        for pair in a.pairs.split(","):
            x, y = pair.split(":")
            da, db = (os.path.join(a.root, a.role, org, n) for n in (x, y))
            if not (os.path.exists(os.path.join(da, "card.json")) and os.path.exists(os.path.join(db, "card.json"))):
                rows.append({"org": org, "A": x, "B": y, "missing": True})
                continue
            ra, rb = load_run(da), load_run(db)
            c = common_delta(ra, rb, n_boot=a.n_boot, seed=0)
            ch = changes(ra, rb)
            un = union_delta(ra, rb, n_boot=a.n_boot)
            row = {"org": org, "A": x, "B": y, "wt_grows_A": c["wt_grows"]["A"], "wt_grows_B": c["wt_grows"]["B"],
                   "n_mapped_conditions": len(ra["conditions"]), "n_conditions_both_grow": c["n_conditions_both_grow"],
                   "n_pairs": c["n_shared_finite_pairs"], "mcc_A": c["mcc"]["A"], "mcc_B": c["mcc"]["B"],
                   "delta": c["paired_mcc_difference_B_minus_A"]["point"], "delta_ci95": c["paired_mcc_difference_B_minus_A"]["ci95"],
                   **ch, **un}
            rows.append(row)
            print(f"{org:8s} {x:>10s} -> {y:<10s} WT {row['wt_grows_A']:>2}->{row['wt_grows_B']:<2}/{row['n_mapped_conditions']:<2} "
                  f"pairs={row['n_pairs']:>6} MCC {row['mcc_A']:.4f}->{row['mcc_B']:.4f} d={row['delta']:+.4f} "
                  f"[{row['delta_ci95'][0]:+.4f},{row['delta_ci95'][1]:+.4f}] to-no-growth {ch['to_no_growth']} ({ch['to_no_growth_agree']} agree) "
                  f"to-growth {ch['to_growth']} ({ch['to_growth_agree']} agree)" +
                  f" | union genes {un['union_genes']}: MCC {un['union_mcc_A']:.4f}->{un['union_mcc_B']:.4f} d={un['union_delta']:+.4f} "
                  f"[{un['union_delta_ci95'][0]:+.4f},{un['union_delta_ci95'][1]:+.4f}]",
                  flush=True)
    if a.out:
        json.dump(rows, open(a.out, "w"), indent=1)


if __name__ == "__main__":
    main()
