"""Matched, retrospective comparisons; independent validation must be designed separately."""
from __future__ import annotations

import json
from pathlib import Path
import numpy as np

from . import metrics as M


def load_run(directory):
    directory = Path(directory)
    with np.load(directory / "matrices.npz", allow_pickle=False) as z:
        run = {k: z[k].copy() for k in ["sim_growth", "fitness", "wt_growth", "browser_genes", "conditions"]}
    card = json.loads((directory / "card.json").read_text())
    run["params"] = card["protocol"]["params"]
    run["benchmark"] = card["benchmark"]
    for key in ("browser_genes", "conditions"):
        if len(set(run[key])) != len(run[key]):
            raise ValueError(f"Duplicate {key}: comparison would be ambiguous")
    shape = (len(run["browser_genes"]), len(run["conditions"]))
    if run["sim_growth"].shape != shape or run["fitness"].shape != shape or run["wt_growth"].shape != (shape[1],):
        raise ValueError("Run matrices do not match their gene/condition labels")
    return run


def aligned_pairs(a, b):
    if a["benchmark"] != b["benchmark"]:
        raise ValueError("Cannot compare different benchmark organisms")
    for key in ("growth_threshold", "fitness_threshold"):
        if a["params"][key] != b["params"][key]:
            raise ValueError(f"Different {key}; rerun with a common protocol")
    st = a["params"]["growth_threshold"]
    genes = sorted(set(a["browser_genes"]) & set(b["browser_genes"]))
    conditions = sorted(set(a["conditions"]) & set(b["conditions"]))
    gene_maps = [{g: i for i, g in enumerate(r["browser_genes"])} for r in (a, b)]
    cond_maps = [{c: i for i, c in enumerate(r["conditions"])} for r in (a, b)]
    grows = [c for c in conditions if all(np.isfinite(r["wt_growth"][ix[c]]) and
             r["wt_growth"][ix[c]] >= st for r, ix in zip((a, b), cond_maps))]
    arrays = []
    for r, gi, ci in zip((a, b), gene_maps, cond_maps):
        index = np.ix_([gi[g] for g in genes], [ci[c] for c in grows])
        arrays.append((r["sim_growth"][index], r["fitness"][index]))
    (sa, fa), (sb, fb) = arrays
    if not np.array_equal(fa, fb, equal_nan=True):
        raise ValueError("Experimental fitness differs on aligned observations; inputs cannot be pooled")
    valid = np.isfinite(sa) & np.isfinite(sb) & np.isfinite(fa)
    fit = np.where(valid, fa, np.nan)
    return sa, sb, fit, genes, grows, conditions


def compare_runs(a, b, *, n_boot=500, seed=0):
    sa, sb, fit, genes, grows, common_conditions = aligned_pairs(a, b)
    st, ft = a["params"]["growth_threshold"], a["params"]["fitness_threshold"]
    fn = lambda s, f: M.mcc(s, f, st, ft)
    ma, mb = fn(sa, fit), fn(sb, fit)
    rng = np.random.default_rng(seed)
    differences = []
    for _ in range(n_boot if len(genes) else 0):
        ix = rng.integers(0, len(genes), len(genes))
        differences.append(fn(sb[ix], fit[ix]) - fn(sa[ix], fit[ix]))
    values = np.asarray(differences)
    values = values[np.isfinite(values)]
    ci = np.percentile(values, [2.5, 97.5]).tolist() if len(values) else [float("nan"), float("nan")]
    return {
        "evaluation_role": "retrospective_development",
        "n_common_genes": len(genes), "n_common_conditions": len(common_conditions),
        "n_conditions_both_grow": len(grows), "conditions_both_grow": grows,
        "n_shared_finite_pairs": int(np.isfinite(fit).sum()),
        "wt_grows": {label: int((np.isfinite(r["wt_growth"]) & (r["wt_growth"] >= st)).sum())
                     for label, r in [("A", a), ("B", b)]},
        "mcc": {"A": ma, "B": mb},
        "confusion": {"A": M.confusion(sa, fit, st, ft), "B": M.confusion(sb, fit, st, ft)},
        "paired_mcc_difference_B_minus_A": {"point": mb-ma, "ci95": ci,
            "finite_bootstrap_samples": len(values), "requested_bootstrap_samples": n_boot},
        "uncertainty_note": "Paired percentile bootstrap over genes; conditions are fixed. This does not account for adaptive patch selection or establish held-out performance.",
    }
