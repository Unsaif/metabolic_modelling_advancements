"""Independent bounds checks for Harvetta 1.03d (v0.4) and, for comparison, Harvey 1.03d.

1. Load the model file (scipy.io) and MATLAB's setup bounds (step 1) and global bounds (step 2); check reaction order.
2. sha256 of the MATLAB global bounds (lb bytes + ub bytes, little-endian float64), raw and with -0.0 normalised,
   against the lp_bounds_sha256 recorded by the Python runs (Harvetta run A and C; Harvey v0.3 and run B).
3. My own transcription of runIEM_HH's 'unified reaction constraints' block, applied to MATLAB's setup bounds,
   compared bit for bit with MATLAB's global bounds; counts against the provenance (set_irreversible/closed/bile duct).
4. Every difference between MATLAB's setup bounds and the bounds stored in the model, grouped by reaction family,
   bound and ratio (this is what re-applying the current Toolbox code changes); implied older GFR = 90 / kidney ratio.
5. Whether any setup change sits on an entry that the global step overwrites.
Imports nothing from gembench or scripts/. Writes bounds_v04.json here.
"""
import os
import sys
from collections import Counter, OrderedDict

import numpy as np

from vcommon import ROOT, dump, load_bounds, load_json, load_model_bounds, sha_bounds

MODELS = OrderedDict([
    ("Harvetta", dict(model="external/COBRA.models/mat/Harvetta_1_03d.mat",
                      setup="results/wbm_iem/matlab_reference/harvetta_step1/matlab_setup_bounds_Harvetta_1_03d.mat",
                      glob="results/wbm_iem/matlab_reference/harvetta_step2/matlab_global_bounds_Harvetta_1_03d.mat",
                      runs={"A (max)": "results/wbm_iem/Harvetta_1_03d_iem_results_v0.4.json",
                            "C (min)": "results/wbm_iem/Harvetta_1_03d_iem_results_v0.4_min.json"})),
    ("Harvey", dict(model="external/COBRA.models/mat/Harvey_1_03d.mat",
                    setup="results/wbm_iem/matlab_reference/step1/matlab_setup_bounds_Harvey_1_03d.mat",
                    glob="results/wbm_iem/matlab_reference/step2/matlab_global_bounds_Harvey_1_03d.mat",
                    runs={"v0.3 (max)": "results/wbm_iem/Harvey_1_03d_iem_results_v0.3.json",
                          "B (min)": "results/wbm_iem/Harvey_1_03d_iem_results_v0.4_min.json"})),
])

# Transcribed by me from runIEM_HH_harvetta_ref.m (unified reaction constraints, lines ~117-160 of the stripped file).
R_IRR = ['_ARGSL', '_GACMTRc', '_FUM', '_FUMm', '_HMR_7698', '_UAG4E', '_UDPG4E', '_GALT', '_G6PDH2c', '_G6PDH2r',
         '_G6PDH2rer', '_GLUTCOADHm', '_r0541', '_ACOAD8m', '_RE2410C', '_RE2410N']
R_EXCL = ['_FUMt', '_FUMAC', '_FUMS', 'BBB']
R_CLOSE = ['_r0784', '_r0463']
BILE = ['12dhchol', '3dhcdchol', '3dhchol', '3dhdchol', '3dhlchol', '7dhcdchol', '7dhchol', 'cdca24g', 'cdca3g', 'cholate',
        'dca24g', 'dca3g', 'dchac', 'dgchol', 'gchola', 'hca24g', 'hca6g', 'hdca24g', 'hdca6g', 'hyochol', 'icdchol',
        'isochol', 'lca24g', 'tchola', 'tdchola', 'tdechola', 'thyochol', 'uchol']
RNEW = {f"BileDuct_EX_{b}[bd]_[luSI]" for b in BILE}


def global_step(rxns, lb, ub):
    lb, ub = lb.copy(), ub.copy()
    mic = {r for r in rxns if "Micro_" in r}
    irr = ({r for r in rxns if any(p in r for p in R_IRR)} - {r for r in rxns if any(p in r for p in R_EXCL)}) - mic
    m_irr = np.array([r in irr for r in rxns])
    lb[m_irr] = 0.0
    close = {r for r in rxns if any(p in r for p in R_CLOSE)} - mic
    m_close = np.array([r in close for r in rxns])
    lb[m_close] = 0.0
    ub[m_close] = 0.0
    m_bile = np.array([r in RNEW for r in rxns])
    ub[m_bile] = 100.0
    return lb, ub, dict(set_irreversible=int(m_irr.sum()), closed=int(m_close.sum()), bile_duct_ub_100=int(m_bile.sum()))


def bit_equal(a, b):
    return np.ascontiguousarray(a, dtype="<f8").view(np.uint64) == np.ascontiguousarray(b, dtype="<f8").view(np.uint64)


def norm0(a):
    a = a.copy()
    a[a == 0] = 0.0
    return a


def family(r, bound, old, new):
    if r.startswith("Kidney_EX_"):
        return "kidney"
    if r.startswith("BBB_") and r.endswith("[CSF]exp"):
        return "csf_export"
    if r.startswith("BBB_") and r.endswith("[CSF]upt"):
        return "bbb_uptake"
    if r.startswith("Diet_EX_"):
        if bound == "lb" and old == 0 and new == -0.1:
            return "diet_opened_0.1"
        if old == 0 and ((bound == "lb" and new == -12.0) or (bound == "ub" and new == -8.0)):
            return "diet_10pm20pct"
        return "diet_other"
    return "other"


def classify(rxns, s_lb, s_ub, n_lb, n_ub):
    groups = OrderedDict()
    for bound, s, n in (("lb", s_lb, n_lb), ("ub", s_ub, n_ub)):
        for i in np.where(~bit_equal(s, n))[0]:
            r, a, b = rxns[i], float(s[i]), float(n[i])
            fam = family(r, bound, a, b)
            ratio = round(b / a, 4) if a != 0 else None
            g = groups.setdefault(fam, dict(n=0, by_bound=Counter(), ratios=Counter(), examples=[], exact_ratios=[]))
            g["n"] += 1
            g["by_bound"][bound] += 1
            g["ratios"][str(ratio)] += 1
            if a != 0:
                g["exact_ratios"].append(b / a)
            if len(g["examples"]) < 4:
                g["examples"].append((r, bound, a, b))
    for g in groups.values():
        ex = g.pop("exact_ratios")
        g["by_bound"] = dict(g["by_bound"])
        g["ratios"] = dict(g["ratios"])
        if ex:
            g["ratio_min"], g["ratio_max"] = min(ex), max(ex)
    return groups


def one(model, cfg):
    rep = OrderedDict()
    var, rx, s_lb, s_ub, m = load_model_bounds(os.path.join(ROOT, cfg["model"]))
    rx1, lb1, ub1 = load_bounds(os.path.join(ROOT, cfg["setup"]))
    rx2, lb2, ub2 = load_bounds(os.path.join(ROOT, cfg["glob"]))
    rep["model_variable"] = var
    rep["n_rxns"] = dict(model=len(rx), setup=len(rx1), glob=len(rx2), unique=len(set(rx)))
    rep["same_order"] = dict(setup=rx1 == rx, glob=rx2 == rx)
    rep["n_mets"] = int(np.asarray(m["mets"]).size)
    fields = list(m.dtype.names)
    rep["has_coupling"] = {k: (k in fields) for k in ("C", "d", "dsense", "ctrs")}
    if "C" in fields:
        rep["coupling_shape"] = list(np.asarray(m["C"]).shape) if not hasattr(m["C"], "shape") else list(m["C"].shape)
    rep["n_bileduct_exits"] = sum(r.startswith("BileDuct_EX_") for r in rx)
    rep["recorded_lp_bounds_sha256"] = {}
    for k, p in cfg["runs"].items():
        hs = {r["provenance"]["lp_bounds_sha256"] for r in load_json(p)}
        rep["recorded_lp_bounds_sha256"][k] = sorted(hs)
    h_raw, h_norm = sha_bounds(lb2, ub2), sha_bounds(norm0(lb2), norm0(ub2))
    rep["matlab_global_sha256"] = dict(raw=h_raw, neg_zero_normalised=h_norm,
                                       raw_matches=[k for k, v in rep["recorded_lp_bounds_sha256"].items() if v == [h_raw]],
                                       normalised_matches=[k for k, v in rep["recorded_lp_bounds_sha256"].items() if v == [h_norm]])
    rep["matlab_setup_sha256"] = sha_bounds(lb1, ub1)
    lbg, ubg, cnt = global_step(rx, lb1, ub1)
    rep["my_global_step_counts"] = cnt
    rep["my_global_from_matlab_setup_vs_matlab_global"] = dict(
        lb_bit_equal=int(bit_equal(lbg, lb2).sum()), ub_bit_equal=int(bit_equal(ubg, ub2).sum()), n=len(rx),
        sha256=sha_bounds(lbg, ubg))
    rep["matlab_setup_to_global_changes"] = dict(lb=int((~bit_equal(lb1, lb2)).sum()), ub=int((~bit_equal(ub1, ub2)).sum()))
    # setup vs stored
    rep["setup_vs_stored"] = dict(n_lb_changed=int((~bit_equal(s_lb, lb1)).sum()), n_ub_changed=int((~bit_equal(s_ub, ub1)).sum()),
                                  n_lb_changed_numeric=int((s_lb != lb1).sum()), n_ub_changed_numeric=int((s_ub != ub1).sum()),
                                  groups=classify(rx, s_lb, s_ub, lb1, ub1))
    kid = rep["setup_vs_stored"]["groups"].get("kidney")
    if kid:
        rep["implied_old_gfr_ml_min"] = dict(from_min_ratio=90.0 / kid["ratio_max"], from_max_ratio=90.0 / kid["ratio_min"])
    csf = rep["setup_vs_stored"]["groups"].get("csf_export")
    if csf:
        rep["csf_ratio_vs_0.52/0.35"] = dict(min=csf["ratio_min"], max=csf["ratio_max"], expected=0.52 / 0.35)
    # setup changes on entries that the global step overwrites
    lbt, ubt, _ = global_step(rx, np.full(len(rx), 7.0), np.full(len(rx), 7.0))
    tl, tu = lbt != 7.0, ubt != 7.0
    rep["setup_changes_on_overwritten_entries"] = dict(
        n_lb_written=int(tl.sum()), n_ub_written=int(tu.sum()),
        n_setup_lb_changes_there=int(((~bit_equal(s_lb, lb1)) & tl).sum()),
        n_setup_ub_changes_there=int(((~bit_equal(s_ub, ub1)) & tu).sum()))
    # shipped bounds + global step (what a 'shipped' setup would give)
    lbs, ubs, _ = global_step(rx, s_lb, s_ub)
    rep["shipped_plus_global_sha256"] = sha_bounds(lbs, ubs)
    # carnitine diet uptake and a few named bounds
    named = {}
    for r in ("Diet_EX_crn[d]", "Diet_EX_ura[d]", "Diet_EX_glc_D[d]"):
        if r in rx:
            i = rx.index(r)
            named[r] = dict(stored=(float(s_lb[i]), float(s_ub[i])), setup=(float(lb1[i]), float(ub1[i])), glob=(float(lb2[i]), float(ub2[i])))
    rep["named_bounds"] = named
    rep["negative_zeros"] = dict(stored_lb=int(np.sum((s_lb == 0) & np.signbit(s_lb))), stored_ub=int(np.sum((s_ub == 0) & np.signbit(s_ub))),
                                 glob_lb=int(np.sum((lb2 == 0) & np.signbit(lb2))), glob_ub=int(np.sum((ub2 == 0) & np.signbit(ub2))))
    return rep


def main():
    out = OrderedDict()
    for model, cfg in MODELS.items():
        out[model] = one(model, cfg)
        r = out[model]
        print(f"== {model}: {r['n_rxns']} order {r['same_order']} mets {r['n_mets']} coupling {r['has_coupling']} "
              f"bile-duct exits {r['n_bileduct_exits']}")
        print(f"   recorded {r['recorded_lp_bounds_sha256']}")
        print(f"   MATLAB global sha {r['matlab_global_sha256']}")
        print(f"   my global step {r['my_global_step_counts']} vs MATLAB global {r['my_global_from_matlab_setup_vs_matlab_global']}")
        print(f"   setup->global changes {r['matlab_setup_to_global_changes']}")
        sv = r["setup_vs_stored"]
        print(f"   setup vs stored: lb {sv['n_lb_changed']} ub {sv['n_ub_changed']} (numeric {sv['n_lb_changed_numeric']}/{sv['n_ub_changed_numeric']})")
        for k, g in sv["groups"].items():
            print(f"      {k:16s} n {g['n']:5d} {g['by_bound']} ratios {g['ratios']} {g.get('ratio_min', '')} {g.get('ratio_max', '')}")
            for e in g["examples"][:2]:
                print(f"            {e}")
        print(f"   implied old GFR {r.get('implied_old_gfr_ml_min')}; CSF ratio {r.get('csf_ratio_vs_0.52/0.35')}")
        print(f"   setup changes on overwritten entries {r['setup_changes_on_overwritten_entries']}")
        print(f"   shipped+global sha {r['shipped_plus_global_sha256']}")
        print(f"   named {r['named_bounds']}; -0 {r['negative_zeros']}")
    dump("bounds_v04.json", out)


if __name__ == "__main__":
    sys.exit(main())
