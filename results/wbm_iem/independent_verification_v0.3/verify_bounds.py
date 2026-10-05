"""Independent bounds check (claim 5), without gembench.wbm_constraints (or any gembench module).

1. Load the model file (my own scipy.io reader) and the two MATLAB bound files; check reaction order.
2. sha256 of (lb as little-endian float64 bytes) + (ub as little-endian float64 bytes) of the MATLAB global
   bounds, raw and with -0.0 normalised to +0.0, against 'lp_bounds_sha256' of the v0.3 results.
3. My own implementation of the "unified reaction constraints" block of runIEM_HH.m (executed copy, lines 140-184)
   applied to (a) the MATLAB setup bounds -> compare with the MATLAB global bounds bit for bit, and hash;
   (b) the bounds stored in the model file -> hash, against 'lp_bounds_sha256' of the v0.2b results (shipped setup).
4. Classify every difference between the MATLAB setup bounds and the stored model bounds into the changes documented
   in docs/studies/wbm-iem-v0.3-plan.md section 1.2 (kidney x0.694, CSF export x1.486, 3 BBB uptakes, 31 diet uptakes
   at 0.1, 5 bile acids at 10 +/- 20%).
Writes bounds.json next to this script.
"""
import hashlib
import json
import os
from collections import Counter, OrderedDict

import numpy as np
import scipy.io as sio

ROOT = "/home/claude/mma"
HERE = os.path.dirname(os.path.abspath(__file__))
MODEL = os.path.join(ROOT, "external/COBRA.models/mat/Harvey_1_03d.mat")
SETUP = os.path.join(ROOT, "results/wbm_iem/matlab_reference/step1/matlab_setup_bounds_Harvey_1_03d.mat")
GLOBAL = os.path.join(ROOT, "results/wbm_iem/matlab_reference/step2/matlab_global_bounds_Harvey_1_03d.mat")


def cellstr(a):
    out = []
    for x in np.asarray(a, dtype=object).ravel():
        x = np.asarray(x)
        out.append(str(x.ravel()[0]) if x.size else "")
    return out


def load_model_bounds(path):
    d = sio.loadmat(path, squeeze_me=False, struct_as_record=True)
    keys = [k for k in d if not k.startswith("__")]
    assert len(keys) == 1, keys
    m = d[keys[0]][0, 0]
    return keys[0], cellstr(m["rxns"]), np.asarray(m["lb"], dtype=np.float64).ravel(), \
        np.asarray(m["ub"], dtype=np.float64).ravel()


def load_bounds(path):
    d = sio.loadmat(path, squeeze_me=False)
    return cellstr(d["rxns"]), np.asarray(d["lb"], dtype=np.float64).ravel(), np.asarray(d["ub"], dtype=np.float64).ravel()


def sha(lb, ub):
    return hashlib.sha256(np.ascontiguousarray(lb, dtype="<f8").tobytes() +
                          np.ascontiguousarray(ub, dtype="<f8").tobytes()).hexdigest()


def norm0(a):
    a = a.copy()
    a[a == 0] = 0.0          # -0.0 == 0 is True; assigning +0.0 drops the sign
    return a


def n_negzero(a):
    return int(np.sum((a == 0) & np.signbit(a)))


def bitwise_equal(a, b):
    return np.ascontiguousarray(a, dtype="<f8").view(np.uint64) == np.ascontiguousarray(b, dtype="<f8").view(np.uint64)


RNEW = ['BileDuct_EX_12dhchol[bd]_[luSI]', 'BileDuct_EX_3dhcdchol[bd]_[luSI]', 'BileDuct_EX_3dhchol[bd]_[luSI]',
        'BileDuct_EX_3dhdchol[bd]_[luSI]', 'BileDuct_EX_3dhlchol[bd]_[luSI]', 'BileDuct_EX_7dhcdchol[bd]_[luSI]',
        'BileDuct_EX_7dhchol[bd]_[luSI]', 'BileDuct_EX_cdca24g[bd]_[luSI]', 'BileDuct_EX_cdca3g[bd]_[luSI]',
        'BileDuct_EX_cholate[bd]_[luSI]', 'BileDuct_EX_dca24g[bd]_[luSI]', 'BileDuct_EX_dca3g[bd]_[luSI]',
        'BileDuct_EX_dchac[bd]_[luSI]', 'BileDuct_EX_dgchol[bd]_[luSI]', 'BileDuct_EX_gchola[bd]_[luSI]',
        'BileDuct_EX_hca24g[bd]_[luSI]', 'BileDuct_EX_hca6g[bd]_[luSI]', 'BileDuct_EX_hdca24g[bd]_[luSI]',
        'BileDuct_EX_hdca6g[bd]_[luSI]', 'BileDuct_EX_hyochol[bd]_[luSI]', 'BileDuct_EX_icdchol[bd]_[luSI]',
        'BileDuct_EX_isochol[bd]_[luSI]', 'BileDuct_EX_lca24g[bd]_[luSI]', 'BileDuct_EX_tchola[bd]_[luSI]',
        'BileDuct_EX_tdchola[bd]_[luSI]', 'BileDuct_EX_tdechola[bd]_[luSI]', 'BileDuct_EX_thyochol[bd]_[luSI]',
        'BileDuct_EX_uchol[bd]_[luSI]']   # copied from runIEM_HH_ref.m line 181 (Rnew)


def unified_constraints(rxns, lb, ub):
    """runIEM_HH.m lines 140-182, transcribed: strfind substring matches, setdiff, ismember."""
    lb, ub = lb.copy(), ub.copy()
    R = ['_ARGSL', '_GACMTRc', '_FUM', '_FUMm', '_HMR_7698', '_UAG4E', '_UDPG4E', '_GALT', '_G6PDH2c', '_G6PDH2r',
         '_G6PDH2rer', '_GLUTCOADHm', '_r0541', '_ACOAD8m', '_RE2410C', '_RE2410N']
    R2 = ['_FUMt', '_FUMAC', '_FUMS', 'BBB']
    all2 = {r for r in rxns if any(p in r for p in R)}
    all4 = {r for r in rxns if any(p in r for p in R2)}
    mic = {r for r in rxns if "Micro_" in r}
    irr = (all2 - all4) - mic
    mask_irr = np.array([r in irr for r in rxns])
    lb[mask_irr] = 0.0
    X = {r for r in rxns if "_r0784" in r or "_r0463" in r} - mic
    mask_x = np.array([r in X for r in rxns])
    lb[mask_x] = 0.0
    ub[mask_x] = 0.0
    rnew = set(RNEW)
    mask_b = np.array([r in rnew for r in rxns])
    ub[mask_b] = 100.0
    return lb, ub, {"set_irreversible": int(mask_irr.sum()), "closed": int(mask_x.sum()), "bile_duct_ub_100": int(mask_b.sum())}


def classify(rxns, stored_lb, stored_ub, lb, ub):
    gfr_ratio = 90.0 / 129.7489655616        # documented: GFR 90 instead of 129.75 ml/min
    csf_ratio = 0.52 / 0.35                  # documented: 0.52 instead of 0.35 ml/min
    cats = Counter()
    examples = {}
    other = []
    for name, s, n in (("lb", stored_lb, lb), ("ub", stored_ub, ub)):
        idx = np.where(s != n)[0]
        for i in idx:
            r, a, b = rxns[i], s[i], n[i]
            ratio = b / a if a != 0 else None
            if r.startswith("Kidney_EX_") and ratio is not None and abs(ratio - gfr_ratio) < 1e-9:
                c = f"kidney x{gfr_ratio:.4f} ({name})"
            elif r.startswith("BBB_") and r.endswith("[CSF]exp") and name == "ub" and ratio is not None and abs(ratio - csf_ratio) < 1e-9:
                c = f"CSF export x{csf_ratio:.4f} (ub)"
            elif r in ("BBB_TRP_L[CSF]upt", "BBB_LKYNR[CSF]upt", "BBB_KYNATE[CSF]upt") and name == "lb":
                c = "BBB uptake constrained (lb)"
            elif r.startswith("Diet_EX_") and name == "lb" and a == 0 and b == -0.1:
                c = "diet uptake opened at 0.1 (lb)"
            elif r.startswith("Diet_EX_") and a == 0 and ((name == "lb" and b == -12.0) or (name == "ub" and b == -8.0)):
                c = f"diet 10 +/-20% ({name})"
            else:
                c = "other"
                other.append({"rxn": r, "bound": name, "stored": float(a), "matlab_setup": float(b)})
            cats[c] += 1
            examples.setdefault(c, {"rxn": r, "bound": name, "stored": float(a), "matlab_setup": float(b)})
    return dict(cats), examples, other


def main():
    rep = OrderedDict()
    var, rx_m, lb_s, ub_s = load_model_bounds(MODEL)
    rx_1, lb_1, ub_1 = load_bounds(SETUP)
    rx_2, lb_2, ub_2 = load_bounds(GLOBAL)
    rep["model_variable"] = var
    rep["n_rxns"] = {"model": len(rx_m), "matlab_setup": len(rx_1), "matlab_global": len(rx_2)}
    rep["same_order_setup"] = rx_1 == rx_m
    rep["same_order_global"] = rx_2 == rx_m
    rep["n_unique_rxns"] = len(set(rx_m))
    rep["negative_zeros"] = {"stored_lb": n_negzero(lb_s), "stored_ub": n_negzero(ub_s), "setup_lb": n_negzero(lb_1),
                             "setup_ub": n_negzero(ub_1), "global_lb": n_negzero(lb_2), "global_ub": n_negzero(ub_2)}
    rep["nan_count"] = {k: int(np.isnan(v).sum()) for k, v in
                        (("stored_lb", lb_s), ("stored_ub", ub_s), ("global_lb", lb_2), ("global_ub", ub_2))}
    results = {}
    for label in ("v0.2b", "v0.3"):
        recs = json.load(open(os.path.join(ROOT, f"results/wbm_iem/Harvey_1_03d_iem_results_{label}.json")))
        hs = {r["provenance"]["lp_bounds_sha256"] for r in recs}
        assert len(hs) == 1
        results[label] = hs.pop()
    rep["recorded_lp_bounds_sha256"] = results
    # 2. MATLAB global bounds against v0.3
    h_raw, h_norm = sha(lb_2, ub_2), sha(norm0(lb_2), norm0(ub_2))
    rep["matlab_global_sha256"] = {"raw": h_raw, "neg_zero_normalised": h_norm,
                                   "raw_matches_v0.3": h_raw == results["v0.3"],
                                   "normalised_matches_v0.3": h_norm == results["v0.3"]}
    # 3a. my unified constraints on the MATLAB setup bounds
    lb_g, ub_g, counts = unified_constraints(rx_m, lb_1, ub_1)
    rep["unified_constraints_counts_on_setup"] = counts
    rep["my_global_from_matlab_setup_vs_matlab_global"] = {
        "lb_bitwise_equal": int(bitwise_equal(lb_g, lb_2).sum()), "ub_bitwise_equal": int(bitwise_equal(ub_g, ub_2).sum()),
        "lb_numeric_equal": int((lb_g == lb_2).sum()), "ub_numeric_equal": int((ub_g == ub_2).sum()),
        "sha256": sha(lb_g, ub_g), "sha256_matches_v0.3": sha(lb_g, ub_g) == results["v0.3"]}
    # how many entries the unified constraints change, and MATLAB setup -> MATLAB global differences
    rep["matlab_setup_to_global"] = {"lb_changed": int((lb_1 != lb_2).sum()), "ub_changed": int((ub_1 != ub_2).sum())}
    # Do any setup changes (MATLAB setup vs stored) sit on entries that the unified constraints overwrite?
    # If not, the global-bounds hash covers every setup change, and the setup bounds at overwritten entries equal
    # the stored ones in MATLAB.
    lb_t, ub_t, _ = unified_constraints(rx_m, np.full(len(rx_m), 7.0), np.full(len(rx_m), 7.0))
    touched_lb, touched_ub = lb_t != 7.0, ub_t != 7.0
    rep["setup_changes_on_overwritten_entries"] = {
        "n_lb_entries_written_by_unified_constraints": int(touched_lb.sum()),
        "n_ub_entries_written_by_unified_constraints": int(touched_ub.sum()),
        "n_setup_lb_changes_there": int(((lb_s != lb_1) & touched_lb).sum()),
        "n_setup_ub_changes_there": int(((ub_s != ub_1) & touched_ub).sum())}
    # 3b. my unified constraints on the stored model bounds (v0.2b = shipped setup + toolbox bile-duct list)
    lb_b, ub_b, counts_b = unified_constraints(rx_m, lb_s, ub_s)
    hb, hbn = sha(lb_b, ub_b), sha(norm0(lb_b), norm0(ub_b))
    rep["shipped_plus_unified"] = {"counts": counts_b, "sha256": hb, "sha256_neg_zero_normalised": hbn,
                                   "matches_v0.2b": hb == results["v0.2b"], "normalised_matches_v0.2b": hbn == results["v0.2b"]}
    # 4. classification of MATLAB setup vs stored
    cats, ex, other = classify(rx_m, lb_s, ub_s, lb_1, ub_1)
    rep["matlab_setup_vs_stored"] = {"n_lb_changed": int((lb_s != lb_1).sum()), "n_ub_changed": int((ub_s != ub_1).sum()),
                                     "categories": cats, "examples": ex, "unclassified": other[:50],
                                     "n_unclassified": len(other)}
    # same classification between the MATLAB global bounds and the stored bounds after my unified constraints
    cats2, _, other2 = classify(rx_m, lb_b, ub_b, lb_2, ub_2)
    rep["matlab_global_vs_shipped_plus_unified"] = {"n_lb_changed": int((lb_b != lb_2).sum()),
                                                    "n_ub_changed": int((ub_b != ub_2).sum()), "categories": cats2,
                                                    "n_unclassified": len(other2), "unclassified": other2[:20]}
    with open(os.path.join(HERE, "bounds.json"), "w") as fh:
        json.dump(rep, fh, indent=1)
        fh.write("\n")
    print(json.dumps(rep, indent=1)[:8000])


if __name__ == "__main__":
    main()
