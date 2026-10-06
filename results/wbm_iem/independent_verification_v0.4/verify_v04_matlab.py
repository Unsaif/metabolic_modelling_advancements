"""Independent comparison of the MATLAB runIEM_HH reference runs with the Python maxima, for Harvetta (v0.4, run M
vs run A) and Harvey (v0.3 step 2 vs the v0.3 results).

MATLAB calls are recomputed from MATLAB's printed values exactly as runIEM_HH does (H_D = str2num(H) - str2num(D);
Increased if H_D < -1e-6, Decreased if H_D > 1e-6, else unchanged; NaN or [] -> unchanged). Python calls follow the
protocol (values with |f| <= 1e-6 set to 0, then the 1e-6 threshold). Imports nothing from gembench or scripts/.
Writes matlab_v04.json here.
"""
import os
import sys
from collections import Counter, OrderedDict

from vcommon import (DIRS, TOL, biomarker, call_protocol, dump, index_results, load_json, matlab_call, max_ok, mfin,
                     parse_matlab_iem, protocol_rows, ROOT, PROTOCOL)

RUNS = OrderedDict([
    ("Harvetta", dict(mat="results/wbm_iem/matlab_reference/harvetta_step2/matlab_iem_results_Harvetta_1_03d.mat",
                      py="results/wbm_iem/Harvetta_1_03d_iem_results_v0.4.json")),
    ("Harvey", dict(mat="results/wbm_iem/matlab_reference/step2/matlab_iem_results_Harvey_1_03d.mat",
                    py="results/wbm_iem/Harvey_1_03d_iem_results_v0.3.json")),
])


def compare(model, cfg, rows):
    mat, scal = parse_matlab_iem(os.path.join(ROOT, cfg["mat"]))
    by = index_results(load_json(cfg["py"]))
    rep = OrderedDict(n_iem_blocks=len(mat), stored=scal)
    prot_iems = [r["iem"] for r in rows]
    rep["iems_missing_in_matlab"] = sorted(set(prot_iems) - set(mat))
    rep["iems_extra_in_matlab"] = sorted(set(mat) - set(prot_iems))
    flat, mismatch = [], []
    counts = Counter(r["iem"] for r in rows)
    for iem, n in counts.items():
        if len(mat[iem]["rows"]) != n:
            mismatch.append((iem, "count", len(mat[iem]["rows"]), n))
    for r in rows:
        m = mat[r["iem"]]["rows"][r["pos"]]
        if m["reaction"] != r["reaction"] or m["label"] != r["label"]:
            mismatch.append((r["iem"], r["pos"], m["reaction"], r["reaction"], m["label"], r["label"]))
        b = biomarker(by, r)
        ok = max_ok(b)
        mexp = "Increased" if "Incre" in m["label"] else ("Decreased" if "Decre" in m["label"] else "Unchanged")
        flat.append(OrderedDict(
            iem=r["iem"], pos=r["pos"], reaction=r["reaction"], expected=r["expected"], matlab_expected=mexp,
            m_h=m["healthy"], m_d=m["disease"], m_hs=m["healthy_str"], m_ds=m["disease_str"], label_prefix=m["label_prefix"],
            m_call=matlab_call(m["healthy"], m["disease"]),
            p_h=b.get("healthy"), p_d=b.get("disease"), p_ok=ok, p_status=(b.get("status_healthy"), b.get("status_disease")),
            p_call=call_protocol(b["healthy"], b["disease"]) if ok else "NA"))
    rep["protocol_mismatches"] = mismatch
    rep["n_biomarkers"] = len(flat)
    rep["expected_label_mismatch"] = [(f["iem"], f["reaction"]) for f in flat if f["expected"] != f["matlab_expected"]]
    # runIEM_HH's own tally
    tally = Counter((f["expected"], f["m_call"]) for f in flat)
    upup, dodo = tally[("Increased", "Increased")], tally[("Decreased", "Decreased")]
    updo = tally[("Increased", "Decreased")]
    rep["runiem"] = OrderedDict(correct=upup + dodo, n=len(flat), accuracy=(upup + dodo) / len(flat),
                                precision=upup / (upup + updo) if upup + updo else None,
                                tally={f"{a}->{b}": n for (a, b), n in sorted(tally.items())},
                                n_unique_biomarkers=len({f["reaction"] for f in flat}),
                                matches_stored_accuracy=abs((upup + dodo) / len(flat) - scal.get("Accuracy", -1)) < 1e-12)
    # non-finite MATLAB values
    nonfin = [f for f in flat if not (mfin(f["m_h"]) and mfin(f["m_d"]))]
    rep["n_nonfinite"] = len(nonfin)
    rep["nonfinite_side"] = dict(Counter(("H" if not mfin(f["m_h"]) else "") + ("D" if not mfin(f["m_d"]) else "") for f in nonfin))
    rep["nonfinite_strings"] = dict(Counter((f["m_hs"], f["m_ds"]) if not (mfin(f["m_h"]) and mfin(f["m_d"])) else None for f in nonfin))
    absent = [f for f in nonfin if f["m_h"] is None and f["m_d"] is None]
    nan_any = [f for f in nonfin if f not in absent]
    rep["absent_or_empty"] = [(f["iem"], f["reaction"], f["m_hs"], f["m_ds"], f["label_prefix"], f["p_status"]) for f in absent]
    rep["n_without_optimum"] = len(nan_any)
    rep["without_optimum_by_iem"] = dict(sorted(Counter(f["iem"] for f in nan_any).items(), key=lambda x: (-x[1], x[0])))
    rep["without_optimum"] = [dict(iem=f["iem"], reaction=f["reaction"], expected=f["expected"], m_h=f["m_hs"], m_d=f["m_ds"],
                                   p_h=f["p_h"], p_d=f["p_d"], p_status=f["p_status"], p_call=f["p_call"], m_call=f["m_call"],
                                   python_correct=f["p_call"] == f["expected"]) for f in nan_any]
    rep["without_optimum_python_all_optimal"] = all(f["p_ok"] for f in nan_any)
    rep["without_optimum_python_correct"] = sum(f["p_call"] == f["expected"] for f in nan_any)
    rep["without_optimum_python_wrong"] = [(f["iem"], f["reaction"], f["expected"], f["p_call"], f["p_h"], f["p_d"])
                                           for f in nan_any if f["p_call"] != f["expected"]]
    rep["without_optimum_same_call"] = [(f["iem"], f["reaction"], f["m_call"], f["p_call"]) for f in nan_any if f["m_call"] == f["p_call"]]
    # calls where both MATLAB values are finite
    fin = [f for f in flat if mfin(f["m_h"]) and mfin(f["m_d"])]
    fin_dir = [f for f in fin if f["expected"] in DIRS]
    rep["n_matlab_finite"] = len(fin)
    rep["n_matlab_finite_directional"] = len(fin_dir)
    rep["finite_python_not_ok"] = [(f["iem"], f["reaction"]) for f in fin if not f["p_ok"]]
    rep["finite_same_call"] = sum(f["m_call"] == f["p_call"] for f in fin)
    rep["finite_different_call"] = [(f["iem"], f["reaction"], f["m_hs"], f["m_ds"], f["m_call"], f["p_h"], f["p_d"], f["p_call"])
                                    for f in fin if f["m_call"] != f["p_call"]]
    rep["finite_matlab_correct"] = sum(f["m_call"] == f["expected"] for f in fin_dir)
    rep["finite_python_correct"] = sum(f["p_call"] == f["expected"] for f in fin_dir)
    rep["accuracy_on_finite"] = rep["finite_matlab_correct"] / len(fin_dir) if fin_dir else None
    # all calls, both present
    present = [f for f in flat if f["p_call"] != "NA"]
    diffs = [f for f in present if f["m_call"] != f["p_call"]]
    rep["n_present_python"] = len(present)
    rep["n_same_call_all"] = len(present) - len(diffs)
    rep["n_diff_call_all"] = len(diffs)
    rep["diffs_all"] = [(f["iem"], f["reaction"], f["expected"], f["m_hs"], f["m_ds"], f["m_call"], f["p_call"]) for f in diffs]
    rep["diffs_with_matlab_healthy_nonfinite"] = sum(not mfin(f["m_h"]) for f in diffs)
    rep["diffs_with_matlab_disease_nonfinite"] = sum(not mfin(f["m_d"]) for f in diffs)
    # values
    vals = []
    for f in present:
        for side, m, p, s in (("healthy", f["m_h"], f["p_h"], f["m_hs"]), ("disease", f["m_d"], f["p_d"], f["m_ds"])):
            if mfin(m) and p is not None:
                scale = max(abs(m), abs(p))
                vals.append(dict(iem=f["iem"], reaction=f["reaction"], side=side, matlab=s, python=p, abs_diff=abs(m - p),
                                 rel=abs(m - p) / scale if scale > 0 else 0.0, scale=scale, minmag=min(abs(m), abs(p))))
    big = [v for v in vals if v["scale"] > 1e-3]
    bigmin = [v for v in vals if v["minmag"] > 1e-3]
    small = [v for v in vals if v["scale"] <= 1e-3]
    rep["values"] = OrderedDict(
        n_pairs=len(vals), n_pairs_scale_over_1e3=len(big), max_rel_scale_over_1e3=max(v["rel"] for v in big),
        n_pairs_minmag_over_1e3=len(bigmin), max_rel_minmag_over_1e3=max(v["rel"] for v in bigmin),
        n_over_0p008pct=sum(v["rel"] > 8e-5 for v in big), n_over_0p03pct=sum(v["rel"] > 3e-4 for v in big),
        n_over_0p1pct=sum(v["rel"] > 1e-3 for v in big),
        largest=sorted(big, key=lambda v: -v["rel"])[:6],
        n_pairs_scale_le_1e3=len(small), max_abs_diff_scale_le_1e3=max((v["abs_diff"] for v in small), default=None),
        largest_small=sorted(small, key=lambda v: -v["abs_diff"])[:3])
    # header rows: IEM-flux maxima
    py_rec = {r["iem"]: r for r in load_json(cfg["py"])}
    vm = []
    for iem, blk in mat.items():
        try:
            hv = float(blk["header"][0][1])
        except ValueError:
            hv = None
        pv = py_rec[iem]["vmax_healthy"]
        if hv is not None and pv is not None:
            vm.append((iem, hv, pv, abs(hv - pv) / max(abs(hv), abs(pv), 1e-12)))
        else:
            vm.append((iem, blk["header"][0][1], pv, None))
    rep["vmax_healthy"] = dict(n=len(vm), n_rel_over_1e4=sum(1 for x in vm if x[3] is not None and x[3] > 1e-4),
                               worst=sorted([x for x in vm if x[3] is not None], key=lambda x: -x[3])[:4],
                               nonfinite=[x for x in vm if x[3] is None])
    # named rows of interest
    rep["named"] = {f"{f['iem']} {f['reaction']}": dict(m_h=f["m_hs"], m_d=f["m_ds"], m_call=f["m_call"], p_h=f["p_h"], p_d=f["p_d"],
                                                       p_call=f["p_call"], expected=f["expected"])
                    for f in flat if (f["iem"], f["reaction"]) in {("ASNSD", "DM_asn_L[bc]"), ("ASNSD", "DM_asn_L[csf]"),
                                                                    ("BTD", "EX_3hpp[u]"), ("LNS", "DM_fol[bc]"), ("PC", "DM_glc_D[bc]"),
                                                                    ("MMA", "DM_crn[bc]"), ("FED", "DM_chsterol[bc]"),
                                                                    ("HPC", "EX_25aics[u]")}}
    return rep, flat


def main():
    rows = protocol_rows(load_json(os.path.relpath(PROTOCOL, ROOT)))
    out = OrderedDict()
    flats = {}
    for model, cfg in RUNS.items():
        out[model], flats[model] = compare(model, cfg, rows)
    # cross-model: biomarkers without a MATLAB optimum in both
    a = {(f["iem"], f["reaction"]) for f in flats["Harvetta"] if not (mfin(f["m_h"]) and mfin(f["m_d"]))}
    b = {(f["iem"], f["reaction"]) for f in flats["Harvey"] if not (mfin(f["m_h"]) and mfin(f["m_d"]))}
    out["nonfinite_in_both_models"] = sorted(a & b)
    out["nonfinite_harvey_only"] = sorted(b - a)
    dump("matlab_v04.json", out)
    for model, r in out.items():
        if not isinstance(r, dict) or "runiem" not in r:
            continue
        print(f"== {model}: blocks {r['n_iem_blocks']} biomarkers {r['n_biomarkers']} mismatches {r['protocol_mismatches']} "
              f"label mismatches {r['expected_label_mismatch']}")
        print(f"   runIEM_HH {r['runiem']['correct']}/{r['runiem']['n']} = {r['runiem']['accuracy']:.4f} (stored {r['stored']}) "
              f"match {r['runiem']['matches_stored_accuracy']} unique {r['runiem']['n_unique_biomarkers']}")
        print(f"   nonfinite {r['n_nonfinite']} sides {r['nonfinite_side']} strings {r['nonfinite_strings']}")
        print(f"   absent {r['absent_or_empty']}")
        print(f"   without optimum {r['n_without_optimum']} by IEM {r['without_optimum_by_iem']}")
        print(f"   python all optimal there {r['without_optimum_python_all_optimal']}, correct {r['without_optimum_python_correct']}, "
              f"wrong {r['without_optimum_python_wrong']}, same call {r['without_optimum_same_call']}")
        print(f"   finite {r['n_matlab_finite']} directional {r['n_matlab_finite_directional']} same call {r['finite_same_call']} "
              f"diff {r['finite_different_call']} matlab correct {r['finite_matlab_correct']} python correct {r['finite_python_correct']} "
              f"acc {r['accuracy_on_finite']:.4f}; python not ok {r['finite_python_not_ok']}")
        print(f"   all present {r['n_present_python']} same {r['n_same_call_all']} diff {r['n_diff_call_all']} "
              f"(healthy nonfinite {r['diffs_with_matlab_healthy_nonfinite']}, disease nonfinite {r['diffs_with_matlab_disease_nonfinite']})")
        v = r["values"]
        print(f"   values: pairs {v['n_pairs']} scale>1e-3 {v['n_pairs_scale_over_1e3']} max rel {v['max_rel_scale_over_1e3']:.3e} "
              f"(>0.008% {v['n_over_0p008pct']}, >0.03% {v['n_over_0p03pct']}, >0.1% {v['n_over_0p1pct']}); minmag>1e-3 "
              f"{v['n_pairs_minmag_over_1e3']} max rel {v['max_rel_minmag_over_1e3']:.3e}; small {v['n_pairs_scale_le_1e3']} max abs "
              f"{v['max_abs_diff_scale_le_1e3']}")
        for x in v["largest"]:
            print(f"      {x['iem']:6s} {x['reaction']:18s} {x['side']:7s} MATLAB {x['matlab']:>12s} Python {x['python']:.10g} rel {x['rel']:.3e}")
        print(f"   vmax healthy: {r['vmax_healthy']}")
        for k, x in r["named"].items():
            print(f"   named {k}: {x}")
        for x in r["without_optimum"]:
            print(f"      NaN: {x['iem']:6s} {x['reaction']:18s} exp {x['expected']:9s} MATLAB h={x['m_h']} d={x['m_d']} | Python "
                  f"h={x['p_h']:.6g} d={x['p_d']:.6g} {x['p_call']} {'ok' if x['python_correct'] else 'WRONG'}")
    print("nonfinite in both models:", out["nonfinite_in_both_models"])
    print("nonfinite Harvey only:", out["nonfinite_harvey_only"])


if __name__ == "__main__":
    sys.exit(main())
