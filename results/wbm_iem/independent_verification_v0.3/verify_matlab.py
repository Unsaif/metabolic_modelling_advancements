"""Independent comparison of the MATLAB reference run (runIEM_HH.m, Gurobi) with the Python runs (claim 4).

Written by the independent verifier; imports nothing from gembench or the project's comparison scripts.
Reads results/wbm_iem/matlab_reference/step2/matlab_iem_results_Harvey_1_03d.mat (IEMSol_<IEM> cell arrays),
the protocol and the Python results files. Writes matlab_comparison.json next to this script.

MATLAB parsing (from runIEM_HH.m, lines 2703-2755 of the executed copy):
  rows 5, 7, 9, ... (1-based) are 'Healthy:<rxn>', the next row 'Disease:<rxn>'; value string in column 2,
  'Disease - Reported:<label>' in column 3.
  H_D = str2num(healthy) - str2num(disease); Increased if H_D < -1e-6, Decreased if H_D > 1e-6, else unchanged.
  str2num of 'NaN' is NaN and of 'NA' or '' is [] (empty); any comparison with NaN or [] is false -> unchanged.
  in vivo: 'Incre' in column 3 -> up, 'Decre' -> down, else unchanged.
  Accuracy = (UpUp + DoDo) / (all rows).
"""
import json
import math
import os
from collections import Counter, OrderedDict

import numpy as np
import scipy.io as sio

ROOT = "/home/claude/mma"
HERE = os.path.dirname(os.path.abspath(__file__))
MAT = os.path.join(ROOT, "results/wbm_iem/matlab_reference/step2/matlab_iem_results_Harvey_1_03d.mat")
RUNS = OrderedDict([
    ("v0.2", "results/wbm_iem/Harvey_1_03d_iem_results_v0.2.json"),
    ("v0.2b", "results/wbm_iem/Harvey_1_03d_iem_results_v0.2b.json"),
    ("v0.3", "results/wbm_iem/Harvey_1_03d_iem_results_v0.3.json"),
])


def cell_text(c):
    a = np.asarray(c)
    if a.size == 0:
        return ""
    if a.dtype.kind in "US":
        return str(a.ravel()[0])
    raise TypeError(f"unexpected cell content {a!r}")


def str2num(s):
    """Return a float, or None for MATLAB's empty result ([])."""
    s = s.strip()
    if s == "":
        return None
    if s in ("NaN", "nan"):
        return float("nan")
    if s in ("Inf", "inf"):
        return float("inf")
    if s in ("-Inf", "-inf"):
        return float("-inf")
    try:
        return float(s)
    except ValueError:
        return None          # 'NA', 'ND', ... -> eval fails -> []


def matlab_call(h, d):
    if h is None or d is None:
        return "Unchanged"
    hd = h - d               # NaN propagates; NaN comparisons are false
    if hd < -1e-6:
        return "Increased"
    if hd > 1e-6:
        return "Decreased"
    return "Unchanged"


def invivo(text):
    return "Increased" if "Incre" in text else "Decreased" if "Decre" in text else "Unchanged"


def py_call(h, d):
    diff = d - h
    return "Increased" if diff > 1e-6 else "Decreased" if diff < -1e-6 else "Unchanged"


def fin(v):
    return v is not None and isinstance(v, (int, float)) and math.isfinite(v)


def parse(path):
    d = sio.loadmat(path, squeeze_me=False)
    out = OrderedDict()
    for key in sorted(k for k in d if k.startswith("IEMSol_")):
        cells = np.asarray(d[key], dtype=object)
        header = [(cell_text(cells[i, 0]), cell_text(cells[i, 1])) for i in range(4)]
        rows = []
        for j in range(4, cells.shape[0], 2):      # MATLAB j = 5:2:size(IEM,1) (0-based 4, 6, ...)
            hname, dname = cell_text(cells[j, 0]), cell_text(cells[j + 1, 0])
            assert hname.startswith("Healthy:") and dname.startswith("Disease:"), (key, j, hname, dname)
            rxn = hname[len("Healthy:"):]
            assert dname[len("Disease:"):] == rxn, (key, j)
            lab = cell_text(cells[j, 2])
            # runIEM_HH reads the label from row j (the Healthy row). checkIEM_WBM writes 'Disease - Reported:'
            # there, except for a biomarker reaction absent from the model ('Healthy - Reported:', HPC EX_25aics[u]).
            prefix = next((p for p in ("Disease - Reported:", "Healthy - Reported:") if lab.startswith(p)), None)
            assert prefix is not None, (key, j, lab)
            hs, ds = cell_text(cells[j, 1]), cell_text(cells[j + 1, 1])
            rows.append({"reaction": rxn, "label": lab[len(prefix):], "label_prefix": prefix, "healthy_str": hs,
                         "disease_str": ds, "healthy": str2num(hs), "disease": str2num(ds)})
        out[key[len("IEMSol_"):]] = {"header": header, "rows": rows}
    return out


def main():
    protocol = json.load(open(os.path.join(ROOT, "data/iem/iem_protocol_v0.2.json")))
    mat = parse(MAT)
    rep = OrderedDict()
    rep["n_iems"] = len(mat)
    prot_iems = [p["iem"] for p in protocol]
    rep["iems_missing_in_matlab"] = sorted(set(prot_iems) - set(mat))
    rep["iems_extra_in_matlab"] = sorted(set(mat) - set(prot_iems))
    # 1. match to the protocol by IEM and position
    mismatch = []
    for p in protocol:
        rows = mat.get(p["iem"], {"rows": []})["rows"]
        if len(rows) != len(p["biomarkers"]):
            mismatch.append((p["iem"], "count", len(rows), len(p["biomarkers"])))
        for k, (rxn, text) in enumerate(p["biomarkers"]):
            if k < len(rows) and (rows[k]["reaction"] != rxn or rows[k]["label"] != text):
                mismatch.append((p["iem"], k, rows[k]["reaction"], rxn, rows[k]["label"], text))
    rep["protocol_mismatches"] = mismatch
    flat = []
    for p in protocol:
        for k, (rxn, text) in enumerate(p["biomarkers"]):
            m = mat[p["iem"]]["rows"][k]
            flat.append({"iem": p["iem"], "call_index": p["call_index"], "pos": k, "reaction": rxn,
                         "expected": invivo(m["label"]), **m, "matlab_call": matlab_call(m["healthy"], m["disease"])})
    rep["n_matlab_biomarkers"] = len(flat)
    nonfin = [f for f in flat if not (fin(f["healthy"]) and fin(f["disease"]))]
    rep["n_matlab_nonfinite"] = len(nonfin)
    rep["matlab_nonfinite"] = [{"iem": f["iem"], "reaction": f["reaction"], "healthy": f["healthy_str"],
                                "disease": f["disease_str"]} for f in nonfin]
    rep["n_nonfinite_by_side"] = dict(Counter(("H" if not fin(f["healthy"]) else "") + ("D" if not fin(f["disease"]) else "")
                                              for f in nonfin))
    # 2. runIEM_HH's own tallies
    cnt = Counter((f["expected"], f["matlab_call"]) for f in flat)
    UpUp, DoDo = cnt[("Increased", "Increased")], cnt[("Decreased", "Decreased")]
    UpDo, DoUp = cnt[("Increased", "Decreased")], cnt[("Decreased", "Increased")]
    rep["runiem_tally"] = {f"{a}->{b}": n for (a, b), n in sorted(cnt.items())}
    rep["runiem_accuracy"] = {"correct": UpUp + DoDo, "n": len(flat), "accuracy": (UpUp + DoDo) / len(flat),
                              "precision": UpUp / (UpUp + UpDo), "fdr": UpDo / (UpUp + UpDo),
                              "n_unique_biomarkers": len({f["reaction"] for f in flat}), "n_diseases": len({f["iem"] for f in flat})}
    stored = sio.loadmat(MAT, squeeze_me=True)
    rep["runiem_accuracy"]["stored_Accuracy"] = float(stored["Accuracy"])
    rep["runiem_accuracy"]["stored_Precision"] = float(stored["Precision"])
    rep["runiem_accuracy"]["stored_NumBiomarkers"] = int(stored["NumBiomarkers"])
    # 3. compare with the Python runs
    rep["python"] = OrderedDict()
    for label, path in RUNS.items():
        res = {(r["iem"], r["call_index"]): r for r in json.load(open(os.path.join(ROOT, path)))}
        rows = []
        for f in flat:
            b = res[(f["iem"], f["call_index"])]["biomarkers"][f["pos"]]
            assert b["reaction"] == f["reaction"]
            ok = b["status_healthy"] == b["status_disease"] == "Optimal" and fin(b["healthy"]) and fin(b["disease"])
            rows.append({**f, "py_healthy": b["healthy"], "py_disease": b["disease"],
                         "py_call": py_call(b["healthy"], b["disease"]) if ok else "NA"})
        present = [r for r in rows if r["py_call"] != "NA"]
        diffs = [r for r in present if r["py_call"] != r["matlab_call"]]
        both = [r for r in present if fin(r["healthy"]) and fin(r["disease"]) and r["expected"] != "Unchanged"]
        # value agreement (healthy and disease separately, both finite)
        vals = []
        for r in present:
            for side in ("healthy", "disease"):
                m, p = r[side], r["py_" + side]
                if fin(m) and fin(p):
                    scale = max(abs(m), abs(p))
                    vals.append({"iem": r["iem"], "reaction": r["reaction"], "side": side, "matlab": m,
                                 "matlab_str": r[side + "_str"], "python": p, "abs_diff": abs(m - p),
                                 "rel_diff": abs(m - p) / scale if scale > 0 else 0.0, "scale": scale,
                                 "min_mag": min(abs(m), abs(p))})
        big = [v for v in vals if v["scale"] > 1e-3]
        big_min = [v for v in vals if v["min_mag"] > 1e-3]
        small = [v for v in vals if v["scale"] <= 1e-3]
        out = OrderedDict(
            n_present_in_both=len(present), n_same_call=len(present) - len(diffs), n_different_call=len(diffs),
            differences=[{"iem": r["iem"], "reaction": r["reaction"], "expected": r["expected"], "matlab_healthy": r["healthy_str"],
                          "matlab_disease": r["disease_str"], "matlab_call": r["matlab_call"], "py_healthy": r["py_healthy"],
                          "py_disease": r["py_disease"], "py_call": r["py_call"]} for r in diffs],
            n_differences_with_matlab_healthy_nonfinite=sum(not fin(r["healthy"]) for r in diffs),
            n_differences_with_matlab_disease_nonfinite=sum(not fin(r["disease"]) for r in diffs),
            n_differences_with_both_matlab_finite=sum(fin(r["healthy"]) and fin(r["disease"]) for r in diffs),
            not_present=[(r["iem"], r["reaction"], r["healthy_str"], r["disease_str"]) for r in rows if r["py_call"] == "NA"],
            n_both_finite_with_direction=len(both),
            n_matlab_correct_there=sum(r["matlab_call"] == r["expected"] for r in both),
            n_python_correct_there=sum(r["py_call"] == r["expected"] for r in both),
            n_values_compared=len(vals),
            n_values_scale_over_1e3=len(big),
            n_values_over_0p1pct_rel_among_scale_over_1e3=sum(v["rel_diff"] > 1e-3 for v in big),
            values_over_0p1pct=[v for v in big if v["rel_diff"] > 1e-3],
            n_values_min_magnitude_over_1e3=len(big_min),
            n_values_over_0p1pct_rel_among_min_mag_over_1e3=sum(v["rel_diff"] > 1e-3 for v in big_min),
            largest_rel_diffs_scale_over_1e3=sorted(big, key=lambda v: -v["rel_diff"])[:6],
            n_values_scale_le_1e3=len(small),
            max_abs_diff_scale_le_1e3=max((v["abs_diff"] for v in small), default=None),
        )
        rep["python"][label] = out
    # header rows: IEM-flux maxima (vmax) agreement, v0.3
    res3 = {r["iem"]: r for r in json.load(open(os.path.join(ROOT, RUNS["v0.3"])))}
    vm = []
    for iem, block in mat.items():
        hv = str2num(block["header"][0][1])
        pv = res3[iem]["vmax_healthy"]
        if fin(hv) and fin(pv):
            vm.append((iem, hv, pv, abs(hv - pv) / max(abs(hv), abs(pv), 1e-12)))
        else:
            vm.append((iem, block["header"][0][1], pv, None))
    rep["vmax_healthy_v0.3"] = {"n": len(vm), "n_rel_diff_over_1e-4": sum(1 for x in vm if x[3] is not None and x[3] > 1e-4),
                                "nonfinite": [x for x in vm if x[3] is None],
                                "worst": sorted([x for x in vm if x[3] is not None], key=lambda x: -x[3])[:5]}
    with open(os.path.join(HERE, "matlab_comparison.json"), "w") as fh:
        json.dump(rep, fh, indent=1, default=str)
        fh.write("\n")
    print(json.dumps({k: v for k, v in rep.items() if k != "python"}, indent=1, default=str)[:6000])
    for label, o in rep["python"].items():
        print(f"== {label}: present {o['n_present_in_both']} same {o['n_same_call']} diff {o['n_different_call']} "
              f"(MATLAB healthy nonfinite {o['n_differences_with_matlab_healthy_nonfinite']}, disease nonfinite "
              f"{o['n_differences_with_matlab_disease_nonfinite']}, both finite {o['n_differences_with_both_matlab_finite']}); "
              f"both finite+direction {o['n_both_finite_with_direction']}: MATLAB correct {o['n_matlab_correct_there']}, "
              f"Python correct {o['n_python_correct_there']}; not present {o['not_present']}")
        print(f"   values compared {o['n_values_compared']}; scale>1e-3: {o['n_values_scale_over_1e3']}, "
              f">0.1% rel: {o['n_values_over_0p1pct_rel_among_scale_over_1e3']}; min-mag>1e-3: {o['n_values_min_magnitude_over_1e3']}, "
              f">0.1%: {o['n_values_over_0p1pct_rel_among_min_mag_over_1e3']}; scale<=1e-3: {o['n_values_scale_le_1e3']} "
              f"max abs diff {o['max_abs_diff_scale_le_1e3']}")
        for v in o["largest_rel_diffs_scale_over_1e3"]:
            print(f"     {v['iem']:7s} {v['reaction']:18s} {v['side']:7s} MATLAB {v['matlab_str']:>12s} Python {v['python']:.8g} rel {v['rel_diff']:.3e}")
        if label == "v0.3":
            for dd in o["differences"]:
                print(f"     diff: {dd['iem']:6s} {dd['reaction']:18s} exp {dd['expected']:9s} MATLAB h={dd['matlab_healthy']} d={dd['matlab_disease']} "
                      f"call {dd['matlab_call']:9s} | Python h={dd['py_healthy']:.6g} d={dd['py_disease']:.6g} call {dd['py_call']}")


if __name__ == "__main__":
    main()
