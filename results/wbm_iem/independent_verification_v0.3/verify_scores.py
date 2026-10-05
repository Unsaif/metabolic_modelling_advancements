"""Independent recomputation of the IEM scores (claims 1-3) from the raw results files and the protocol.

Written by the independent verifier; imports nothing from gembench or the project's comparison scripts.
Reads: data/iem/iem_protocol_v0.2.json and results/wbm_iem/Harvey_1_03d_iem_results_{v0.2,v0.2b,v0.3}.json
Writes: scores.json (next to this script)

Definitions used here (my own, from the task statement and runIEM_HH.m):
  expected label  : 'Increased' if 'Incre' in the protocol label text, 'Decreased' if 'Decre', else 'Unchanged'
                    (the case-sensitive substring rule at the end of runIEM_HH.m).
  scored          : expected in {Increased, Decreased} and status_healthy == status_disease == 'Optimal'
                    and both healthy and disease values are finite numbers.
  call (tau = 0)  : d - h > 1e-6 -> Increased; < -1e-6 -> Decreased; else Unchanged (recomputed from the stored values).
  call (tau > 0)  : r = (d - h) / max(|h|, |d|) (r = 0 if both are 0); r > tau -> Increased; r < -tau -> Decreased.
"""
import json
import math
import os
import sys
from collections import Counter, OrderedDict

ROOT = "/home/claude/mma"
HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = OrderedDict([
    ("v0.2", "results/wbm_iem/Harvey_1_03d_iem_results_v0.2.json"),
    ("v0.2b", "results/wbm_iem/Harvey_1_03d_iem_results_v0.2b.json"),
    ("v0.3", "results/wbm_iem/Harvey_1_03d_iem_results_v0.3.json"),
])
TAUS = [0.0, 0.001, 0.01, 0.05, 0.10]
TOL = 1e-6


def expected_label(text):
    if "Incre" in text:
        return "Increased"
    if "Decre" in text:
        return "Decreased"
    return "Unchanged"


def is_finite_number(v):
    return isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v)


def call_abs(h, d):
    diff = d - h
    if diff > TOL:
        return "Increased"
    if diff < -TOL:
        return "Decreased"
    return "Unchanged"


def call_rel(h, d, tau):
    if tau == 0.0:
        return call_abs(h, d)
    scale = max(abs(h), abs(d))
    r = (d - h) / scale if scale > 0 else 0.0
    if r > tau:
        return "Increased"
    if r < -tau:
        return "Decreased"
    return "Unchanged"


def opposite(a, b):
    return {a, b} == {"Increased", "Decreased"}


def load_rows(protocol, results, label):
    """One row per protocol biomarker, matched by (iem, call_index) and by position; names are checked."""
    by_key = {}
    for r in results:
        key = (r["iem"], r["call_index"])
        if key in by_key:
            raise SystemExit(f"{label}: duplicate record {key}")
        by_key[key] = r
    problems = []
    if set(by_key) != {(p["iem"], p["call_index"]) for p in protocol}:
        problems.append(f"{label}: record keys differ from protocol keys")
    rows = []
    for p in protocol:
        rec = by_key.get((p["iem"], p["call_index"]))
        bms = rec["biomarkers"] if rec else []
        if rec and len(bms) != len(p["biomarkers"]):
            problems.append(f"{label} {p['iem']}: {len(bms)} biomarker records for {len(p['biomarkers'])} protocol biomarkers")
        for k, (rxn, text) in enumerate(p["biomarkers"]):
            b = bms[k] if k < len(bms) else None
            if b is not None and b["reaction"] != rxn:
                problems.append(f"{label} {p['iem']}: position {k} is {b['reaction']}, protocol {rxn}")
            exp = expected_label(text)
            row = {"iem": p["iem"], "call_index": p["call_index"], "pos": k, "reaction": rxn, "label": text,
                   "expected": exp, "iem_status": rec["status"] if rec else "missing"}
            if b is None:
                row.update(scored=False, call="NA", stored_predicted="absent", healthy=None, disease=None,
                           status_healthy="missing", status_disease="missing", stored_expected=None, stored_correct=None)
            else:
                h, d = b["healthy"], b["disease"]
                ok = (exp in ("Increased", "Decreased") and b["status_healthy"] == "Optimal"
                      and b["status_disease"] == "Optimal" and is_finite_number(h) and is_finite_number(d))
                both_values = is_finite_number(h) and is_finite_number(d) and b["status_healthy"] == "Optimal" \
                    and b["status_disease"] == "Optimal"
                row.update(scored=ok, call=call_abs(h, d) if both_values else "NA", stored_predicted=b["predicted"],
                           healthy=h, disease=d, status_healthy=b["status_healthy"], status_disease=b["status_disease"],
                           stored_expected=b["expected"], stored_correct=b["correct"])
            rows.append(row)
    return rows, problems


def summarise(rows):
    scored = [r for r in rows if r["scored"]]
    correct = [r for r in scored if r["call"] == r["expected"]]
    errs = Counter()
    for r in scored:
        if r["call"] != r["expected"]:
            errs["no change predicted" if r["call"] == "Unchanged" else
                 "opposite direction" if opposite(r["call"], r["expected"]) else "other"] += 1
    per_iem = OrderedDict()
    for r in rows:
        e = per_iem.setdefault(r["iem"], {"n": 0, "scored": 0, "correct": 0})
        e["n"] += 1
        e["scored"] += r["scored"]
        e["correct"] += r["scored"] and r["call"] == r["expected"]
    all_ok = [i for i, e in per_iem.items() if e["scored"] == e["n"] and e["correct"] == e["n"]]
    all_ok_among_scored = [i for i, e in per_iem.items() if e["scored"] and e["correct"] == e["scored"]]
    unscored = [{"iem": r["iem"], "reaction": r["reaction"], "expected": r["expected"], "status_healthy": r["status_healthy"],
                 "status_disease": r["status_disease"], "healthy": r["healthy"], "disease": r["disease"]}
                for r in rows if not r["scored"]]
    # consistency with the stored fields
    stored_pred_mismatch = [(r["iem"], r["reaction"], r["call"], r["stored_predicted"]) for r in rows
                            if r["call"] != "NA" and r["call"] != r["stored_predicted"]]
    stored_exp_mismatch = [(r["iem"], r["reaction"], r["expected"], r["stored_expected"]) for r in rows
                           if r["stored_expected"] is not None and r["stored_expected"] != r["expected"]]
    stored_correct_mismatch = [(r["iem"], r["reaction"]) for r in rows
                               if r["scored"] and r["stored_correct"] is not (r["call"] == r["expected"])]
    stored_correct_on_unscored = [(r["iem"], r["reaction"], r["stored_correct"]) for r in rows
                                  if not r["scored"] and r["stored_correct"] is not None]
    by_exp = Counter((r["expected"], r["call"]) for r in scored)
    return OrderedDict(
        n_protocol_biomarkers=len(rows), n_scored=len(scored), n_correct=len(correct),
        accuracy=len(correct) / len(scored) if scored else None, coverage=len(scored) / len(rows),
        errors=dict(errs), expected_vs_call={f"{a}->{b}": n for (a, b), n in sorted(by_exp.items())},
        n_iems=len(per_iem), n_iems_every_biomarker_scored_and_correct=len(all_ok),
        n_iems_all_correct_among_scored=len(all_ok_among_scored),
        iem_status=dict(Counter(r["iem_status"] for r in {(x["iem"]): x for x in rows}.values())),
        unscored=unscored, stored_predicted_mismatch=stored_pred_mismatch,
        stored_expected_mismatch=stored_exp_mismatch, stored_correct_mismatch=stored_correct_mismatch,
        stored_correct_set_on_unscored=stored_correct_on_unscored,
        iems_not_all_correct=[i for i in per_iem if i not in all_ok])


def changes(rows_a, rows_b):
    out = []
    for a, b in zip(rows_a, rows_b):
        assert (a["iem"], a["pos"], a["reaction"]) == (b["iem"], b["pos"], b["reaction"])
        if a["call"] != b["call"] or a["stored_predicted"] != b["stored_predicted"]:
            out.append({"iem": a["iem"], "reaction": a["reaction"], "expected": a["expected"],
                        "from": a["call"], "to": b["call"], "stored_from": a["stored_predicted"],
                        "stored_to": b["stored_predicted"], "healthy_from": a["healthy"], "disease_from": a["disease"],
                        "healthy_to": b["healthy"], "disease_to": b["disease"]})
    return out


def effect_size(rows):
    scored = [r for r in rows if r["scored"]]
    out = OrderedDict(n_scored=len(scored))
    for tau in TAUS:
        out[f"tau={tau:g}"] = sum(call_rel(r["healthy"], r["disease"], tau) == r["expected"] for r in scored)
    small = []
    for r in scored:
        if r["call"] == r["expected"]:
            scale = max(abs(r["healthy"]), abs(r["disease"]))
            rel = abs(r["disease"] - r["healthy"]) / scale if scale > 0 else 0.0
            if rel <= 0.10:
                small.append({"iem": r["iem"], "reaction": r["reaction"], "healthy": r["healthy"], "disease": r["disease"],
                              "relative_change": rel})
    out["protocol_correct_with_relative_change_le_10pct"] = sorted(small, key=lambda x: x["relative_change"])
    return out


def main():
    protocol = json.load(open(os.path.join(ROOT, "data/iem/iem_protocol_v0.2.json")))
    report = OrderedDict(protocol=OrderedDict(
        n_iems=len(protocol), n_biomarkers=sum(len(p["biomarkers"]) for p in protocol),
        expected=dict(Counter(expected_label(t) for p in protocol for _, t in p["biomarkers"]))))
    all_rows, problems = OrderedDict(), []
    for label, path in RUNS.items():
        results = json.load(open(os.path.join(ROOT, path)))
        rows, prob = load_rows(protocol, results, label)
        problems += prob
        all_rows[label] = rows
        report[label] = summarise(rows)
        report[label]["effect_size"] = effect_size(rows)
    report["structure_problems"] = problems
    labels = list(RUNS)
    report["prediction_changes"] = OrderedDict()
    for a, b in [("v0.2", "v0.2b"), ("v0.2b", "v0.3"), ("v0.2", "v0.3")]:
        report["prediction_changes"][f"{a}->{b}"] = changes(all_rows[a], all_rows[b])
    with open(os.path.join(HERE, "scores.json"), "w") as fh:
        json.dump(report, fh, indent=1, allow_nan=False)
        fh.write("\n")
    # per-biomarker table for the report
    with open(os.path.join(HERE, "per_biomarker_calls.tsv"), "w") as fh:
        head = ["iem", "pos", "reaction", "expected"]
        for l in labels:
            head += [f"{l}_healthy", f"{l}_disease", f"{l}_call", f"{l}_scored"]
        fh.write("\t".join(head) + "\n")
        for i in range(len(all_rows[labels[0]])):
            r0 = all_rows[labels[0]][i]
            line = [r0["iem"], str(r0["pos"]), r0["reaction"], r0["expected"]]
            for l in labels:
                r = all_rows[l][i]
                line += [repr(r["healthy"]), repr(r["disease"]), r["call"], str(int(r["scored"]))]
            fh.write("\t".join(line) + "\n")
    for l in labels:
        s = report[l]
        print(f"{l}: scored {s['n_scored']}/{s['n_protocol_biomarkers']} correct {s['n_correct']} "
              f"acc {s['accuracy']:.4f} errors {s['errors']} IEMs all scored+correct {s['n_iems_every_biomarker_scored_and_correct']} "
              f"(all correct among scored {s['n_iems_all_correct_among_scored']}) status {s['iem_status']}")
        print(f"   unscored: {[(u['iem'], u['reaction'], u['status_healthy'], u['status_disease']) for u in s['unscored']]}")
        print(f"   mismatches stored pred {len(s['stored_predicted_mismatch'])} exp {len(s['stored_expected_mismatch'])} "
              f"correct {len(s['stored_correct_mismatch'])} correct-on-unscored {s['stored_correct_set_on_unscored']}")
        es = s["effect_size"]
        print("   effect size: " + ", ".join(f"{k} {v}" for k, v in es.items() if k.startswith("tau")))
    print("structure problems:", problems)
    for k, v in report["prediction_changes"].items():
        print(f"{k}: {len(v)} changes")
        for c in v:
            print(f"   {c['iem']:7s} {c['reaction']:20s} exp {c['expected']:9s} {c['from']} -> {c['to']} "
                  f"(stored {c['stored_from']} -> {c['stored_to']}) h {c['healthy_from']!r}->{c['healthy_to']!r} "
                  f"d {c['disease_from']!r}->{c['disease_to']!r}")


if __name__ == "__main__":
    sys.exit(main())
