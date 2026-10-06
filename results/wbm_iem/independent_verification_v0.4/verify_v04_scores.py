"""Independent recomputation of the v0.4 scores, error anatomy, shared errors, effect sizes, flux-range rules
and the post hoc 'indeterminate' figures, from the raw results files and the protocol only.

Imports nothing from gembench or the project's scripts. Writes scores_v04.json and per_biomarker_v04.tsv here.
"""
import sys
from collections import Counter, OrderedDict

from vcommon import (DIRS, FILES, TOL, biomarker, call_protocol, call_rel, dump, finite, index_results, load_json,
                     max_ok, min_ok, protocol_rows, rel_change, z, HERE, PROTOCOL)
import json
import os

TAUS = [0.001, 0.01, 0.05, 0.10]


def fluid(rxn):
    if rxn.endswith("[u]"):
        return "urine"
    if rxn.endswith("[bc]"):
        return "blood"
    if rxn.endswith("[csf]"):
        return "csf"
    return "other"


def score_rows(rows, by, zero=True):
    out = []
    for row in rows:
        b = biomarker(by, row)
        r = OrderedDict(row)
        h, d = b.get("healthy"), b.get("disease")
        r.update(h=h, d=d, sh=b.get("status_healthy"), sd=b.get("status_disease"),
                 stored_pred=b.get("predicted"), stored_exp=b.get("expected"), stored_correct=b.get("correct"))
        ok = max_ok(b)
        r["values_ok"] = ok
        r["scored"] = ok and row["expected"] in DIRS
        r["call"] = call_protocol(h, d, zero) if ok else "NA"
        r["call_nozero"] = call_protocol(h, d, False) if ok else "NA"
        r["correct"] = bool(r["scored"] and r["call"] == row["expected"])
        cat = sub = None
        if r["scored"] and not r["correct"]:
            if r["call"] == "Unchanged":
                cat = "no_change"
                if abs(d - h) <= TOL and h > TOL:
                    sub = "capped"
                elif abs(h) <= TOL and abs(d) <= TOL:
                    sub = "zero"
                else:
                    sub = "other"
            elif r["call"] in DIRS:
                cat = "opposite"
        r["err"], r["err_sub"] = cat, sub
        r["rel"] = rel_change(h, d) if ok else None
        r["rel_zeroed"] = rel_change(h, d, zero=True) if ok else None
        out.append(r)
    return out


def summary(rows):
    scored = [r for r in rows if r["scored"]]
    correct = [r for r in scored if r["correct"]]
    errs = Counter(r["err"] for r in scored if not r["correct"])
    subs = Counter(r["err_sub"] for r in scored if r["err"] == "no_change")
    per_iem = OrderedDict()
    for r in rows:
        e = per_iem.setdefault(r["iem"], dict(n=0, scored=0, correct=0))
        e["n"] += 1
        e["scored"] += r["scored"]
        e["correct"] += r["correct"]
    fully = [i for i, e in per_iem.items() if e["n"] > 0 and e["scored"] == e["n"] and e["correct"] == e["n"]]
    stored_pred_mismatch = [(r["iem"], r["reaction"], r["call"], r["stored_pred"]) for r in rows
                            if r["values_ok"] and r["call"] != r["stored_pred"]]
    stored_exp_mismatch = [(r["iem"], r["reaction"], r["expected"], r["stored_exp"]) for r in rows
                           if r["stored_exp"] != r["expected"]]
    stored_correct_mismatch = [(r["iem"], r["reaction"], r["correct"], r["stored_correct"]) for r in rows
                               if r["scored"] and r["stored_correct"] is not r["correct"]]
    zero_vs_raw = [(r["iem"], r["reaction"], r["h"], r["d"], r["call"], r["call_nozero"]) for r in rows
                   if r["values_ok"] and r["call"] != r["call_nozero"]]
    return OrderedDict(
        n_biomarkers=len(rows), n_directional=sum(r["expected"] in DIRS for r in rows), n_scored=len(scored),
        n_correct=len(correct), accuracy=len(correct) / len(scored) if scored else None,
        errors=dict(errs), no_change_breakdown=dict(subs),
        n_iems=len(per_iem), n_iems_fully_correct=len(fully),
        iems_not_fully_correct=[i for i in per_iem if i not in fully],
        unscored=[(r["iem"], r["reaction"], r["sh"], r["sd"], r["h"], r["d"]) for r in rows if not r["scored"]],
        opposite_errors=[(r["iem"], r["reaction"], r["expected"], r["call"], r["h"], r["d"]) for r in scored if r["err"] == "opposite"],
        no_change_errors=[(r["iem"], r["reaction"], r["expected"], r["err_sub"], r["h"], r["d"]) for r in scored if r["err"] == "no_change"],
        stored_predicted_mismatch=stored_pred_mismatch, stored_expected_mismatch=stored_exp_mismatch,
        stored_correct_mismatch=stored_correct_mismatch, zeroing_changes_call=zero_vs_raw)


def effect_sizes(rows):
    scored = [r for r in rows if r["scored"]]
    out = OrderedDict(n_scored=len(scored), protocol=sum(r["correct"] for r in scored))
    for zero in (False, True):
        tag = "zeroed" if zero else "raw"
        for tau in TAUS:
            calls = [(r, call_rel(r["h"], r["d"], tau, zero)) for r in scored]
            corr = [r for r, c in calls if c == r["expected"]]
            lost = [(r["iem"], r["reaction"]) for r, c in calls if r["correct"] and c != r["expected"]]
            gained = [(r["iem"], r["reaction"]) for r, c in calls if not r["correct"] and c == r["expected"]]
            # variant: protocol call AND |rel| > tau (both thresholds)
            both = sum(r["correct"] and abs(rel_change(r["h"], r["d"], zero)) > tau for r in scored)
            out[f"tau={tau:g} ({tag})"] = dict(correct=len(corr), lost=lost, gained=gained, protocol_and_rel=both)
    small = sorted([dict(iem=r["iem"], reaction=r["reaction"], h=r["h"], d=r["d"], rel=abs(r["rel"]))
                    for r in scored if r["correct"] and abs(r["rel"]) < 0.05], key=lambda x: x["rel"])
    out["correct_below_5pct"] = small
    out["n_correct_below_5pct"] = len(small)
    out["n_correct_below_0.1pct"] = sum(x["rel"] < 0.001 for x in small)
    out["correct_5_to_10pct"] = sorted([dict(iem=r["iem"], reaction=r["reaction"], h=r["h"], d=r["d"], rel=abs(r["rel"]))
                                        for r in scored if r["correct"] and 0.05 <= abs(r["rel"]) <= 0.10], key=lambda x: x["rel"])
    return out


def sign(a, b, zero=True):
    if zero:
        a, b = z(a), z(b)
    diff = b - a
    return 1 if diff > TOL else (-1 if diff < -TOL else 0)


def r1(hmax, dmax, hmin, dmin, zero=True):
    smax, smin = sign(hmax, dmax, zero), sign(hmin, dmin, zero)
    up = smax > 0 or smin > 0
    down = smax < 0 or smin < 0
    if up and down:
        return "Conflicting"
    if up:
        return "Increased"
    if down:
        return "Decreased"
    return "Unchanged"


def r2(prot_call, hmin, dmin, zero=True):
    if prot_call != "Unchanged":
        return prot_call
    s = sign(hmin, dmin, zero)
    return "Increased" if s > 0 else ("Decreased" if s < 0 else "Unchanged")


def range_rules(rows_max, by_min, zero=True):
    sc = []
    n_dir = 0
    excluded = []
    for r in rows_max:
        if r["expected"] not in DIRS:
            continue
        n_dir += 1
        bm = biomarker(by_min, r)
        if not (r["values_ok"] and min_ok(bm)):
            excluded.append((r["iem"], r["reaction"], r["sh"], r["sd"], bm.get("status_healthy_min"), bm.get("status_disease_min")))
            continue
        hmn, dmn = bm["healthy_min"], bm["disease_min"]
        prot = call_protocol(r["h"], r["d"], zero)
        c1 = r1(r["h"], r["d"], hmn, dmn, zero)
        c2 = r2(prot, hmn, dmn, zero)
        sc.append(dict(iem=r["iem"], reaction=r["reaction"], expected=r["expected"], hmax=r["h"], dmax=r["d"],
                       hmin=hmn, dmin=dmn, prot=prot, r1=c1, r2=c2, capped=r["err_sub"] == "capped",
                       err_sub=r["err_sub"], smax=sign(r["h"], r["d"], zero), smin=sign(hmn, dmn, zero),
                       fluid=fluid(r["reaction"])))
    out = OrderedDict(n_directional=n_dir, n_range_scored=len(sc), excluded=excluded)
    out["maxima_only_correct"] = sum(x["prot"] == x["expected"] for x in sc)
    for rule in ("r1", "r2"):
        tr = Counter()
        lists = {"right_to_wrong": [], "wrong_to_right": [], "wrong_to_other_wrong": []}
        for x in sc:
            a, b = x["prot"] == x["expected"], x[rule] == x["expected"]
            if a and not b:
                tr["right_to_wrong"] += 1
                lists["right_to_wrong"].append(x)
            elif not a and b:
                tr["wrong_to_right"] += 1
                lists["wrong_to_right"].append(x)
            elif not a and not b and x[rule] != x["prot"]:
                tr["wrong_to_other_wrong"] += 1
                lists["wrong_to_other_wrong"].append(x)
            elif a and b:
                tr["right_to_right"] += 1
            else:
                tr["wrong_same_wrong"] += 1
        out[rule] = OrderedDict(
            correct=sum(x[rule] == x["expected"] for x in sc), transitions=dict(tr),
            calls=dict(Counter(x[rule] for x in sc)),
            right_to_wrong_new_calls=dict(Counter(x[rule] for x in lists["right_to_wrong"])),
            right_to_wrong_pattern=dict(Counter(
                f"smax={x['smax']} smin={x['smin']} " + ("dmin_zero" if abs(x["dmin"]) <= TOL else "dmin_nonzero") + " " +
                ("hmin_pos" if x["hmin"] > TOL else "hmin_not_pos") for x in lists["right_to_wrong"])),
            right_to_wrong=[(x["iem"], x["reaction"], x["expected"], x["prot"], x[rule], x["hmax"], x["dmax"], x["hmin"], x["dmin"])
                            for x in lists["right_to_wrong"]],
            wrong_to_right=[(x["iem"], x["reaction"], x["expected"], x["prot"], x[rule]) for x in lists["wrong_to_right"]],
            wrong_to_other_wrong=[(x["iem"], x["reaction"], x["expected"], x["prot"], x[rule], x["hmax"], x["dmax"], x["hmin"], x["dmin"], x["capped"])
                                  for x in lists["wrong_to_other_wrong"]],
            conflicting_total=sum(x[rule] == "Conflicting" for x in sc),
            by_fluid={f: dict(n=sum(x["fluid"] == f for x in sc),
                              maxima_only=sum(x["fluid"] == f and x["prot"] == x["expected"] for x in sc),
                              rule=sum(x["fluid"] == f and x[rule] == x["expected"] for x in sc))
                      for f in sorted({x["fluid"] for x in sc})})
    capped = [x for x in sc if x["capped"]]
    out["capped"] = OrderedDict(
        n=len(capped),
        zero_minima_both=sum(abs(x["hmin"]) <= TOL and abs(x["dmin"]) <= TOL for x in capped),
        differing_minima=sum(abs(z(x["dmin"]) - z(x["hmin"])) > TOL for x in capped),
        differing=[(x["iem"], x["reaction"], x["expected"], x["hmax"], x["dmax"], x["hmin"], x["dmin"], x["r2"]) for x in capped
                   if abs(z(x["dmin"]) - z(x["hmin"])) > TOL],
        nonzero_equal_minima=[(x["iem"], x["reaction"], x["hmin"], x["dmin"]) for x in capped
                              if not (abs(x["hmin"]) <= TOL and abs(x["dmin"]) <= TOL) and abs(z(x["dmin"]) - z(x["hmin"])) <= TOL],
        r1_correct=sum(x["r1"] == x["expected"] for x in capped), r2_correct=sum(x["r2"] == x["expected"] for x in capped))
    zero_err = [x for x in sc if x["err_sub"] == "zero"]
    out["zero_no_change_errors"] = [(x["iem"], x["reaction"], x["hmin"], x["dmin"], x["r1"], x["r2"]) for x in zero_err]
    return out, sc


def main():
    protocol = load_json(os.path.relpath(PROTOCOL, "/home/claude/mma"))
    rows = protocol_rows(protocol)
    rep = OrderedDict()
    rep["protocol"] = dict(n_iems=len(protocol), n_biomarkers=len(rows),
                           expected=dict(Counter(r["expected"] for r in rows)),
                           iem_names_unique=len({p["iem"] for p in protocol}) == len(protocol))
    by = {k: index_results(load_json(v)) for k, v in FILES.items()}
    for k, b in by.items():
        if set(b) != {(p["iem"], p["call_index"]) for p in protocol}:
            raise SystemExit(f"{k}: record keys differ from protocol")
    scored = {}
    for k in ("harvey_max_v0.2", "harvey_max_v0.2b", "harvey_max_v0.3", "harvetta_max_A"):
        scored[k] = score_rows(rows, by[k])
        rep[k] = summary(scored[k])
    # ---- consistency of the four zeroing variants (raw vs zeroed) for anatomy
    # ---- shared errors between Harvey v0.3 and Harvetta A
    H, T = scored["harvey_max_v0.3"], scored["harvetta_max_A"]
    pairs = list(zip(H, T))
    for a, b in pairs:
        assert (a["iem"], a["pos"], a["reaction"]) == (b["iem"], b["pos"], b["reaction"])
    key = lambda r: (r["iem"], r["reaction"])
    sh = OrderedDict()
    sh["opposite_both"] = [key(a) for a, b in pairs if a["err"] == "opposite" and b["err"] == "opposite"]
    sh["opposite_harvey_only"] = [(key(a), a["err"] or "correct" if b["err"] is None else b["err"], b["err"] or ("correct" if b["correct"] else "unscored"))
                                  for a, b in pairs if a["err"] == "opposite" and b["err"] != "opposite"]
    sh["opposite_harvetta_only"] = [(key(b), a["err"] or ("correct" if a["correct"] else "unscored"), a["err_sub"])
                                    for a, b in pairs if b["err"] == "opposite" and a["err"] != "opposite"]
    sh["no_change_both"] = [key(a) for a, b in pairs if a["err"] == "no_change" and b["err"] == "no_change"]
    sh["capped_both"] = [key(a) for a, b in pairs if a["err_sub"] == "capped" and b["err_sub"] == "capped"]
    sh["capped_harvetta_not_capped_harvey"] = [(key(b), a["err"], a["err_sub"], a["correct"]) for a, b in pairs
                                               if b["err_sub"] == "capped" and a["err_sub"] != "capped"]
    sh["capped_harvey_not_capped_harvetta"] = [(key(a), b["err"], b["err_sub"], b["correct"]) for a, b in pairs
                                               if a["err_sub"] == "capped" and b["err_sub"] != "capped"]
    sh["no_change_harvetta_only"] = [(key(b), a["err"], a["correct"]) for a, b in pairs if b["err"] == "no_change" and a["err"] != "no_change"]
    sh["no_change_harvey_only"] = [(key(a), b["err"], b["correct"]) for a, b in pairs if a["err"] == "no_change" and b["err"] != "no_change"]
    sh["counts"] = {k: len(v) for k, v in sh.items() if isinstance(v, list)}
    # union-based errors present in both, regardless of type
    sh["any_error_both"] = sum(1 for a, b in pairs if a["err"] and b["err"])
    sh["same_call_both_models"] = sum(1 for a, b in pairs if a["call"] == b["call"])
    rep["shared_errors_harvey_v0.3_vs_harvetta"] = sh
    # ---- values of the named biomarkers in both models
    named = [("ASNSD", "DM_asn_L[bc]"), ("ASNSD", "DM_asn_L[csf]"), ("MMA", "DM_crn[bc]"), ("BTD", "EX_3hpp[u]"),
             ("MMA", "DM_3hpp[bc]"), ("GACR", "DM_gln_L[bc]"), ("GACR", "EX_glu_L[u]"), ("HPC", "EX_C05770[u]"),
             ("HPC", "DM_C05770[bc]"), ("HPC", "DM_C05770[csf]"), ("LNS", "DM_fol[bc]"), ("PC", "DM_glc_D[bc]"),
             ("CYP21D", "DM_aldstrn[bc]"), ("CYP21D", "DM_crtsl[bc]")]
    vals = OrderedDict()
    for a, b in pairs:
        if any(a["iem"] == i and (a["reaction"] == x) for i, x in named) or a["err"] or b["err"]:
            vals[f"{a['iem']} {a['reaction']}"] = dict(
                expected=a["expected"], harvey=dict(h=a["h"], d=a["d"], call=a["call"], err=a["err"], sub=a["err_sub"], rel=a["rel"]),
                harvetta=dict(h=b["h"], d=b["d"], call=b["call"], err=b["err"], sub=b["err_sub"], rel=b["rel"]))
    rep["values_errors_and_named"] = vals
    # ---- effect sizes
    rep["effect_size"] = OrderedDict((k, effect_sizes(scored[k])) for k in ("harvey_max_v0.3", "harvetta_max_A", "harvey_max_v0.2", "harvey_max_v0.2b"))
    # ---- flux-range rules
    rng = OrderedDict()
    for model, kmax, kmin in (("Harvey", "harvey_max_v0.3", "harvey_min_B"), ("Harvetta", "harvetta_max_A", "harvetta_min_C")):
        res, sc = range_rules(scored[kmax], by[kmin], zero=True)
        res_raw, _ = range_rules(scored[kmax], by[kmin], zero=False)
        res["raw_variant_same_counts"] = (res_raw["r1"]["correct"], res_raw["r2"]["correct"], res_raw["maxima_only_correct"]) == \
                                         (res["r1"]["correct"], res["r2"]["correct"], res["maxima_only_correct"])
        res["raw_variant_counts"] = dict(r1=res_raw["r1"]["correct"], r2=res_raw["r2"]["correct"], maxima_only=res_raw["maxima_only_correct"])
        # min-file sanity: the min files carry no maxima; check that every min record is 'not_run' for maxima
        mins = by[kmin]
        res["min_file_maxima_statuses"] = dict(Counter((b["status_healthy"], b["status_disease"]) for r in mins.values() for b in r["biomarkers"]))
        res["min_file_min_statuses"] = dict(Counter((b["status_healthy_min"], b["status_disease_min"]) for r in mins.values() for b in r["biomarkers"]))
        res["n_min_pairs_optimal_finite"] = sum(min_ok(b) for r in mins.values() for b in r["biomarkers"])
        rng[model] = res
    rep["flux_ranges"] = rng
    # ---- post hoc indeterminate
    ph = OrderedDict()
    for k in ("harvey_max_v0.3", "harvetta_max_A"):
        s = [r for r in scored[k] if r["scored"]]
        capped = [r for r in s if r["err_sub"] == "capped"]
        eqpos = [r for r in s if abs(r["d"] - r["h"]) <= TOL and r["h"] > TOL]
        rest = [r for r in s if r["err_sub"] != "capped"]
        ph[k] = dict(n_capped=len(capped), n_equal_positive_maxima_among_scored=len(eqpos), n_rest=len(rest),
                     correct=sum(r["correct"] for r in rest), accuracy=sum(r["correct"] for r in rest) / len(rest))
    rep["post_hoc_indeterminate"] = ph
    # ---- Harvey v0.2 -> v0.2b changes (Paper 2 table)
    ch = []
    for a, b in zip(scored["harvey_max_v0.2"], scored["harvey_max_v0.2b"]):
        if a["call"] != b["call"]:
            ch.append(dict(iem=a["iem"], reaction=a["reaction"], expected=a["expected"], v02=(a["h"], a["d"], a["call"]),
                           v02b=(b["h"], b["d"], b["call"]), rel_v02=a["rel"], rel_v02b=b["rel"],
                           effect=("right->wrong" if a["correct"] and not b["correct"] else
                                   "wrong->right" if b["correct"] and not a["correct"] else "other")))
    rep["harvey_v0.2_to_v0.2b_changes"] = ch
    ch2 = [(a["iem"], a["reaction"]) for a, b in zip(scored["harvey_max_v0.2b"], scored["harvey_max_v0.3"]) if a["call"] != b["call"]]
    rep["harvey_v0.2b_to_v0.3_call_changes"] = ch2
    dump("scores_v04.json", rep)
    with open(os.path.join(HERE, "per_biomarker_v04.tsv"), "w") as fh:
        fh.write("iem\tpos\treaction\texpected\tharvey_h\tharvey_d\tharvey_call\tharvey_err\tharvetta_h\tharvetta_d\tharvetta_call\tharvetta_err\n")
        for a, b in pairs:
            fh.write("\t".join(map(str, [a["iem"], a["pos"], a["reaction"], a["expected"], a["h"], a["d"], a["call"],
                                         a["err_sub"] or a["err"] or ("ok" if a["correct"] else "unscored"),
                                         b["h"], b["d"], b["call"], b["err_sub"] or b["err"] or ("ok" if b["correct"] else "unscored")])) + "\n")
    # ---- print
    for k in ("harvey_max_v0.2", "harvey_max_v0.2b", "harvey_max_v0.3", "harvetta_max_A"):
        s = rep[k]
        print(f"{k}: scored {s['n_scored']}/{s['n_biomarkers']} correct {s['n_correct']} acc {s['accuracy']:.4f} "
              f"errors {s['errors']} no-change {s['no_change_breakdown']} fully correct {s['n_iems_fully_correct']}/{s['n_iems']}")
        print(f"   unscored {s['unscored']}")
        print(f"   stored mismatches pred {len(s['stored_predicted_mismatch'])} exp {len(s['stored_expected_mismatch'])} "
              f"correct {len(s['stored_correct_mismatch'])}; zeroing changes call: {s['zeroing_changes_call']}")
    print("shared:", sh["counts"])
    for kk in ("opposite_both", "opposite_harvey_only", "opposite_harvetta_only", "capped_harvetta_not_capped_harvey",
               "capped_harvey_not_capped_harvetta", "no_change_harvetta_only", "no_change_harvey_only"):
        print(f"  {kk}: {sh[kk]}")
    for k, e in rep["effect_size"].items():
        print(f"effect {k}: protocol {e['protocol']}, " + ", ".join(f"{t}: {v['correct']} (both {v['protocol_and_rel']}, lost {len(v['lost'])}, gained {len(v['gained'])})"
                                                           for t, v in e.items() if t.startswith("tau")))
        print(f"   correct below 5%: {e['n_correct_below_5pct']}; below 0.1%: {e['n_correct_below_0.1pct']}")
        for x in e["correct_below_5pct"]:
            print(f"      {x['iem']:6s} {x['reaction']:18s} h={x['h']:.8g} d={x['d']:.8g} rel={x['rel']:.3g}")
        for x in e["correct_5_to_10pct"]:
            print(f"      (5-10%) {x['iem']:6s} {x['reaction']:18s} h={x['h']:.8g} d={x['d']:.8g} rel={x['rel']:.3g}")
    for m, res in rng.items():
        print(f"ranges {m}: directional {res['n_directional']} scored {res['n_range_scored']} max-only {res['maxima_only_correct']} "
              f"R1 {res['r1']['correct']} {res['r1']['transitions']} R2 {res['r2']['correct']} {res['r2']['transitions']} "
              f"raw same {res['raw_variant_same_counts']}")
        print(f"   R1 right->wrong new calls {res['r1']['right_to_wrong_new_calls']} pattern {res['r1']['right_to_wrong_pattern']}")
        print(f"   R2 wrong->other {res['r2']['wrong_to_other_wrong']}")
        print(f"   capped {res['capped']}")
        print(f"   zero errors {res['zero_no_change_errors']}; excluded {res['excluded']}")
        print(f"   min statuses {res['min_file_min_statuses']} optimal-finite min pairs {res['n_min_pairs_optimal_finite']}")
        print(f"   by fluid R1 {res['r1']['by_fluid']}")
    print("post hoc:", ph)
    print("v0.2->v0.2b:", json.dumps(ch, default=str, indent=0))
    print("v0.2b->v0.3 call changes:", ch2)


if __name__ == "__main__":
    sys.exit(main())
