"""Collect the verifier's recomputed numbers (from scores_v04.json, matlab_v04.json, bounds_v04.json,
provenance_v04.json and caps_v04.json in this folder) into recomputed_numbers_v04.json, with each documented value
next to the recomputed one. Run after the verify_v04_*.py scripts. No gembench or project-script imports.
"""
import json
import os
import sys
from collections import OrderedDict

HERE = os.path.dirname(os.path.abspath(__file__))


def j(name):
    with open(os.path.join(HERE, name)) as fh:
        return json.load(fh)


def main():
    s, m, b, p, c = j("scores_v04.json"), j("matlab_v04.json"), j("bounds_v04.json"), j("provenance_v04.json"), j("caps_v04.json")
    H, T = s["harvey_max_v0.3"], s["harvetta_max_A"]
    eh, et = s["effect_size"]["harvey_max_v0.3"], s["effect_size"]["harvetta_max_A"]
    rh, rt = s["flux_ranges"]["Harvey"], s["flux_ranges"]["Harvetta"]
    mh, mt = m["Harvey"], m["Harvetta"]
    bh, bt = b["Harvey"], b["Harvetta"]
    g = lambda grp, k: bt["setup_vs_stored"]["groups"][grp][k]
    gh = lambda grp, k: bh["setup_vs_stored"]["groups"][grp][k]
    named = s["values_errors_and_named"]
    out = OrderedDict()
    out["completeness_provenance"] = {k: dict(n_records=v["n_records"], statuses=v["statuses"], n_fingerprints=v["n_fingerprints"],
                                              fingerprint=v["fingerprint"], fp_equals_sha256_provenance=v["fingerprint_equals_sha256_of_provenance"],
                                              senses=v["senses"], lp_bounds_sha256=v["lp_bounds_sha256"],
                                              model_sha_ok=v["model_sha256_matches_file"], protocol_sha_ok=v["protocol_sha256_matches_file"],
                                              wall_h=round(v["session_wall_h"], 3), sum_time_h=round(v["sum_time_h"], 3),
                                              inferred_start_utc=v["inferred_start_utc"], end_utc=v["summary_mtime_utc"])
                                      for k, v in p["runs"].items()}
    out["completeness_provenance"]["B_lp_bounds_equals_v0.3"] = p["runs"]["B_harvey_min"]["lp_bounds_equals_v0.3"]
    out["freeze_and_time_order"] = dict(
        freeze_valid=p["freeze"]["freeze_verify_tool"]["valid"], fingerprint_recomputed=p["freeze"]["recomputed"],
        fingerprint_match=p["freeze"]["match"], n_files_ok=p["freeze"]["n_sha_ok"], created_at=p["freeze"]["created_at"],
        freeze_head=p["freeze"]["repository"]["head"], freeze_dirty=p["freeze"]["repository"]["dirty"],
        plan_commit=p["git"]["commits"]["plan_md"][-1], freeze_commit=p["git"]["commits"]["freeze_manifest"][-1],
        first_in_git=p["git"]["first_appearance_any_ref"],
        plan_ancestor=p["git"]["plan_ancestor_of_first_results"], freeze_ancestor=p["git"]["freeze_ancestor_of_first_results"],
        frozen_changed_to_HEAD=p["freeze"]["frozen_paths_changed_freeze_commit_to_HEAD"],
        external_record_created_at="2026-10-05T13:39:50.613661Z (Claude project doc studies/2026-10-05-wbm-iem-v0.4-plan-freeze.md, fingerprint f3dde4e2..., read with the Projects tool)",
        matlab_jobs=p["matlab_jobs"])
    out["accuracy_and_anatomy"] = {
        k: dict(scored=v["n_scored"], correct=v["n_correct"], accuracy=round(v["accuracy"], 4), errors=v["errors"],
                no_change=v["no_change_breakdown"], iems_fully_correct=v["n_iems_fully_correct"], unscored=v["unscored"])
        for k, v in (("Harvey_v0.2", s["harvey_max_v0.2"]), ("Harvey_v0.2b", s["harvey_max_v0.2b"]), ("Harvey_v0.3", H), ("Harvetta_v0.4", T))}
    out["shared_errors"] = s["shared_errors_harvey_v0.3_vs_harvetta"]["counts"]
    out["shared_error_lists"] = {k: s["shared_errors_harvey_v0.3_vs_harvetta"][k] for k in ("opposite_both", "opposite_harvetta_only", "opposite_harvey_only")}
    out["named_values"] = {k: named[k] for k in ("ASNSD DM_asn_L[bc]", "MMA DM_crn[bc]", "BTD EX_3hpp[u]", "MMA DM_3hpp[bc]")}
    out["effect_sizes"] = dict(
        Harvey=[eh["protocol"]] + [eh[f"tau={t} (raw)"]["correct"] for t in ("0.001", "0.01", "0.05", "0.1")],
        Harvetta=[et["protocol"]] + [et[f"tau={t} (raw)"]["correct"] for t in ("0.001", "0.01", "0.05", "0.1")],
        zeroed_variant_identical=all(eh[f"tau={t} (raw)"]["correct"] == eh[f"tau={t} (zeroed)"]["correct"] and
                                     et[f"tau={t} (raw)"]["correct"] == et[f"tau={t} (zeroed)"]["correct"] for t in ("0.001", "0.01", "0.05", "0.1")),
        any_gains=any(len(e[k]["gained"]) for e in (eh, et) for k in e if k.startswith("tau")),
        correct_below_5pct=dict(Harvey=eh["n_correct_below_5pct"], Harvetta=et["n_correct_below_5pct"]),
        harvetta_below_5pct=[(x["iem"], x["reaction"], x["rel"]) for x in et["correct_below_5pct"]],
        harvey_below_0p1pct=[(x["iem"], x["reaction"], x["rel"]) for x in eh["correct_below_5pct"] if x["rel"] < 1e-3])
    out["matlab"] = {k: dict(runiem=f"{v['runiem']['correct']}/{v['runiem']['n']} = {v['runiem']['accuracy']:.4f}",
                             matches_stored=v["runiem"]["matches_stored_accuracy"], n_nonfinite=v["n_nonfinite"],
                             nonfinite_side=v["nonfinite_side"], absent=v["absent_or_empty"], n_without_optimum=v["n_without_optimum"],
                             without_optimum_by_iem=v["without_optimum_by_iem"], python_correct_there=v["without_optimum_python_correct"],
                             python_wrong_there=v["without_optimum_python_wrong"], finite=v["n_matlab_finite"],
                             same_call_finite=v["finite_same_call"], matlab_correct_finite=v["finite_matlab_correct"],
                             acc_finite=round(v["accuracy_on_finite"], 4), same_call_all=v["n_same_call_all"],
                             value_pairs_over_1e3=v["values"]["n_pairs_scale_over_1e3"],
                             max_rel_over_1e3=v["values"]["max_rel_scale_over_1e3"],
                             largest=[(x["iem"], x["reaction"], x["side"], x["rel"]) for x in v["values"]["largest"][:3]])
                     for k, v in (("Harvey", mh), ("Harvetta", mt))}
    out["bounds"] = {k: dict(n_rxns=v["n_rxns"]["model"], global_sha256_raw=v["matlab_global_sha256"]["raw"],
                             matches_python_runs=v["matlab_global_sha256"]["raw_matches"], my_global_step=v["my_global_step_counts"],
                             my_global_vs_matlab=v["my_global_from_matlab_setup_vs_matlab_global"],
                             setup_vs_stored_lb_ub=(v["setup_vs_stored"]["n_lb_changed"], v["setup_vs_stored"]["n_ub_changed"]),
                             groups={gname: dict(n=gv["n"], by_bound=gv["by_bound"], ratios=gv["ratios"]) for gname, gv in v["setup_vs_stored"]["groups"].items()},
                             implied_old_gfr=v.get("implied_old_gfr_ml_min"), csf_ratio=v.get("csf_ratio_vs_0.52/0.35"),
                             setup_changes_on_overwritten_entries=v["setup_changes_on_overwritten_entries"],
                             n_bileduct_exits=v["n_bileduct_exits"], diet_crn=v["named_bounds"].get("Diet_EX_crn[d]"))
                     for k, v in (("Harvey", bh), ("Harvetta", bt))}
    out["flux_ranges"] = {k: dict(scored=v["n_range_scored"], maxima_only=v["maxima_only_correct"], R1=v["r1"]["correct"],
                                  R2=v["r2"]["correct"], R1_transitions=v["r1"]["transitions"], R2_transitions=v["r2"]["transitions"],
                                  R1_losses_pattern=v["r1"]["right_to_wrong_pattern"], R1_wrong_to_other=v["r1"]["wrong_to_other_wrong"],
                                  R2_wrong_to_other=v["r2"]["wrong_to_other_wrong"], capped=v["capped"],
                                  by_fluid_R1=v["r1"]["by_fluid"], raw_variant_same=v["raw_variant_same_counts"])
                          for k, v in (("Harvey", rh), ("Harvetta", rt))}
    sp = j("state_pairing_v04.json")
    out["minima_runs_share_states_with_maxima_runs"] = {k: dict(same_iem_reaction_sets=v["same_iem_reaction_sets"], n_iems=v["n_iems"],
                                                                max_rel_diff_vmax=v["max_rel_diff_vmax"]) for k, v in sp.items()}
    out["post_hoc_indeterminate"] = s["post_hoc_indeterminate"]
    out["harvey_v0.2_to_v0.2b"] = s["harvey_v0.2_to_v0.2b_changes"]
    out["caps_single_bound_check"] = {k: dict(n_capped=v["n_capped"], n_equal_to_a_bound_of_same_metabolite=v["n_with_same_metabolite_bound"],
                                              unmatched=[(r["iem"], r["reaction"], r["value"]) for r in v["rows"] if r["n_same_metabolite_bound_matches"] == 0])
                                      for k, v in c.items()}
    out["references"] = OrderedDict([
        ("1_Thiele2020", "Mol Syst Biol 2020;16(5):e8982, doi:10.15252/msb.20198982 - citation correct. Text (Europe PMC full text and Leiden "
                         "repository copy): 'a total of 215 out of 252 (85.3%) ... by Harvetta and 214 out of 252 (84.9%) by Harvey'; Recon3D "
                         "103 of 205 (50.2%). The documents assign 85.3% to Harvey and 84.9% to Harvetta (swapped)."),
        ("2_Heirendt2019", "Nat Protoc 2019;14(3):639-702 - correct (ORBilu record)."),
        ("3_Ebrahim2013", "BMC Syst Biol 2013;7:74, doi:10.1186/1752-0509-7-74 - correct (eScholarship record)."),
        ("4_HuangfuHall2018", "Math Program Comput 2018;10(1):119-142, doi:10.1007/s12532-017-0130-5 - correct; it is the article HiGHS asks users to "
                              "cite, but it describes the parallel dual simplex; the runs used the interior point method (IPX, Schork & Gondzio, "
                              "Math Program Comput 2020;12:603-635)."),
        ("5_Shlomi2009", "Mol Syst Biol 2009;5:263, doi:10.1038/msb.2009.22, PMID 19401675 - citation correct (Europe PMC record). The 'recall 0.10 "
                         "on a clinical database' claim could not be checked against the full text (PMC reCAPTCHA; Europe PMC, Springer and "
                         "archive.org rate-limited; NCBI E-utilities and Crossref blocked)."),
    ])
    with open(os.path.join(HERE, "recomputed_numbers_v04.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=str)
        fh.write("\n")
    print(json.dumps(out, indent=1, default=str)[:6000])


if __name__ == "__main__":
    sys.exit(main())
