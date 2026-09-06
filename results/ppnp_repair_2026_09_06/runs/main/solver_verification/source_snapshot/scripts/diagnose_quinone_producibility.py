"""Post-run structural/numerical checks for the frozen Putida quinone study.

No reactions are selected for a patch and no new organisms are evaluated. Ideal
quinone sources below are diagnostic boundary reactions, not biological models.
The prespecified runner, its inputs, and its result matrices are never modified.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
from pathlib import Path
import sys

import cobra
import highspy
import numpy as np
import scipy.sparse as sp
import swiglpk

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from gembench import wbm as W  # noqa: E402
from gembench.fitness_browser import base_medium, carbon_source_conditions, load_organism  # noqa: E402
from gembench.media import apply_medium  # noqa: E402
from gembench.protocols.carbon_fitness_generic import GenericParams, complete_medium_transport  # noqa: E402
from scripts.run_quinone_biomass_sensitivity import apply_biomass_intervention, prepare  # noqa: E402


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def safe(value):
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {k: safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [safe(v) for v in value]
    return value


def write_json(path, value):
    with path.open("x") as fh:
        fh.write(json.dumps(safe(value), indent=2, allow_nan=False) + "\n")


def matrix_model(model):
    """Translate the ordinary stoichiometric LP for an independent HiGHS solve."""
    S = sp.csc_matrix(cobra.util.array.create_stoichiometric_matrix(model, array_type="lil"))
    objective = cobra.util.solver.linear_reaction_coefficients(model)
    return W.WBM(name=model.id, rxns=np.array([r.id for r in model.reactions]),
                 mets=np.array([m.id for m in model.metabolites]), S=S,
                 b=np.zeros(len(model.metabolites)), csense=np.full(len(model.metabolites), "E"),
                 C=sp.csc_matrix((0, len(model.reactions))), d=np.array([]), dsense=np.array([]), ctrs=np.array([]),
                 lb=np.array([r.lower_bound for r in model.reactions]),
                 ub=np.array([r.upper_bound for r in model.reactions]),
                 c=np.array([objective.get(r, 0.0) for r in model.reactions]), osense="max")


def pool_certificate(model):
    """Exact row-sum certificate for the combined q8/q8h2 pool."""
    matrix = matrix_model(model)
    q8 = list(matrix.mets).index("q8_c")
    q8h2 = list(matrix.mets).index("q8h2_c")
    row = (matrix.S.getrow(q8) + matrix.S.getrow(q8h2)).toarray().ravel()
    nonzero = [{"reaction": str(matrix.rxns[i]), "coefficient": float(row[i])}
               for i in np.flatnonzero(row != 0)]
    incident = sorted({r.id for mid in ("q8_c", "q8h2_c")
                       for r in model.metabolites.get_by_id(mid).reactions})
    return {"metabolite_rows": ["q8_c", "q8h2_c"], "row_weights": [1, 1],
            "n_incident_reactions": len(incident), "incident_reactions": incident,
            "nonzero_summed_row_coefficients": nonzero,
            "max_abs_summed_row_coefficient": float(np.max(np.abs(row))),
            "exactly_conserved_pool": len(nonzero) == 0}


def solve_case(model, solver, tolerance):
    matrix = matrix_model(model)
    if solver == "glpk":
        model.solver = "glpk"
        model.solver.configuration.tolerances.feasibility = tolerance
        solution = model.optimize()
        status = model.solver.status
        values = solution.fluxes.reindex(matrix.rxns).to_numpy() if status == "optimal" else None
        objective = solution.objective_value if status == "optimal" else None
        effective_tolerance = model.solver.configuration.tolerances.feasibility
    else:
        answer = W.solve_highs(matrix, method="simplex", threads=0, feas_tol=tolerance,
                               opt_tol=tolerance, time_limit=60.0, extra_options={"parallel": "off"})
        status = answer.status
        values = answer.x if status == "Optimal" else None
        objective = answer.objective if status == "Optimal" else None
        effective_tolerance = tolerance
    record = {"solver": solver, "feasibility_tolerance": effective_tolerance,
              "status": status, "objective": objective}
    if values is not None:
        certificate = W.certify(matrix, values)
        record["primal_certificate"] = certificate
        fluxes = dict(zip(matrix.rxns, values))
        record["growth_flux"] = float(fluxes["Growth"])
        record["diagnostic_fluxes"] = {rid: float(value) for rid, value in fluxes.items() if rid.startswith("DIAG_")}
        residual = matrix.S @ values
        record["quinone_balance_residuals"] = {mid: float(residual[list(matrix.mets).index(mid)])
                                                for mid in ("q8_c", "q8h2_c")}
    return record


def add_boundary(model, metabolite, kind):
    rid = f"DIAG_{kind}_{metabolite}"
    reaction = cobra.Reaction(rid, lower_bound=0, upper_bound=1000)
    reaction.add_metabolites({model.metabolites.get_by_id(metabolite): 1 if kind == "SOURCE" else -1})
    model.add_reactions([reaction])
    return rid


def independent_score(sim, fit, wt, growth_threshold, fitness_threshold):
    grows = np.isfinite(wt) & (wt >= growth_threshold)
    good = np.isfinite(sim) & np.isfinite(fit) & grows[None, :]
    pred, observed = sim >= growth_threshold, fit >= fitness_threshold
    counts = {"tp": int((good & pred & observed).sum()), "tn": int((good & ~pred & ~observed).sum()),
              "fp": int((good & pred & ~observed).sum()), "fn": int((good & ~pred & observed).sum())}
    tp, tn, fp, fn = [counts[k] for k in ("tp", "tn", "fp", "fn")]
    denominator = (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)
    mcc = (tp * tn - fp * fn) / math.sqrt(denominator) if denominator else None
    return {"n_conditions_wt_grows": int(grows.sum()), "n_finite_pairs": int(good.sum()),
            "confusion": counts, "mcc": mcc}


def independent_summary(study):
    declared = json.loads((study / "summary.json").read_text())
    params = declared["study"]["protocol"]
    threshold, fitness_threshold = params["growth_threshold"], params["fitness_threshold"]
    base = dict(np.load(study / "no_quinone_demand/matrices.npz", allow_pickle=False))
    arms, changed = [], []
    for arm in declared["runs"]:
        name = arm["arm"]
        saved = dict(np.load(study / name / "matrices.npz", allow_pickle=False))
        for key in ("browser_genes", "conditions", "fitness"):
            assert np.array_equal(base[key], saved[key], equal_nan=True) if saved[key].dtype.kind == "f" else np.array_equal(base[key], saved[key])
        score = independent_score(saved["sim_growth"], saved["fitness"], saved["wt_growth"], threshold, fitness_threshold)
        assert score["n_conditions_wt_grows"] == arm["condition_level"]["n_conditions_wt_grows"]
        recorded_score = arm["gene_level_conditions_where_wt_grows"]
        if score["n_finite_pairs"]:
            assert score["confusion"] == recorded_score["confusion"]
            assert abs(score["mcc"] - recorded_score["mcc"]["point"]) < 1e-12
        common = (base["wt_growth"] >= threshold) & (saved["wt_growth"] >= threshold)
        common_wt = np.where(common, base["wt_growth"], 0.0)
        a = independent_score(base["sim_growth"], base["fitness"], common_wt, threshold, fitness_threshold)
        b = independent_score(saved["sim_growth"], saved["fitness"], common_wt, threshold, fitness_threshold)
        arms.append({"arm": name, **score, "conditions_both_grow": saved["conditions"][common].tolist(),
                     "matched_baseline_mcc": a["mcc"], "matched_arm_mcc": b["mcc"],
                     "max_wt_flux": float(np.max(saved["wt_growth"])),
                     "wt_rates_in_0.0005_to_0.002": int(((saved["wt_growth"] >= 0.0005) & (saved["wt_growth"] <= 0.002)).sum())})
        if name != "no_quinone_demand":
            if b["mcc"] is not None:
                comparison = next(c for c in declared["comparisons"] if c["arm"] == name)
                assert abs(comparison["paired_mcc_difference_B_minus_A"]["point"] - (b["mcc"] - a["mcc"])) < 1e-12
            difference = ((base["sim_growth"] >= threshold) != (saved["sim_growth"] >= threshold)) & common[None, :]
            for i, j in zip(*np.where(difference)):
                changed.append({"arm": name, "gene": str(base["browser_genes"][i]), "condition": str(base["conditions"][j]),
                                "before": float(base["sim_growth"][i, j]), "after": float(saved["sim_growth"][i, j]),
                                "fitness": float(base["fitness"][i, j])})
    assert changed == json.loads((study / "changed_predictions.json").read_text())
    assert len(changed) == declared["n_changed_gene_condition_predictions"]
    prior_paths = {
        "incomplete_attempt_baseline": ROOT / "results/quinone_biomass_2026_09_06_attempt01/no_quinone_demand/matrices.npz",
        "archived_cycle6": ROOT / "results/carbon_fitness_multi/Putida/Pseudomonas_putida_KT2440_xml_gapfilled__gapfilled__patched-v0.1+model-v0.3+gpr-v0.2+0.3+0.4+medium__nodroprich/matrices.npz",
        "development_sprint_cycle7": ROOT / "results/development_sprint_2026_09_06/runs/Putida_cycle7/artifacts/Putida/Pseudomonas_putida_KT2440_xml_gapfilled__gapfilled__patched-v0.1+model-v0.4+gpr-v0.2+0.3+0.4+medium__nodroprich/matrices.npz",
    }
    comparisons = []
    for name, path in prior_paths.items():
        if not path.exists():
            comparisons.append({"comparison": name, "available": False})
            continue
        old = dict(np.load(path, allow_pickle=False))
        same_labels = all(np.array_equal(base[k], old[k]) for k in ("browser_genes", "conditions"))
        same_sets = all(set(base[k]) == set(old[k]) for k in ("browser_genes", "conditions"))
        for data in (base, old):
            for key in ("browser_genes", "conditions"):
                assert len(set(data[key])) == len(data[key]), f"Duplicate {key} in comparison {name}"
        record = {"comparison": name, "available": True, "path": str(path.relative_to(ROOT)),
                  "sha256": sha256(path), "same_gene_condition_labels_in_order": same_labels,
                  "same_gene_condition_label_sets": same_sets}
        if same_sets:
            old_genes = {gene: i for i, gene in enumerate(old["browser_genes"])}
            old_conditions = {condition: i for i, condition in enumerate(old["conditions"])}
            gi = np.array([old_genes[gene] for gene in base["browser_genes"]])
            ci = np.array([old_conditions[condition] for condition in base["conditions"]])
            old_wt = old["wt_growth"][ci]
            old_ko = old["sim_growth"][np.ix_(gi, ci)]
            old_fit = old["fitness"][np.ix_(gi, ci)]
            finite_wt = np.isfinite(base["wt_growth"]) & np.isfinite(old_wt)
            finite_ko = np.isfinite(base["sim_growth"]) & np.isfinite(old_ko)
            record.update(axis_alignment="exact gene and condition identifiers",
                          same_wt_finite_mask=bool(np.array_equal(np.isfinite(base["wt_growth"]), np.isfinite(old_wt))),
                          same_ko_finite_mask=bool(np.array_equal(np.isfinite(base["sim_growth"]), np.isfinite(old_ko))),
                          n_jointly_finite_ko=int(finite_ko.sum()),
                          max_abs_wt_difference=float(np.max(np.abs(base["wt_growth"][finite_wt] - old_wt[finite_wt]))) if finite_wt.any() else None,
                          max_abs_ko_difference=float(np.max(np.abs(base["sim_growth"][finite_ko] - old_ko[finite_ko]))) if finite_ko.any() else None,
                          n_changed_binary_ko_predictions=int((finite_ko & ((base["sim_growth"] >= threshold) != (old_ko >= threshold))).sum()),
                          same_experimental_fitness=bool(np.array_equal(base["fitness"], old_fit, equal_nan=True)))
        comparisons.append(record)
    return {"status": "recomputed_from_saved_matrices", "method": "Independent confusion counts and closed-form MCC; no new phenotype simulation",
            "arms": arms, "changed_predictions": changed, "baseline_comparisons": comparisons}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", type=Path, default=ROOT / "results/quinone_biomass_2026_09_06")
    parser.add_argument("--out", type=Path, help="New diagnostic output directory; defaults to <study>/diagnostics")
    args = parser.parse_args()
    out = args.out or args.study / "diagnostics"
    if out.exists():
        raise FileExistsError(f"Diagnostic output already exists: {out}; choose a fresh --out directory")
    out.mkdir(parents=True)
    manifest_path = args.study / "prespecified_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    changed_inputs = [path for path, digest in manifest["input_sha256"].items() if sha256(ROOT / path) != digest]
    if changed_inputs:
        raise ValueError(f"Frozen main-study inputs changed: {changed_inputs}")
    params = GenericParams(**manifest["plan"]["protocol"])
    with gzip.open(ROOT / "models/gapfilled/Putida.xml.gz", "rt") as fh:
        template = cobra.io.read_sbml_model(fh)
    curated = cobra.io.read_sbml_model(str(ROOT / "models/bigg/iJN1463.xml"))
    template_coefficient = template.reactions.Growth.get_coefficient("mql8_c")
    curated_coefficient = curated.reactions.BIOMASS_KT2440_WT3.get_coefficient("q8h2_c")
    fb = load_organism("Putida")
    conditions = carbon_source_conditions(fb)
    model, gene_map, _, _ = prepare(fb)
    model.solver = "glpk"
    for exchange in model.exchanges:
        exchange.bounds = (0, 1000)
    added = complete_medium_transport(model, sorted({c.media for c in conditions if c.bigg_ids}), params.medium_completion_exclude)
    glucose = next(c for c in conditions if c.name == "D-Glucose" and c.media == "MOPS minimal media_noCarbon")
    apply_medium(model, base_medium(glucose.media), close_all=True)
    for exchange in glucose.exchanges:
        model.reactions.get_by_id(exchange).lower_bound = params.carbon_uptake
    certificate = pool_certificate(model)
    assert certificate["exactly_conserved_pool"], "The closed-pool interpretation must not be assumed if the row certificate fails"
    configs = {c["id"]: c for c in manifest["plan"]["configurations"]}
    cases = [("glucose_growth_without_quinone_demand", model.copy())]
    for name in ("ubiquinol_template_amount", "ubiquinol_curated_amount", "menaquinol_restored_control"):
        variant = model.copy()
        apply_biomass_intervention(variant, configs[name], template_coefficient, curated_coefficient)
        cases.append(("glucose_growth_" + name, variant))
        if name.startswith("ubiquinol"):
            rescue = variant.copy()
            add_boundary(rescue, "q8h2_c", "SOURCE")
            cases.append(("glucose_growth_" + name + "_ideal_source_rescue", rescue))
            forced = variant.copy()
            forced.reactions.Growth.lower_bound = 0.1
            cases.append(("glucose_growth_" + name + "_forced_0.1", forced))
    for mid in ("q8_c", "q8h2_c", "mql8_c"):
        demand = model.copy()
        demand.reactions.Growth.lower_bound = 0
        objective = add_boundary(demand, mid, "DEMAND")
        demand.objective = objective
        cases.append(("glucose_max_" + mid + "_demand_no_growth_requirement", demand))
    # These three genes are numerical checks of saved changed calls, not candidates for patch selection.
    changes = json.loads((args.study / "changed_predictions.json").read_text())
    for change in changes:
        config = configs[change["arm"]]
        condition = next(c for c in conditions if c.key == change["condition"])
        check = model.copy()
        apply_medium(check, base_medium(condition.media), close_all=True)
        for exchange in condition.exchanges:
            check.reactions.get_by_id(exchange).lower_bound = params.carbon_uptake
        apply_biomass_intervention(check, config, template_coefficient, curated_coefficient)
        gene = next(g for g, browser_gene in gene_map.model_to_browser.items() if browser_gene == change["gene"])
        check.genes.get_by_id(gene).knock_out()
        cases.append(("saved_changed_call_" + change["gene"], check))
    results = []
    for name, case in cases:
        record = {"case": name, "quinone_pool_certificate": pool_certificate(case), "solves": []}
        for solver, tol in (("glpk", 1e-7), ("glpk", 1e-9), ("highs", 1e-9)):
            answer = solve_case(case.copy(), solver, tol)
            record["solves"].append(answer)
            print(name, solver, tol, answer["status"], answer["objective"], flush=True)
        results.append(record)
    report = {"created_at": datetime.now(timezone.utc).isoformat(),
              "role": "Post-run structural and numerical diagnostic; not prespecified model selection or external validation",
              "script_sha256": sha256(__file__), "main_manifest_sha256": sha256(manifest_path),
              "source_metadata_omission": {"path": "models/gapfilled/Putida_gapfill.json", "sha256_now": sha256(ROOT / "models/gapfilled/Putida_gapfill.json"),
                  "note": "Read by load_model for source text but omitted from the main manifest; hash captured only during this post-run diagnostic. Numerical model bytes were frozen."},
              "versions": {"cobra": cobra.__version__, "glpk": swiglpk.glp_version(), "highs": highspy.Highs().version()},
              "condition": glucose.key, "carbon_uptake": params.carbon_uptake, "medium_completion_added": added,
              "coefficients": {"template_menaquinol": template_coefficient, "curated_ubiquinol": curated_coefficient},
              "baseline_quinone_pool_certificate": certificate,
              "structural_conclusion": "The q8_c and q8h2_c balance rows sum to zero in the prepared baseline. Adding only -epsilon q8h2_c to Growth makes their sum -epsilon*v_Growth=0, forcing zero growth for any nonzero epsilon. This is a model-structure property, not a statement that P. putida lacks quinone synthesis.",
              "source_rescue_caveat": "An ideal source bypasses missing pool synthesis for diagnosis and is not an accepted biological reaction or patch.",
              "cases": results}
    write_json(out / "producibility.json", report)
    write_json(out / "independent_summary.json", independent_summary(args.study))


if __name__ == "__main__":
    main()
