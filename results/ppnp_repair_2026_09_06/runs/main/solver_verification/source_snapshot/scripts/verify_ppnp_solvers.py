"""Cross-check seven saved PpnP physical arms with GLPK and HiGHS.

No fitness values are parsed. Reconstruct every saved medium and gene identity
before solving. Reuse the guarded ordinary-LP translator and solver checks from
the earlier independent quinone audit; never change frozen experiment inputs.
Selected fluxes need not match another optimum. Their bounds and fully covered
balance rows are checked, while their reported whole-model residual cannot be
reconstructed from a partial flux vector. Preserve source snapshots and failures
in a fresh verification directory.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from fractions import Fraction
from importlib.metadata import version
import json
import math
from pathlib import Path
import platform
import shutil
import sys

import cobra
import swiglpk

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from gembench.study import verify_study  # noqa: E402
from scripts.verify_quinone_repair_solvers import (  # noqa: E402
    FEASIBILITY_TOLERANCE, OBJECTIVE_ABS_TOLERANCE, OBJECTIVE_REL_TOLERANCE,
    RESIDUAL_TOLERANCE, assert_ordinary_lp, check_case, close_objective,
    file_record, metadata_conditions, read_json, reconstruct, resolve_targets,
    unchanged, validate_arm, write_json,
)

ARM_IDS = ('parent', 'adenosine_forward', 'inosine_forward', 'both_forward',
           'parent_rhcys_closed', 'both_forward_rhcys_closed', 'both_reversible')
SINGLE_LOCI = ('PP_2458', 'PP_4248', 'PP_1777', 'PP_4976', 'PP_0591', 'PP_2460', 'PP_3254', 'PP_5317')
MECHANISM_IDS = ('wild_type', 'rbk', 'rbk_ppnp', 'rbk_ppm', 'rbk_ahcy', 'rbk_ada', 'rbk_insh', 'rbk_added_closed')
ADDED_IDS = ('PUNP1', 'PUNP5')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def finite_number(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def unique_aliases(model, mapping, loci):
    require(mapping['org_id'] == 'Putida', 'Wrong gene-map organism')
    result = {}
    for locus in loci:
        hits = sorted(gid for gid, name in mapping['model_to_browser'].items() if name == locus)
        require(len(hits) <= 1, f'Ambiguous locus identity: {locus}')
        gid = hits[0] if hits else locus
        require(not hits or gid in model.genes, f'Mapped gene absent from saved model: {locus}/{gid}')
        require(gid not in model.genes or mapping['model_to_browser'].get(gid) == locus,
                f'Gene identity collision: {locus}/{gid}')
        result[locus] = gid
    require(len(set(result.values())) == len(result), 'Distinct loci resolve to the same model gene')
    return result


def mechanism_inventory(model, mapping, saved, recipe):
    expected = recipe['mechanistic_cases']
    require([row['id'] for row in saved['cases']] == [row['id'] for row in expected],
            'Saved mechanistic case inventory differs from recipe')
    aliases = unique_aliases(model, mapping, recipe['targeted_gene_loci'])
    validated = []
    for row, specification in zip(saved['cases'], expected):
        require(row['loci'] == specification['loci'], f'Mechanism loci differ: {row["id"]}')
        require(row['close_added_reactions'] == specification.get('close_added_reactions', False),
                f'Mechanism closure differs: {row["id"]}')
        missing = [locus for locus in specification['loci'] if aliases[locus] not in model.genes or
                   not model.genes.get_by_id(aliases[locus]).reactions]
        if missing:
            require(row['status'] == 'unrepresented_gene_or_function' and row['growth'] is None and
                    row.get('unrepresented') == missing, f'Invalid absent/orphan case: {row["id"]}')
            require('quinone_production' not in row and 'selected_fluxes' not in row,
                    f'Unrepresented case fabricated numerical predictions: {row["id"]}')
        else:
            require(row['status'] == 'optimal' and finite_number(row['growth']),
                    f'Invalid saved mechanistic status/value: {row["id"]}')
            require(('quinone_production' in row) == specification.get('probe_quinone', False),
                    f'Mechanistic production inventory differs: {row["id"]}')
            if 'quinone_production' in row:
                validate_probe_inventory(row['quinone_production'], recipe)
        validated.append((specification, row, missing))
    return aliases, validated


def validate_probe_inventory(probes, recipe):
    require([row['metabolite'] for row in probes] == list(recipe['pool_weights']), 'Production probe metabolite inventory differs')
    settings = recipe['physical_settings']
    for row in probes:
        rate = row['maximum_demand_flux']
        require(row['status'] == 'optimal' and finite_number(rate), 'Failed/nonfinite saved production probe')
        require(row['threshold'] == settings['production_threshold'] and row['probe_capacity'] == settings['production_capacity'],
                'Saved production settings differ')
        require(row['producible_at_threshold'] == (rate >= row['threshold']), 'Saved production threshold classification differs')
        require(-RESIDUAL_TOLERANCE <= rate <= row['probe_capacity'] + RESIDUAL_TOLERANCE,
                'Saved production rate exceeds its bounds')


def apply_mechanism(model, aliases, specification):
    """Called only inside a model context so knockout/closure bounds restore."""
    for locus in specification['loci']:
        model.genes.get_by_id(aliases[locus]).knock_out()
    if specification.get('close_added_reactions', False):
        for rid in ADDED_IDS:
            if rid in model.reactions:
                model.reactions.get_by_id(rid).bounds = (0, 0)


def selected_flux_checks(model, row, recipe):
    """Check what a saved partial witness actually permits us to verify."""
    selected = row['selected_fluxes']
    expected = {rid for rid in recipe['witness_reactions'] if rid in model.reactions}
    require(set(selected) == expected, 'Selected witness reaction inventory differs')
    require(all(finite_number(value) for value in selected.values()), 'Nonfinite selected flux')
    require(close_objective(selected['Growth'], row['growth']), 'Selected Growth flux differs from saved objective')
    tolerance = recipe['physical_settings']['residual_tolerance']
    require(tolerance == RESIDUAL_TOLERANCE, 'Unexpected residual tolerance')
    for field in ('mass_balance_max_abs', 'bound_violation_max'):
        require(finite_number(row[field]) and 0 <= row[field] <= tolerance, f'Invalid saved witness residual claim: {field}')
    bound_error = max((max(model.reactions.get_by_id(rid).lower_bound - value,
                           value - model.reactions.get_by_id(rid).upper_bound, 0)
                       for rid, value in selected.items()), default=0.0)
    require(bound_error <= tolerance, 'Selected witness violates reaction/knockout bounds')
    rows = {}
    for met in model.metabolites:
        if met.reactions and {r.id for r in met.reactions} <= set(selected):
            rows[met.id] = float(sum(r.metabolites[met] * selected[r.id] for r in met.reactions))
    require(all(abs(value) <= tolerance for value in rows.values()), 'Selected witness violates a completely covered metabolite balance')
    pool_coefficients = {}
    for reaction in model.reactions:
        total = sum((Fraction.from_float(float(recipe['pool_weights'].get(met.id, 0))) * Fraction.from_float(float(value))
                     for met, value in reaction.metabolites.items()), Fraction(0))
        if total:
            pool_coefficients[reaction.id] = float(total)
    require(set(pool_coefficients) <= set(selected), 'Selected witness cannot reconstruct declared quinone-pool balance')
    pool_residual = sum(coefficient * selected[rid] for rid, coefficient in pool_coefficients.items())
    require(abs(pool_residual) <= tolerance, 'Selected witness violates aggregate quinone-pool balance')
    return {'selected_bound_violation_max': bound_error,
            'complete_metabolite_balance_residuals': rows,
            'aggregate_quinone_pool_residual': pool_residual,
            'reported_full_mass_balance_max_abs': row['mass_balance_max_abs'],
            'reported_full_bound_violation_max': row['bound_violation_max'],
            'scope': 'Selected bounds and completely covered/aggregate rows independently checked. The saved whole-model residual is only a reported claim because the full saved flux vector is unavailable. Independent solves below produce their own full primal certificates.'}


def crosscheck_duplicate_cases(physical, mechanism, recipe):
    by_id = {row['id']: row for row in mechanism['cases']}
    require(close_objective(by_id['wild_type']['growth'], physical['wild_type']['growth']), 'Physical/mechanistic wild type differs')
    rbk_single = next(row for row in physical['targeted_gene_predictions'] if row['locus'] == 'PP_2458')
    require(close_objective(by_id['rbk']['growth'], rbk_single['growth']), 'Physical/mechanistic RBK deletion differs')
    for row in mechanism['cases']:
        if row['status'] == 'optimal':
            require(row['growth'] >= -RESIDUAL_TOLERANCE and row['growth'] <= physical['wild_type']['growth'] + OBJECTIVE_ABS_TOLERANCE,
                    'Mechanistic restriction grows above the WT feasible optimum')


def check_mechanisms(model, mapping, mechanism, recipe, emit):
    aliases, cases = mechanism_inventory(model, mapping, mechanism, recipe)
    for specification, row, missing in cases:
        name = 'mechanism:' + specification['id']
        if missing:
            emit(name, {'skipped': True, 'valid': True, 'reason': 'unrepresented_gene_or_function', 'loci': missing})
            continue
        with model:
            apply_mechanism(model, aliases, specification)
            partial_check = selected_flux_checks(model, row, recipe)
            result = check_case(model, row['growth'], recipe['protocol']['growth_threshold'])
            # check_case's GLPK primal remains available after the independent
            # HiGHS solve. Capture it without an additional optimization.
            independent_selected = None
            if result['solves'][0]['status'].lower() == 'optimal':
                independent_selected = {rid: float(model.reactions.get_by_id(rid).flux)
                                        for rid in recipe['witness_reactions'] if rid in model.reactions}
                require(all(math.isfinite(v) for v in independent_selected.values()), 'Nonfinite independent selected flux')
            result.update(loci=specification['loci'], close_added_reactions=specification.get('close_added_reactions', False),
                          saved_partial_witness_check=partial_check,
                          independent_glpk_selected_fluxes=independent_selected,
                          witness_comparison='Objectives and feasibility are compared; individual selected fluxes are not required to agree across distinct optima.')
            emit(name, result)
            for probe in row.get('quinone_production', []):
                with model:
                    rid = f'DIAG_DEMAND_{probe["metabolite"]}'
                    require(rid not in model.reactions, 'Existing diagnostic reaction would be overwritten')
                    demand = cobra.Reaction(rid, lower_bound=0, upper_bound=probe['probe_capacity'])
                    demand.add_metabolites({model.metabolites.get_by_id(probe['metabolite']): -1})
                    model.add_reactions([demand])
                    model.objective = demand
                    # Mechanistic probes retain current growth bounds, exactly
                    # as the frozen runner; the parent lower bound is zero.
                    emit(name + ':production:' + probe['metabolite'],
                         check_case(model, probe['maximum_demand_flux'], probe['threshold']))


def summarize(cases, n_arms, scope):
    solves = [solve for case in cases for solve in case.get('solves', [])]
    def certificate_max(field):
        values = [solve.get('primal_certificate', {}).get(field) for solve in solves]
        return max((value for value in values if finite_number(value)), default=None)
    return {'valid': all(case['valid'] for case in cases), 'n_arms': n_arms, 'n_cases': len(cases),
            'n_skipped': sum(bool(case.get('skipped')) for case in cases), 'n_solver_runs': len(solves),
            'failed_cases': [case for case in cases if not case['valid']],
            'max_solver_objective_disagreement': max((abs(case['solves'][0]['objective'] - case['solves'][1]['objective'])
                                                     for case in cases if len(case.get('solves', [])) == 2 and
                                                     all(finite_number(s['objective']) for s in case['solves'])), default=None),
            'max_independent_mass_balance_residual': certificate_max('max_abs_S_residual'),
            'max_independent_bound_violation': certificate_max('max_bound_violation'),
            'n_solver_runs_without_full_finite_certificate': sum(any(not finite_number(solve.get('primal_certificate', {}).get(field))
                                                                    for field in ('max_abs_S_residual', 'max_bound_violation'))
                                                                for solve in solves),
            'limitations': scope}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True, help='Fresh verification directory; never overwritten')
    args = parser.parse_args(argv)
    study, out = args.study.resolve(), args.out.resolve()
    required = [study / 'manifest.json']
    required += [study / arm / name for arm in ARM_IDS
                 for name in ('model.xml.gz', 'physical.json', 'gene_map.json', 'intervention.json', 'mechanism.json')]
    require(all(path.is_file() for path in required), 'All seven model/physical/mechanistic arms must exist before numerical verification')
    out.mkdir(parents=True, exist_ok=False)
    cases = []
    try:
        manifest = read_json(study / 'manifest.json')
        recipe = manifest['plan']['recipe']
        require(tuple(row['id'] for row in recipe['configurations']) == ARM_IDS, 'Unexpected PpnP arm inventory')
        require(tuple(recipe['targeted_gene_loci']) == SINGLE_LOCI, 'Unexpected single-gene inventory')
        require(tuple(row['id'] for row in recipe['mechanistic_cases']) == MECHANISM_IDS, 'Unexpected mechanism inventory')
        require(recipe['organism'] == 'Putida' and recipe['evaluation_role'] == 'development', 'Wrong study organism/role')
        require(set(recipe['pool_weights']) == {'q8_c', 'q8h2_c'}, 'Wrong quinone probe inventory')
        verification = verify_study(ROOT, manifest)
        write_json(out / 'main_input_verification_before.json', verification)
        require(verification['valid'], 'Frozen experiment inputs changed')
        sources = {Path(__file__).resolve()}
        for module in tuple(sys.modules.values()):
            filename = getattr(module, '__file__', None)
            if filename:
                path = Path(filename).resolve()
                if path.is_relative_to(ROOT) and path.suffix == '.py' and '.venv' not in path.parts:
                    sources.add(path)
        inputs = set(required) | sources | {ROOT / row['path'] for row in manifest['files']}
        records = [file_record(path) for path in sorted(inputs)]
        for source in sorted(sources):
            destination = out / 'source_snapshot' / source.relative_to(ROOT)
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, destination)
        scope = ('Independent GLPK/HiGHS agreement and primal feasibility for seven saved physical arms, represented single-gene deletions, '
                 'mechanistic restrictions and production probes. No fitness values parsed or energy probes re-solved; no exact optimality, '
                 'enzyme validation or independent biological performance claim. Saved partial flux witnesses do not permit reconstruction '
                 'of their whole-model residual; selected fluxes can differ between equally optimal solutions.')
        provenance = {'created_at': datetime.now(timezone.utc).isoformat(), 'inputs': records,
                      'versions': {name: version(name) for name in ('cobra', 'highspy', 'numpy', 'scipy', 'optlang', 'swiglpk', 'pandas')},
                      'python': sys.version, 'platform': platform.platform(), 'glpk_library': swiglpk.glp_version(),
                      'feasibility_tolerance': FEASIBILITY_TOLERANCE, 'residual_tolerance': RESIDUAL_TOLERANCE,
                      'objective_absolute_tolerance': OBJECTIVE_ABS_TOLERANCE,
                      'objective_relative_tolerance': OBJECTIVE_REL_TOLERANCE, 'scope': scope}
        write_json(out / 'provenance.json', provenance)
        conditions = metadata_conditions(ROOT)
        prepared = []
        # Complete every context, alias and partial-witness preflight before
        # optimizing any model; mismatch never becomes a no-growth result.
        for arm, configuration in zip(ARM_IDS, recipe['configurations']):
            physical = read_json(study / arm / 'physical.json')
            mechanism = read_json(study / arm / 'mechanism.json')
            mapping = read_json(study / arm / 'gene_map.json')
            intervention = read_json(study / arm / 'intervention.json')
            require(intervention['configuration'] == configuration, 'Intervention differs from recipe')
            model = cobra.io.read_sbml_model(str(study / arm / 'model.xml.gz'))
            context = reconstruct(model, physical, recipe, conditions)
            aliases = unique_aliases(model, mapping, recipe['targeted_gene_loci'])
            require(intervention['gene_aliases'] == {'PP_4248': aliases['PP_4248']}, 'Intervention PpnP alias differs')
            resolve_targets(model, mapping, physical, recipe)
            validate_probe_inventory(physical['production_at_growth_lower_bound_zero'], recipe)
            _, planned = mechanism_inventory(model, mapping, mechanism, recipe)
            crosscheck_duplicate_cases(physical, mechanism, recipe)
            selected_checks = {}
            for specification, row, missing in planned:
                if not missing:
                    with model:
                        apply_mechanism(model, aliases, specification)
                        assert_ordinary_lp(model)
                        selected_checks[specification['id']] = selected_flux_checks(model, row, recipe)
            write_json(out / f'{arm}_context.json', {'reconstructed_context': context, 'unique_target_aliases': aliases,
                                                   'saved_partial_witness_checks': selected_checks})
            prepared.append((arm, model, physical, mapping, mechanism))
        unchanged(records)
        with (out / 'cases.jsonl').open('x') as stream:
            for arm, model, physical, mapping, mechanism in prepared:
                def emit(name, result):
                    record = {'arm': arm, 'case': name, **result}
                    cases.append(record)
                    stream.write(json.dumps(record, allow_nan=False) + '\n')
                    stream.flush()
                validate_arm(model, physical, mapping, recipe, emit)
                check_mechanisms(model, mapping, mechanism, recipe, emit)
                print(f'{arm}: physical and mechanistic numerical verification complete', flush=True)
        unchanged(records)
        verification = verify_study(ROOT, manifest)
        write_json(out / 'main_input_verification_after.json', verification)
        summary = summarize(cases, len(prepared), scope)
        summary['valid'] = summary['valid'] and verification['valid']
        write_json(out / 'summary.json', summary)
        print(json.dumps(summary, indent=2), flush=True)
        return 0 if summary['valid'] else 1
    except Exception as error:
        write_json(out / 'failure.json', {'type': type(error).__name__, 'message': str(error), 'completed_cases': len(cases)})
        raise


if __name__ == '__main__':
    raise SystemExit(main())
