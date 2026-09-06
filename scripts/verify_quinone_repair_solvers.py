"""Independently check saved quinone-repair LPs, without loading fitness values.

This checks numerical agreement and primal residuals, not enzyme chemistry,
experimental essentiality, evaluation independence, or exact optimality proofs.
The main experiment's code and outputs are never modified. A fresh directory
preserves source snapshots, inputs, completed cases and any failure.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
from importlib.metadata import version
import json
import math
from pathlib import Path
import platform
import shutil
import sys
from types import SimpleNamespace

import cobra
import pandas as pd
import swiglpk

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from gembench.fitness_browser import base_medium, carbon_source_conditions  # noqa: E402
from gembench.media import apply_medium  # noqa: E402
from gembench.protocols import carbon_fitness_generic as P  # noqa: E402
from gembench.study import verify_study  # noqa: E402
from scripts.diagnose_quinone_producibility import matrix_model, solve_case  # noqa: E402

ARM_IDS = ('baseline', 'demand_only', 'curated_path_only', 'curated_path_template',
           'curated_path_curated', 'sequence_path_template')
FEASIBILITY_TOLERANCE = 1e-9
RESIDUAL_TOLERANCE = 1e-8
OBJECTIVE_ABS_TOLERANCE = 1e-8
OBJECTIVE_REL_TOLERANCE = 1e-8


def read_json(path):
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError(f'Duplicate JSON key: {key}')
            result[key] = value
        return result
    def invalid(value):
        raise ValueError(f'Nonfinite JSON value: {value}')
    return json.loads(Path(path).read_text(), object_pairs_hook=pairs, parse_constant=invalid)


def write_json(path, value):
    with Path(path).open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')


def file_record(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return {'path': str(path), 'size_bytes': path.stat().st_size, 'sha256': digest.hexdigest()}


def unchanged(records):
    changed = [r['path'] for r in records if file_record(Path(r['path'])) != r]
    if changed:
        raise ValueError(f'Checker inputs changed during the run: {changed}')


def coefficients(expression):
    return {term: float(value) for term, value in expression.as_coefficients_dict().items()
            if value != 0}


def assert_ordinary_lp(model):
    """Reject every solver feature omitted by the old matrix translator.

    SBML generally does not retain arbitrary solver constraints. This guard
    certifies the loaded/reconstructed LP, not lost in-memory annotations.
    """
    if model.objective_direction != 'max':
        raise ValueError('The matrix helper supports maximization only')
    expected_variables = {v.name for r in model.reactions
                          for v in (r.forward_variable, r.reverse_variable)}
    if {v.name for v in model.variables} != expected_variables:
        raise ValueError('Unsupported custom solver variables')
    if {c.name for c in model.constraints} != {m.id for m in model.metabolites}:
        raise ValueError('Unsupported custom solver constraints')
    for reaction in model.reactions:
        lb, ub = reaction.bounds
        if not all(math.isfinite(v) for v in (lb, ub)) or lb > ub:
            raise ValueError(f'Unsupported reaction bounds: {reaction.id}')
        if lb >= 0:
            expected = ((lb, ub), (0, 0))
        elif ub <= 0:
            expected = ((0, 0), (-ub, -lb))
        else:
            expected = ((0, ub), (0, -lb))
        for variable, bounds in zip((reaction.forward_variable, reaction.reverse_variable), expected):
            if variable.type != 'continuous' or (variable.lb, variable.ub) != bounds:
                raise ValueError(f'Unsupported split-variable bounds/type: {variable.name}')
    for metabolite in model.metabolites:
        row = metabolite.constraint
        expected = {}
        for reaction in metabolite.reactions:
            value = float(reaction.metabolites[metabolite])
            if not math.isfinite(value):
                raise ValueError(f'Nonfinite stoichiometry: {reaction.id}/{metabolite.id}')
            if value:
                expected[reaction.forward_variable] = value
                expected[reaction.reverse_variable] = -value
        if row.lb != 0 or row.ub != 0 or coefficients(row.expression) != expected:
            raise ValueError(f'Unsupported modified metabolite balance: {metabolite.id}')
    expected = {}
    for reaction, value in cobra.util.solver.linear_reaction_coefficients(model).items():
        if not math.isfinite(value):
            raise ValueError('Nonfinite objective coefficient')
        if value:
            expected[reaction.forward_variable] = value
            expected[reaction.reverse_variable] = -value
    if coefficients(model.objective.expression) != expected:
        raise ValueError('Unsupported non-reaction or asymmetric objective')


def metadata_conditions(root):
    """Only experiment metadata and the fitness column header are parsed."""
    folder = root / 'data/fitness_browser/Putida'
    experiments = pd.read_table(folder / 'experiments.tsv', dtype=str, keep_default_na=False)
    columns = pd.read_table(folder / 'fit_logratios.tsv', nrows=0).columns
    metadata = {'orgId', 'locusId', 'sysName', 'geneName', 'desc'}
    org = SimpleNamespace(experiments=experiments, fitness=SimpleNamespace(
        columns=[c.split(' ')[0] for c in columns if c not in metadata]))
    return carbon_source_conditions(org)


def assert_physical_context(model, saved, recipe, condition, added, missing):
    actual = {
        'condition': condition.key,
        'settings': recipe['physical_settings'],
        'objective': {'expression': str(model.objective.expression), 'direction': model.objective_direction},
        'growth_bounds': list(model.reactions.Growth.bounds),
        'atp_maintenance_bounds': list(model.reactions.ATPM.bounds) if 'ATPM' in model.reactions else None,
        'mandatory_flux_bounds': {r.id: list(r.bounds) for r in model.reactions
                                  if r.lower_bound > 0 or r.upper_bound < 0},
        'final_exchange_bounds': {r.id: list(r.bounds) for r in model.exchanges},
        'declared_carbon_source_ids': condition.bigg_ids,
        'declared_carbon_exchanges': condition.exchanges,
        'declared_carbon_uptake': P.GenericParams(**recipe['protocol']).carbon_uptake,
        'medium_completion_added': added,
        'missing_medium_components': missing,
    }
    differences = [key for key, value in actual.items() if saved.get(key) != value]
    if differences:
        raise ValueError(f'Reconstructed physical context differs: {differences}')
    objective = {r.id: c for r, c in cobra.util.solver.linear_reaction_coefficients(model).items()}
    if objective != {'Growth': 1.0}:
        raise ValueError('WT/KO stored growth requires the unit Growth objective')
    if saved['effective_feasibility_tolerance'] != FEASIBILITY_TOLERANCE:
        raise ValueError('Unexpected saved physical feasibility tolerance')
    assert_ordinary_lp(model)
    return actual


def reconstruct(model, saved, recipe, conditions):
    assert_ordinary_lp(model)  # reject features before any solver replacement
    params = P.GenericParams(**recipe['protocol'])
    model.solver = 'glpk'
    model.solver.configuration.tolerances.feasibility = FEASIBILITY_TOLERANCE
    for reaction in model.exchanges:
        reaction.bounds = (0, 1000)
    for gid in params.knockout_genes:
        model.genes.get_by_id(gid).knock_out()
    media = sorted({c.media for c in conditions if c.bigg_ids})
    added = P.complete_medium_transport(model, media, params.medium_completion_exclude) \
        if params.complete_medium_transport else []
    selected = recipe['physical_condition']
    matches = [c for c in conditions if c.name == selected['name'] and c.media == selected['media']]
    if len(matches) != 1:
        raise ValueError('Declared physical condition did not resolve uniquely')
    condition = matches[0]
    missing = apply_medium(model, base_medium(condition.media), close_all=True)
    for rid in condition.exchanges:
        model.reactions.get_by_id(rid).lower_bound = params.carbon_uptake
    context = assert_physical_context(model, saved, recipe, condition, added, missing)
    context['completion_media_names'] = media
    return context


def close_objective(left, right):
    return (isinstance(left, (int, float)) and isinstance(right, (int, float))
            and math.isfinite(left) and math.isfinite(right)
            and math.isclose(left, right, abs_tol=OBJECTIVE_ABS_TOLERANCE,
                             rel_tol=OBJECTIVE_REL_TOLERANCE))


def check_case(model, expected, threshold):
    assert_ordinary_lp(model)
    matrix = matrix_model(model)
    record = {'expected_objective': expected, 'classification_threshold': threshold,
              'n_reactions': len(matrix.rxns), 'n_metabolites': len(matrix.mets), 'solves': []}
    problems = []
    for solver in ('glpk', 'highs'):
        result = solve_case(model, solver, FEASIBILITY_TOLERANCE)
        record['solves'].append(result)
        if result['status'].lower() != 'optimal':
            problems.append(f'{solver} did not return optimal')
            continue
        certificate = result.get('primal_certificate', {})
        for key in ('max_abs_S_residual', 'max_bound_violation'):
            value = certificate.get(key)
            if value is None or not math.isfinite(value) or value > RESIDUAL_TOLERANCE:
                problems.append(f'{solver} failed residual certificate: {key}={value}')
        objective = result['objective']
        if not close_objective(objective, expected):
            problems.append(f'{solver} objective disagrees with saved result')
        if math.isfinite(objective) and (objective >= threshold) != (expected >= threshold):
            problems.append(f'{solver} threshold classification disagrees with saved result')
    if not close_objective(*(r['objective'] for r in record['solves'])):
        problems.append('GLPK and HiGHS objectives disagree')
    record.update(valid=not problems, problems=problems)
    return record


def resolve_targets(model, mapping, saved, recipe):
    rows = saved['targeted_gene_predictions']
    if [r['locus'] for r in rows] != recipe['targeted_gene_loci']:
        raise ValueError('Targeted gene list/order differs from recipe')
    resolved = []
    for row in rows:
        hits = [gid for gid, locus in mapping['model_to_browser'].items() if locus == row['locus']]
        if len(hits) > 1:
            raise ValueError(f'Ambiguous gene mapping: {row["locus"]}')
        gid = hits[0] if hits else row['locus']
        exists = gid in model.genes
        reactions = sorted(r.id for r in model.genes.get_by_id(gid).reactions) if exists else []
        if row['model_gene_id'] != (gid if exists else None) or row['associated_reactions'] != reactions:
            raise ValueError(f'Targeted gene representation differs: {row["locus"]}')
        skip = 'unrepresented_gene' if not exists else 'unrepresented_function' if not reactions else None
        if skip:
            if row['status'] != skip or row['growth'] is not None:
                raise ValueError(f'Invalid skipped gene result: {row["locus"]}')
        elif row['status'] != 'optimal' or not isinstance(row['growth'], (int, float)) or not math.isfinite(row['growth']):
            raise ValueError(f'Invalid saved targeted gene solve: {row["locus"]}')
        resolved.append((row, gid, skip))
    return resolved


def validate_arm(model, saved, mapping, recipe, emit):
    threshold = P.GenericParams(**recipe['protocol']).growth_threshold
    if saved['wild_type']['status'] != 'optimal':
        raise ValueError('Saved WT growth is not optimal')
    targets = resolve_targets(model, mapping, saved, recipe)
    probes = saved['production_at_growth_lower_bound_zero']
    if [p['metabolite'] for p in probes] != list(recipe['pool_weights']):
        raise ValueError('Saved production probes differ from declared pool')
    emit('wild_type', check_case(model, saved['wild_type']['growth'], threshold))
    for probe in probes:
        settings = recipe['physical_settings']
        if (probe['status'] != 'optimal' or probe['probe_capacity'] != settings['production_capacity']
                or probe['threshold'] != settings['production_threshold']):
            raise ValueError('Saved production probe settings differ from recipe')
        with model:
            model.reactions.Growth.lower_bound = 0
            rid = f'DIAG_DEMAND_{probe["metabolite"]}'
            if rid in model.reactions:
                raise ValueError(f'Existing diagnostic demand: {rid}')
            demand = cobra.Reaction(rid, lower_bound=0, upper_bound=probe['probe_capacity'])
            demand.add_metabolites({model.metabolites.get_by_id(probe['metabolite']): -1})
            model.add_reactions([demand])
            model.objective = demand
            emit(f'production:{probe["metabolite"]}',
                 check_case(model, probe['maximum_demand_flux'], probe['threshold']))
    for row, gid, skip in targets:
        if skip:
            emit(f'gene:{row["locus"]}', {'skipped': True, 'reason': skip, 'valid': True})
            continue
        with model:
            model.genes.get_by_id(gid).knock_out()
            record = check_case(model, row['growth'], threshold)
            record.update(model_gene_id=gid, associated_reactions=row['associated_reactions'],
                          interpretation='Conditional model prediction; no experimental validation',
                          wild_type_grows_at_threshold=saved['wild_type']['growth'] >= threshold)
            emit(f'gene:{row["locus"]}', record)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True, help='Fresh output directory; never overwritten')
    args = parser.parse_args(argv)
    study, out = args.study.resolve(), args.out.resolve()
    required = [study / 'manifest.json', study / 'physical_summary.json']
    required += [study / arm / name for arm in ARM_IDS for name in ('model.xml.gz', 'physical.json', 'gene_map.json')]
    if any(not p.is_file() for p in required):
        raise ValueError('All six saved physical arms must exist before numerical validation')
    out.mkdir(parents=True, exist_ok=False)
    cases = []
    try:
        manifest = read_json(study / 'manifest.json')
        recipe = manifest['plan']['recipe']
        verification = verify_study(ROOT, manifest)
        write_json(out / 'main_input_verification_before.json', verification)
        if not verification['valid']:
            raise ValueError('Original experiment inputs have changed')
        if tuple(a['id'] for a in recipe['configurations']) != ARM_IDS or set(recipe['pool_weights']) != {'q8_c', 'q8h2_c'}:
            raise ValueError('Checker is restricted to the declared six quinone-repair arms')
        sources = {Path(__file__).resolve()}
        for module in tuple(sys.modules.values()):
            filename = getattr(module, '__file__', None)
            if filename:
                path = Path(filename).resolve()
                if path.is_relative_to(ROOT) and path.suffix == '.py' and '.venv' not in path.parts:
                    sources.add(path)
        inputs = set(required) | sources | {ROOT / r['path'] for r in manifest['files']}
        records = [file_record(p) for p in sorted(inputs)]
        for path in sorted(sources):
            target = out / 'source_snapshot' / path.relative_to(ROOT)
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(path, target)
        provenance = {'created_at': datetime.now(timezone.utc).isoformat(), 'inputs': records,
                      'versions': {name: version(name) for name in ('cobra', 'highspy', 'numpy', 'scipy', 'optlang', 'swiglpk', 'pandas')},
                      'python': sys.version, 'platform': platform.platform(), 'glpk_library': swiglpk.glp_version(),
                      'feasibility_tolerance': FEASIBILITY_TOLERANCE, 'residual_tolerance': RESIDUAL_TOLERANCE,
                      'objective_absolute_tolerance': OBJECTIVE_ABS_TOLERANCE,
                      'objective_relative_tolerance': OBJECTIVE_REL_TOLERANCE,
                      'scope': 'Independent solver agreement and primal feasibility for saved ordinary LPs. No fitness values parsed, enzyme validation, or exact optimality proof.'}
        write_json(out / 'provenance.json', provenance)
        conditions = metadata_conditions(ROOT)
        prepared = []
        # Check every reconstruction and gene representation before optimizing any arm.
        for arm in ARM_IDS:
            saved = read_json(study / arm / 'physical.json')
            mapping = read_json(study / arm / 'gene_map.json')
            model = cobra.io.read_sbml_model(str(study / arm / 'model.xml.gz'))
            context = reconstruct(model, saved, recipe, conditions)
            resolve_targets(model, mapping, saved, recipe)
            write_json(out / f'{arm}_context.json', context)
            prepared.append((arm, model, saved, mapping))
        unchanged(records)
        with (out / 'cases.jsonl').open('x') as stream:
            for arm, model, saved, mapping in prepared:
                def emit(name, result):
                    record = {'arm': arm, 'case': name, **result}
                    cases.append(record)
                    stream.write(json.dumps(record, allow_nan=False) + '\n')
                    stream.flush()
                validate_arm(model, saved, mapping, recipe, emit)
                print(f'{arm}: numerical checks complete', flush=True)
        unchanged(records)
        verification = verify_study(ROOT, manifest)
        write_json(out / 'main_input_verification_after.json', verification)
        valid = verification['valid'] and all(c['valid'] for c in cases)
        summary = {'valid': valid, 'n_arms': len(prepared), 'n_cases': len(cases),
                   'n_skipped': sum(c.get('skipped', False) for c in cases),
                   'n_solver_runs': sum(len(c.get('solves', [])) for c in cases),
                   'failed_cases': [c for c in cases if not c['valid']],
                   'limitations': provenance['scope']}
        write_json(out / 'summary.json', summary)
        print(json.dumps(summary, indent=2), flush=True)
        return 0 if valid else 1
    except Exception as error:
        write_json(out / 'failure.json', {'type': type(error).__name__, 'message': str(error), 'completed_cases': len(cases)})
        raise


if __name__ == '__main__':
    raise SystemExit(main())
