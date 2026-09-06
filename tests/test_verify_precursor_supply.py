"""Synthetic artifact checks: no optimization or production solver helper."""
from copy import deepcopy
import hashlib
import json

import pytest

from scripts import verify_precursor_supply as V


@pytest.fixture
def case():
    data = [
        ('EX_carbon', {'carbon_c': -1.}, (-10., 1000.), ''),
        ('CHRPL', {'carbon_c': -1., '4hbz_c': 1.}, (0., 1000.), 'g_chor'),
        ('HBZOPT', {'4hbz_c': -1., 'q8h2_c': 1.}, (0., 1000.), 'g_ubi'),
        ('UHBZ1t_pp', {'4hbz_e': -1., '4hbz_c': 1.}, (-1000., 1000.), 'g_tr'),
        (V.EX, {'4hbz_e': -1.}, (0., 1000.), ''),
        (V.COUM, {'T4hcinnm_e': -1.}, (0., 1000.), ''),
        ('4HBHYOX', {'4hbz_c': -1., 'carbon_c': 1.}, (0., 1000.), ''),
        ('sink_2ohph_c', {'4hbz_c': -1.}, (0., 1000.), ''),
        ('Growth', {'carbon_c': -1., 'q8h2_c': -.125}, (0., 1000.), ''),
        ('QRED', {'q8_c': -1., 'q8h2_c': 1.}, (0., 1000.), ''),
        ('QOX', {'q8_c': 1., 'q8h2_c': -1.}, (0., 1000.), '')]
    records = {rid: {'stoich': stoich, 'bounds': bounds, 'gpr': rule} for rid, stoich, bounds, rule in data}
    bounds = {rid: record['bounds'] for rid, record in records.items()}
    fluxes = dict.fromkeys(records, 0.)
    fluxes.update(EX_carbon=-1.125, CHRPL=.125, HBZOPT=.125, Growth=1.)
    answers = []
    for solver, cycle in [('glpk', 1.), ('highs', 2.)]:
        values = {**fluxes, 'QRED': cycle, 'QOX': cycle}
        answers.append({'solver': solver, 'status': 'optimal', 'objective': 1., 'growth_flux': 1.,
            'feasibility_tolerance': 1e-9, 'full_primal_fluxes': values,
            'primal_certificate': {'max_abs_S_residual': 0., 'max_bound_violation': 0.,
                                  'n_S_rows_over_1e-6': 0, 'n_bounds_over_1e-6': 0},
            'quinone_balance_residuals': {'q8_c': 0., 'q8h2_c': 0.}})
    expected = {'id': 'wild_type', 'kind': 'growth', 'condition': 'toy'}
    witnesses = ['Growth', 'CHRPL', V.EX, V.COUM, '4HBHYOX', 'sink_2ohph_c']
    row = {'specification': dict(expected), 'solves': answers, 'objective': 1., 'growth': 1.,
           'selected_glpk_fluxes': {rid: fluxes[rid] for rid in witnesses}, 'precursor_balance_residual': 0.}
    recipe = {'solver_tolerance': 1e-9, 'objective_tolerance': 1e-8, 'residual_tolerance': 1e-8,
              'witness_reactions': witnesses}
    return row, expected, records, bounds, {'PP_5317': 'g_chor'}, recipe, .125


def test_feasible_distinct_solver_vectors_and_no_legacy_diagnostic_field(case):
    result = V.verify_case(*case)
    assert result['objective_disagreement'] == 0
    assert len(result['solvers']) == 2
    assert case[0]['solves'][0]['full_primal_fluxes'] != case[0]['solves'][1]['full_primal_fluxes']


@pytest.mark.parametrize('corruption', [
    'missing_reaction', 'extra_reaction', 'nonfinite', 'bool_flux', 'infeasible_flux',
    'false_certificate', 'false_count', 'negative_certificate', 'objective', 'growth',
    'top_objective', 'selected', 'quinone', 'precursor', 'status', 'solver', 'tolerance',
    'specification', 'unsupported_diagnostic'])
def test_inconsistent_artifacts_fail(case, corruption):
    row = case[0]
    solve = row['solves'][1]
    fluxes = solve['full_primal_fluxes']
    if corruption == 'missing_reaction':
        del fluxes['QOX']
    elif corruption == 'extra_reaction':
        fluxes['unapproved'] = 0.
    elif corruption == 'nonfinite':
        fluxes['QOX'] = float('nan')
    elif corruption == 'bool_flux':
        fluxes['QOX'] = True
    elif corruption == 'infeasible_flux':
        fluxes['HBZOPT'] += .01
    elif corruption == 'false_certificate':
        solve['primal_certificate']['max_abs_S_residual'] = 1e-9
    elif corruption == 'false_count':
        solve['primal_certificate']['n_S_rows_over_1e-6'] = True
    elif corruption == 'negative_certificate':
        solve['primal_certificate']['max_bound_violation'] = -1.
    elif corruption == 'objective':
        solve['objective'] = .5
    elif corruption == 'growth':
        solve['growth_flux'] = .5
    elif corruption == 'top_objective':
        row['objective'] = .5
    elif corruption == 'selected':
        row['selected_glpk_fluxes']['CHRPL'] = 0.
    elif corruption == 'quinone':
        solve['quinone_balance_residuals']['q8_c'] = .001
    elif corruption == 'precursor':
        row['precursor_balance_residual'] = .001
    elif corruption == 'status':
        solve['status'] = 'infeasible'
    elif corruption == 'solver':
        solve['solver'] = 'glpk'
    elif corruption == 'tolerance':
        solve['feasibility_tolerance'] = 1e-7
    elif corruption == 'specification':
        row['specification']['close_reactions'] = ['CHRPL']
    elif corruption == 'unsupported_diagnostic':
        case[2]['DIAG_unapproved'] = {'stoich': {}, 'bounds': (0., 0.), 'gpr': ''}
        case[3]['DIAG_unapproved'] = (0., 0.)
        for answer in row['solves']:
            answer['full_primal_fluxes']['DIAG_unapproved'] = 0.
    with pytest.raises(ValueError):
        V.verify_case(*case)


@pytest.mark.parametrize('restriction', ['gene', 'reaction', 'growth', 'exchange'])
def test_full_vectors_are_checked_against_case_specific_bounds(case, restriction):
    row, expected = case[:2]
    if restriction == 'gene':
        expected['loci'] = ['PP_5317']
    elif restriction == 'reaction':
        expected['close_reactions'] = ['CHRPL']
    elif restriction == 'growth':
        expected['growth_bounds'] = [.5, .5]
    elif restriction == 'exchange':
        expected['exchange_bounds'] = [-1., -.1]
    row['specification'] = deepcopy(expected)
    with pytest.raises(ValueError, match='Certificate|feasibility'):
        V.verify_case(*case)


def test_exact_row_sum_detects_artificial_source_and_direction_change(case):
    records = case[2]
    assert V.exact_pool(records)[0] == .125
    records['unapproved'] = {'stoich': {'4hbz_c': 1.}, 'bounds': (0., 1000.), 'gpr': ''}
    with pytest.raises(ValueError, match='row combination'):
        V.exact_pool(records)
    del records['unapproved']
    records['4HBHYOX']['bounds'] = (-1000., 1000.)
    with pytest.raises(ValueError, match='nonnegative'):
        V.exact_pool(records)


def test_gene_logic_and_undeclared_expressions():
    assert V.gpr_active('(a and b) or c', {'a'})
    assert not V.gpr_active('(a and b) or c', {'b', 'c'})
    assert V.gpr_active('', {'a'})
    with pytest.raises(ValueError, match='Unsupported'):
        V.gpr_active('not a', {'a'})


def test_medium_completion_is_uptake_only_and_respects_exclusions():
    base = {'Growth': {'stoich': {}, 'bounds': (0., 1000.), 'gpr': ''},
            'ATPM': {'stoich': {}, 'bounds': (2., 1000.), 'gpr': ''}}
    media = {'MOPS': {'bigg_ids': 'na1;btn;fol;zn2;glc__D', 'aerobic': 'yes', 'trace_components': 'btn'}}
    condition = {'key': 'Glucose | MOPS', 'media': 'MOPS', 'exchanges': ['EX_glc__D_e']}
    records, bounds, context = V.prepare_network(base, dict.fromkeys(['na1_c', 'btn_c', 'fol_c', 'zn2_c', 'glc__D_c'], 'c'),
        set(), [condition], media, condition, {'protocol': {'carbon_uptake': -10.}})
    assert records['MEDt_btn']['bounds'] == (0., 1000.)
    assert records['MEDt_btn']['stoich'] == {'btn_e': -1., 'btn_c': 1.}
    assert bounds['EX_btn_e'] == (-.001, 1000.)
    assert bounds['EX_zn2_e'] == (-.1, 1000.)
    assert bounds['EX_glc__D_e'] == (-10., 1000.)
    assert context['missing_medium'] == ['EX_fol_e', 'EX_o2_e']
    assert base.keys() == {'Growth', 'ATPM'}


def test_mixed_external_compartments_remain_in_full_lp_even_if_context_omits_them():
    base = {'Growth': {'stoich': {}, 'bounds': (0., 1000.), 'gpr': ''},
            'ATPM': {'stoich': {}, 'bounds': (2., 1000.), 'gpr': ''},
            'EX_na1_e': {'stoich': {'na1_e': -1.}, 'bounds': (0., 1000.), 'gpr': ''}}
    media = {'MOPS': {'bigg_ids': 'na1;btn', 'aerobic': 'no', 'trace_components': 'btn'}}
    condition = {'key': 'toy', 'media': 'MOPS', 'exchanges': []}
    records, bounds, context = V.prepare_network(base, {'na1_e': 'C_e', 'na1_c': 'C_c', 'btn_c': 'C_c'},
        {'EX_na1_e'}, [condition], media, condition, {'protocol': {'carbon_uptake': -10.}})
    assert context['medium_completion_added'] == ['EX_btn_e']
    assert 'EX_btn_e' not in context['exchange_bounds']
    assert records['EX_btn_e']['stoich'] == {'btn_e': -1.}
    assert bounds['EX_btn_e'] == (-.001, 1000.)
    assert bounds['MEDt_btn'] == (0., 1000.)


def manifest_for(root):
    path = root / 'input.txt'
    path.write_text('source\n')
    manifest = {'schema_version': 1, 'manifest_type': 'test', 'plan': {'paths': ['input.txt']},
                'files': [{'path': 'input.txt', 'sha256': V.sha(path), 'size_bytes': path.stat().st_size}]}
    manifest['content_fingerprint'] = hashlib.sha256(json.dumps(manifest, sort_keys=True,
        separators=(',', ':'), ensure_ascii=False, allow_nan=False).encode()).hexdigest()
    return manifest


def test_input_hashes_allow_relocation_and_reject_tampering(tmp_path):
    left, right = tmp_path / 'left', tmp_path / 'right'
    left.mkdir()
    right.mkdir()
    manifest = manifest_for(left)
    (right / 'input.txt').write_bytes((left / 'input.txt').read_bytes())
    assert V.verify_manifest(right, manifest)['valid']
    (right / 'input.txt').write_text('change\n')
    with pytest.raises(ValueError, match='Frozen input changed'):
        V.verify_manifest(right, manifest)
    manifest['plan']['paths'].append('input.txt')
    with pytest.raises(ValueError, match='inventory'):
        V.verify_manifest(left, manifest)


@pytest.mark.parametrize('content', ['{"x":1,"x":2}', '{"x":NaN}', '{"x":Infinity}'])
def test_ambiguous_or_nonfinite_json_is_rejected(tmp_path, content):
    path = tmp_path / 'bad.json'
    path.write_text(content)
    with pytest.raises(ValueError):
        V.read_json(path)


def test_snapshot_precedes_first_verification_and_failure_is_preserved(tmp_path, monkeypatch):
    study, out = tmp_path / 'study', tmp_path / 'verification'
    study.mkdir()
    for name in ('manifest.json', 'summary.json', 'contexts.json', 'cases.jsonl.gz', 'input_verification_after_run.json'):
        (study / name).write_text('{}')
    def fail(root, actual_study):
        assert actual_study == study
        assert (out / 'verifier_source.py').read_bytes() == V.Path(V.__file__).read_bytes()
        assert V.read_json(out / 'verifier_provenance.json')['source_sha256'] == V.sha(out / 'verifier_source.py')
        assert len(V.read_json(out / 'result_inputs.json')) == 5
        raise ValueError('synthetic disagreement')
    monkeypatch.setattr(V, 'verify_results', fail)
    with pytest.raises(ValueError, match='synthetic disagreement'):
        V.main(['--study', str(study), '--out', str(out)])
    assert V.read_json(out / 'failure.json')['message'] == 'synthetic disagreement'
    with pytest.raises(FileExistsError):
        V.main(['--study', str(study), '--out', str(out)])
