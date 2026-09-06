"""Adversarial checks for the independent saved-artifact verifier; no solves."""
from copy import deepcopy
import hashlib
import importlib.util
import json
from pathlib import Path

import cobra
import pytest


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('medium_artifact_verifier', ROOT / 'scripts/verify_medium_preparation.py')
V = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(V)


def toy():
    model = cobra.Model('artifact')
    ext = cobra.Metabolite('a_e', compartment='e', formula='C', charge=0)
    cyt = cobra.Metabolite('a_c', compartment='c', formula='C', charge=0)
    exchange = cobra.Reaction('EX_a_e', lower_bound=-10, upper_bound=1000)
    exchange.add_metabolites({ext: -1})
    carrier = cobra.Reaction('T', lower_bound=0, upper_bound=1000)
    carrier.add_metabolites({ext: -1, cyt: 1})
    growth = cobra.Reaction('Growth', lower_bound=0, upper_bound=1000)
    growth.add_metabolites({cyt: -1})
    growth.gene_reaction_rule = 'g1'
    source = cobra.Reaction('sink_a_c', lower_bound=-1, upper_bound=1000)
    source.add_metabolites({cyt: -1})
    model.add_reactions([exchange, carrier, growth, source])
    model.objective = growth
    return model


def policy(**kwargs):
    return V.policy_record({'external': 'e', 'cytoplasm': 'c', **kwargs})


@pytest.fixture(autouse=True)
def no_optimizer(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError('No optimization is allowed in artifact-verifier tests')
    monkeypatch.setattr(cobra.Model, 'optimize', forbidden)
    monkeypatch.setattr(cobra.Model, 'slim_optimize', forbidden)


@pytest.mark.parametrize('raw', ['{"x": 1, "x": 2}', '{"x": NaN}', '{"x": Infinity}'])
def test_untrustworthy_json_rejected(raw):
    with pytest.raises(ValueError):
        V.loads(raw)


def test_manifest_detects_missing_record_tampered_digest_and_path(tmp_path):
    data = b'published input\n'
    (tmp_path / 'a.txt').write_bytes(data)
    record = {'path': 'a.txt', 'size_bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}
    body = {'schema_version': 1, 'manifest_type': 'gembench-study-freeze',
            'plan': {'paths': ['a.txt']}, 'files': [record]}
    manifest = {**body, 'content_fingerprint': V.json_digest(body, unicode=True)}
    assert V.verify_manifest(tmp_path, manifest) == {'a.txt': record['sha256']}
    for change in ('missing', 'duplicate', 'digest', 'escape', 'fingerprint', 'size'):
        broken = deepcopy(manifest)
        if change == 'missing':
            broken['files'] = []
        elif change == 'duplicate':
            broken['files'].append(deepcopy(record))
        elif change == 'digest':
            broken['files'][0]['sha256'] = '0' * 64
        elif change == 'escape':
            broken['files'][0]['path'] = '../a.txt'
        elif change == 'size':
            broken['files'][0]['size_bytes'] = True
        else:
            broken['content_fingerprint'] = '0' * 64
        with pytest.raises((ValueError, FileNotFoundError)):
            V.verify_manifest(tmp_path, broken)
    (tmp_path / 'a.txt').write_bytes(b'changed')
    with pytest.raises(ValueError, match='Changed input'):
        V.verify_manifest(tmp_path, manifest)


def test_exchange_identity_exception_requires_exact_existing_pair():
    model = toy()
    model.reactions.EX_a_e.id = 'EX_alias_e'
    with pytest.raises(ValueError, match='identity/sign'):
        V.boundary_record(model, policy())
    report = V.boundary_record(model, policy(exchange_metabolites={'EX_alias_e': 'a_e'}))
    assert report['exchange_ids'] == ['EX_alias_e']
    assert report['internal_supply_capable_boundaries'] == ['sink_a_c']
    with pytest.raises(ValueError, match='identity/sign'):
        V.boundary_record(model, policy(exchange_metabolites={'EX_alias_e': 'wrong_e'}))
    with pytest.raises(ValueError, match='Unknown exchange'):
        V.boundary_record(model, policy(exchange_metabolites={'EX_absent_e': 'a_e'}))
    model.reactions.EX_alias_e.add_metabolites({model.metabolites.a_e: 2})
    with pytest.raises(ValueError, match='identity/sign'):
        V.boundary_record(model, policy(exchange_metabolites={'EX_alias_e': 'a_e'}))


def test_explicit_reset_reports_every_boundary_and_preserves_internal_source():
    source = toy()
    source.reactions.EX_a_e.id = 'EX_alias_e'
    config = policy(exchange_metabolites={'EX_alias_e': 'a_e'})
    medium = V.Medium('withdrawal', '', {})
    target, report = V.independent_target(source, medium, [], config)
    assert source.reactions.EX_alias_e.bounds == (-10, 1000)
    assert target.reactions.EX_alias_e.bounds == (0, 1000)
    assert target.reactions.sink_a_c.bounds == (-1, 1000)
    assert set(report['boundary_reactions']) == {'EX_alias_e', 'Growth', 'sink_a_c'}
    observed = {**report, 'scope': 'No biological transport validation'}
    V.check_report(observed, report)
    del observed['boundary_reactions']['sink_a_c']
    # Do not alias the expected dict when simulating a tampered artifact.
    _, expected = V.independent_target(source, medium, [], config)
    with pytest.raises(ValueError, match='boundary_reactions'):
        V.check_report(observed, expected)


@pytest.mark.parametrize('field', ['stoichiometry', 'bounds', 'gpr', 'formula', 'charge', 'objective'])
def test_signature_detects_mechanistically_material_changes(field):
    model = toy()
    original = V.json_digest(V.canonical(model))
    if field == 'stoichiometry':
        model.reactions.T.add_metabolites({model.metabolites.a_c: 1})
    elif field == 'bounds':
        model.reactions.T.upper_bound = 9
    elif field == 'gpr':
        model.reactions.Growth.gene_reaction_rule = 'g2'
    elif field == 'formula':
        model.metabolites.a_c.formula = 'C2'
    elif field == 'charge':
        model.metabolites.a_c.charge = 1
    else:
        model.objective = model.reactions.T
    assert V.json_digest(V.canonical(model)) != original


def witness():
    return {'solver': 'glpk', 'status': 'optimal', 'objective': 10.,
            'full_primal_fluxes': {'EX_a_e': -10., 'T': 10., 'Growth': 10., 'sink_a_c': 0.},
            'primal_certificate': {'max_abs_S_residual': 0., 'max_bound_violation': 0.,
                                  'n_S_rows_over_1e-6': 0, 'n_bounds_over_1e-6': 0}}


@pytest.mark.parametrize('failure', ['missing', 'nan', 'imbalance', 'bounds', 'objective', 'certificate'])
def test_saved_primal_failures_are_rejected(failure):
    structure, solve = V.canonical(toy()), witness()
    assert V.primal_certificate(structure, solve, 1e-8)['objective_recomputed'] == 10
    if failure == 'missing':
        del solve['full_primal_fluxes']['T']
    elif failure == 'nan':
        solve['full_primal_fluxes']['T'] = float('nan')
    elif failure == 'imbalance':
        solve['full_primal_fluxes']['T'] = 9
    elif failure == 'bounds':
        structure['reactions']['EX_a_e']['bounds'][0] = -9
    elif failure == 'objective':
        solve['objective'] = 0
    else:
        solve['primal_certificate']['max_abs_S_residual'] = 1
    with pytest.raises(ValueError):
        V.primal_certificate(structure, solve, 1e-8)


def test_infeasible_is_not_zero_or_an_independently_certified_primal():
    structure = V.canonical(toy())
    assert V.primal_certificate(structure, {'status': 'infeasible', 'objective': None}, 1e-8) is None
    for solve in ({'status': 'infeasible', 'objective': 0},
                  {'status': 'infeasible', 'objective': None, 'full_primal_fluxes': {}},
                  {'status': 'undefined', 'objective': None}):
        with pytest.raises(ValueError):
            V.primal_certificate(structure, solve, 1e-8)


def test_partial_prefix_rejects_missing_reordered_duplicate_or_complete_cases():
    inventory = [('A', 'one'), ('A', 'two'), ('B', 'three')]
    def rows(keys):
        return [{'model': label, 'condition': condition} for label, condition in keys]
    assert V.strict_prefix(rows(inventory[:2]), inventory) == set(inventory[:2])
    for keys in ([], inventory, inventory[1:2], inventory[1::-1], inventory[:1] * 2):
        with pytest.raises(ValueError):
            V.strict_prefix(rows(keys), inventory)


def test_failed_checker_snapshots_source_and_refuses_overwrite(tmp_path, monkeypatch):
    out = tmp_path / 'attempt'
    monkeypatch.setattr(V.sys, 'argv', ['verify', '--partial', '--out', str(out)])
    def fail(run):
        assert (out / 'checker_source.py').read_bytes() == Path(V.__file__).read_bytes()
        raise ValueError('Toy artifact refusal')
    monkeypatch.setattr(V, 'verify_partial', fail)
    with pytest.raises(ValueError, match='Toy artifact refusal'):
        V.main()
    before = {p.name: p.read_bytes() for p in out.iterdir()}
    assert json.loads(before['failure.json'])['message'] == 'Toy artifact refusal'
    with pytest.raises(FileExistsError):
        V.main()
    assert before == {p.name: p.read_bytes() for p in out.iterdir()}


def test_cli_normalizes_relative_run_path_before_reading(tmp_path, monkeypatch):
    relative = Path('results/medium_preparation_2026_09_06/runs/main')
    monkeypatch.chdir(ROOT)
    monkeypatch.setattr(V.sys, 'argv', ['verify', '--partial', '--run', str(relative), '--out', str(tmp_path / 'attempt')])
    def check(run):
        assert run == ROOT / relative and run.is_absolute()
        return {'status': 'passed_partial', 'cases_verified': 0, 'finite_primal_vectors_verified': 0,
                'infeasible_solver_reports_retained': 0, 'maxima': {}}
    monkeypatch.setattr(V, 'verify_partial', check)
    V.main()


def test_each_byte_signature_checked_without_conflating_int_float_json():
    old = {'bounds': [0, 1000]}
    new = {'bounds': [0., 1000.]}
    assert old == new and V.json_digest(old) != V.json_digest(new)
    case = {'legacy_signature_sha256': V.json_digest(old), 'algebra_and_chemistry_sha256': V.json_digest(new)}
    changes = V.verify_signatures(case, old, new)
    assert [row['path'] for row in changes] == [['bounds', 0], ['bounds', 1]]
    with pytest.raises(ValueError, match='prepared signature'):
        V.verify_signatures({**case, 'algebra_and_chemistry_sha256': V.json_digest(old)}, old, new)
    with pytest.raises(ValueError, match='numeric LP'):
        V.verify_signatures(case, old, {'bounds': [0., 999.999999999]})


def test_structural_checker_snapshots_both_sources_and_retains_failure(tmp_path, monkeypatch):
    from scripts import verify_medium_structures as structural
    out = tmp_path / 'structural_attempt'
    monkeypatch.setattr(structural.sys, 'argv', ['structural', '--out', str(out)])
    def fail(run, artifact):
        assert (out / 'checker_source.py').read_bytes() == Path(structural.__file__).read_bytes()
        assert (out / 'oracle_source.py').read_bytes() == Path(structural.V.__file__).read_bytes()
        raise ValueError('Toy structural refusal')
    monkeypatch.setattr(structural, 'verify', fail)
    with pytest.raises(ValueError, match='Toy structural refusal'):
        structural.main()
    before = {p.name: p.read_bytes() for p in out.iterdir()}
    with pytest.raises(FileExistsError):
        structural.main()
    assert before == {p.name: p.read_bytes() for p in out.iterdir()}


@pytest.mark.parametrize('change', ['bounds', 'gpr', 'report'])
def test_metadata_comparison_detects_nonmetadata_regressions(change):
    from types import SimpleNamespace
    from scripts.compare_medium_metadata_isolation import compare_prepared
    oracle, expected = V.independent_target(toy(), V.Medium('closed', '', {}), [], policy())
    report = {**expected, 'scope': 'No biological transport validation'}
    old = SimpleNamespace(model=oracle.copy(), report=deepcopy(report))
    new = SimpleNamespace(model=oracle.copy(), report=deepcopy(report))
    # Live dataclass defaults serialize tuples as JSON arrays in artifacts.
    old.report['policy']['completion_exclude'] = tuple(old.report['policy']['completion_exclude'])
    new.report['policy']['completion_exclude'] = tuple(new.report['policy']['completion_exclude'])
    assert compare_prepared(old, new, oracle, expected) == V.json_digest(V.canonical(oracle))
    if change == 'bounds':
        new.model.reactions.T.upper_bound -= 1
    elif change == 'gpr':
        new.model.reactions.Growth.gene_reaction_rule = 'g_other'
    else:
        del new.report['boundary_reactions']['sink_a_c']
    with pytest.raises(ValueError):
        compare_prepared(old, new, oracle, expected)
