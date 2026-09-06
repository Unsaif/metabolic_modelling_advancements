"""Independently verify the additive 156-case medium structural report.

Uses the independent artifact oracle and historical source loader only, never
the new preparation helper, runner or structural auditor. No solver is called.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import logging
from pathlib import Path
import sys
import time
import warnings

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import verify_medium_preparation as V


def verify(run, artifact):
    direct_paths = [run / 'manifest.json', artifact / 'report.json', artifact / 'auditor_source.py',
                    Path(__file__).resolve(), Path(V.__file__).resolve()]
    initial = {str(path.relative_to(ROOT)): V.sha(path) for path in direct_paths}
    manifest, report = V.read(direct_paths[0]), V.read(direct_paths[1])
    frozen = V.verify_manifest(ROOT, manifest)
    plan = manifest['plan']['recipe']
    V.require(plan['study_id'] == 'explicit-medium-preparation-v2', 'Structural report must use the frozen v2 plan')
    V.require(report['input_verification']['valid'] and
              report['input_verification']['content_fingerprint'] == manifest['content_fingerprint'],
              'Structural audit input verification does not identify this declaration')
    for relative, digest in report['input_sha256'].items():
        V.require(V.sha(V.resolve_input(ROOT, relative)) == digest, 'Structural auditor input changed')
    V.require(report['input_sha256'] == {
        'scripts/audit_medium_structures.py': V.sha(artifact / 'auditor_source.py'),
        str((run / 'manifest.json').relative_to(ROOT)): V.sha(run / 'manifest.json')},
        'Structural auditor source snapshot or manifest dependency does not match')
    V.require(len(plan['models']) == len(V.LABELS) and {m['label'] for m in plan['models']} == V.LABELS,
              'Changed or duplicated model panel')
    cases = {(row['model'], row['condition']): row for row in report['cases']}
    withdrawals = {row['label']: row for row in report['withdrawal_checks']}
    V.require(len(cases) == len(report['cases']) and len(withdrawals) == len(report['withdrawal_checks']) == 5
              and set(withdrawals) == V.LABELS, 'Duplicate cases or missing withdrawal control')
    source = V.load_source()
    inventory, models, serialization = [], [], []
    boundary_records = 0
    for item in plan['models']:
        label = item['label']
        base, provenance = source.build_parent(label)
        V.require(all(frozen.get(path) == digest for path, digest in provenance['input_sha256'].items()),
                  'A parent reconstruction dependency was not frozen')
        conditions = source.metadata_conditions(item['organism'])
        media = sorted({condition.media for condition in conditions})
        policy = V.policy_record(item['compartments'])
        V.require(policy['exchange_metabolites'] == ({'EX_AEP_e': '2ameph_e'} if label == 'Putida_iJN1463' else {}),
                  'Unexpected exchange identity binding')
        legacy = base.copy()
        added = V.complete_medium_transport(legacy, media, policy['completion_exclude'])
        for reaction in legacy.exchanges:
            reaction.bounds = (0., 1000.)
        count = 0
        for condition in conditions:
            key = (label, condition.key)
            inventory.append(key)
            V.require(key in cases, 'Missing declared structural case')
            case = cases[key]
            medium = V.medium_for(condition)
            old = legacy.copy()
            V.apply_medium(old, medium)
            new, expected = V.independent_target(base, medium, media, policy)
            changed = V.verify_signatures({'legacy_signature_sha256': case['legacy_sha256'],
                'algebra_and_chemistry_sha256': case['new_sha256']}, V.canonical(old), V.canonical(new))
            if changed:
                serialization.append({'model': label, 'condition': condition.key,
                                      'equal_numeric_value_type_changes': changed})
            V.check_report(case['new_boundary_report'], expected)
            boundary_records += len(expected['boundary_reactions'])
            count += 1
        reference = next((c for c in conditions if 'EX_glc__D_e' in c.exchanges), conditions[0])
        medium = V.medium_for(reference)
        rich = V.Medium('completion_withdrawal_start', 'Synthetic bound-reset stress test',
                        {**medium.uptakes, **{rid: -.001 for rid in added}})
        target = V.Medium('completion_withdrawal_target', 'Remove all newly completed components',
                          {rid: lower for rid, lower in medium.uptakes.items() if rid not in added})
        old_reuse, old_fresh = legacy.copy(), legacy.copy()
        V.apply_medium(old_reuse, rich)
        V.apply_medium(old_reuse, target)
        V.apply_medium(old_fresh, target)
        new_start, _ = V.independent_target(base, rich, media, policy)
        new_reuse, reuse_report = V.independent_target(new_start, target, media, policy)
        new_fresh, fresh_report = V.independent_target(base, target, media, policy)
        a, b = V.canonical(new_reuse), V.canonical(new_fresh)
        V.require(a == b, 'Independent fresh/reused explicit preparation differs')
        observed = withdrawals[label]
        V.require(observed['strict_reuse_signature_sha256'] == V.json_digest(a) and
                  observed['strict_fresh_signature_sha256'] == V.json_digest(b), 'Withdrawal signature differs')
        V.check_report(observed['strict_reuse_report'], reuse_report)
        V.check_report(observed['strict_fresh_report'], fresh_report)
        boundary_records += len(reuse_report['boundary_reactions']) + len(fresh_report['boundary_reactions'])
        retained = {rid: list(old_reuse.reactions.get_by_id(rid).bounds) for rid in added
                    if old_reuse.reactions.get_by_id(rid).lower_bound < 0}
        V.require(observed['legacy_retained_uptakes'] == retained and
                  observed['legacy_reuse_vs_fresh'] == V.differences(V.canonical(old_fresh), V.canonical(old_reuse)),
                  'Legacy history difference does not match independent reconstruction')
        closed = {rid: list(new_reuse.reactions.get_by_id(rid).bounds) for rid in added}
        V.require(all(bounds == [0., 1000.] for bounds in closed.values()), 'Withdrawal leaves an added nutrient open')
        models.append({'label': label, 'fresh_cases': count, 'legacy_retained_uptakes': retained,
                       'explicit_withdrawn_bounds': closed, 'fresh_and_reuse_exactly_equal': True,
                       'parent_source_sha256': provenance['source_sha256']})
    V.require([(r['model'], r['condition']) for r in report['cases']] == inventory and
              len(inventory) == report['n_cases'] == plan['case_count_from_source_inventory'] == 156 and
              report['all_exactly_equal'] is True, 'Structural inventory or completion claim differs')
    V.require(frozen == V.verify_manifest(ROOT, manifest), 'Inputs changed during independent verification')
    V.require(initial == {str(path.relative_to(ROOT)): V.sha(path) for path in direct_paths},
              'Artifacts or checker sources changed during verification')
    return {'status': 'passed', 'cases_verified': len(inventory), 'withdrawal_controls_verified': 5,
            'boundary_reports_verified': len(inventory) + 10, 'boundary_records_verified': boundary_records,
            'models': models, 'numeric_serialization_differences': serialization,
            'frozen_input_files_verified': len(frozen), 'study_fingerprint': manifest['content_fingerprint'],
            'frozen_input_sha256': frozen, 'direct_input_sha256': initial,
            'method': 'Independently reconstructed legacy and explicit preparation; each saved byte signature checked separately, numeric dictionaries compared exactly, and all fresh/reuse boundary reports checked. No optimizer or new preparation code imported.',
            'scope': 'Structural conformance and synthetic history reset only. The original numerical failures remain failures; this does not certify the unsaved curated solutions, biological transport, or mutant fitness.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, default=ROOT / 'results/medium_preparation_2026_09_06/runs/main_v2')
    parser.add_argument('--artifact', type=Path, default=ROOT / 'results/medium_preparation_2026_09_06/structural_completion')
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    run, artifact = args.run.resolve(), args.artifact.resolve()
    V.require(run.is_relative_to(ROOT) and artifact.is_relative_to(ROOT), 'Inputs must be in the repository')
    out = args.out or artifact / 'independent_verification'
    out.mkdir(parents=True, exist_ok=False)
    for source, name in [(Path(__file__).resolve(), 'checker_source.py'), (Path(V.__file__).resolve(), 'oracle_source.py')]:
        (out / name).write_bytes(source.read_bytes())
    (out / 'attempt.json').write_text(json.dumps({'started_utc': datetime.now(timezone.utc).isoformat(),
        'scope': 'Independent additive checker; sources snapshotted before reading the artifact.'}, indent=2) + '\n')
    started = time.time()
    logging.getLogger('cobra').setLevel(logging.ERROR)
    warnings.filterwarnings('ignore', category=FutureWarning, module=r'cobra\.medium\.boundary_types')
    def forbidden(*args, **kwargs):
        raise AssertionError('Independent structural verification must not optimize')
    V.cobra.Model.optimize = forbidden
    V.cobra.Model.slim_optimize = forbidden
    try:
        result = verify(run, artifact)
        result['elapsed_seconds'] = time.time() - started
        (out / 'result.json').write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + '\n')
        print(json.dumps({key: result[key] for key in ('status', 'cases_verified', 'withdrawal_controls_verified',
                                                     'boundary_reports_verified', 'boundary_records_verified')}))
    except Exception as error:
        (out / 'failure.json').write_text(json.dumps({'type': type(error).__name__, 'message': str(error),
            'elapsed_seconds': time.time() - started, 'scope': 'Checker attempt retained; primary artifacts unchanged.'}, indent=2) + '\n')
        raise


if __name__ == '__main__':
    main()
