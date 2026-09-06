"""Compare the v3 metadata-isolation revision against frozen v2 and an oracle.

All 156 declared conditions and five withdrawal controls; no optimization.
Copies and hashes this comparator, both implementations and their tests before
model preparation. The old experiment declarations/results remain unchanged.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import importlib
import json
import logging
from pathlib import Path
import sys
import time
import warnings

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import verify_medium_preparation as V


SOURCES = ['scripts/compare_medium_metadata_isolation.py', 'scripts/verify_medium_preparation.py',
           'gembench/medium_preparation_v2.py', 'gembench/medium_preparation_v3.py',
           'tests/test_medium_preparation_v2.py', 'tests/test_medium_preparation_v3.py',
           'tests/test_verify_medium_preparation.py',
           V.SOURCE_LOADER]


def compare_prepared(old, new, oracle_model, oracle_report):
    a, b, c = (V.canonical(model) for model in (old.model, new.model, oracle_model))
    V.require(a == b == c, 'Metadata-isolation revision or independent oracle differs in LP/GPR/chemistry')
    V.require(V.json_digest(a) == V.json_digest(b) == V.json_digest(c),
              'Metadata-only revision changed exact numeric serialization')
    # The independent artifact oracle uses the saved JSON schema; dataclass
    # tuple defaults become arrays when the live reports are serialized.
    V.check_report(V.loads(json.dumps(old.report, allow_nan=False)), oracle_report)
    V.check_report(V.loads(json.dumps(new.report, allow_nan=False)), oracle_report)
    V.require(old.report == new.report, 'Metadata-isolation revision changed the preparation report')
    return V.json_digest(b)


def compare(plan, source):
    old_api = importlib.import_module('gembench.medium_preparation_v2')
    new_api = importlib.import_module('gembench.medium_preparation_v3')
    cases, withdrawals, inventory = [], [], []
    reports = 0
    for item in plan['models']:
        base, provenance = source.build_parent(item['label'])
        conditions = source.metadata_conditions(item['organism'])
        media = sorted({condition.media for condition in conditions})
        completion = [V.base_medium(name) for name in media]
        policy = V.policy_record(item['compartments'])
        old_policy = old_api.CompartmentPolicy(**item['compartments'])
        new_policy = new_api.CompartmentPolicy(**item['compartments'])
        def prepare(api, model, medium, api_policy):
            return api.prepare_medium(model, medium, policy=api_policy,
                                      completion_media=completion, missing_policy='report')
        for condition in conditions:
            inventory.append((item['label'], condition.key))
            medium = V.medium_for(condition)
            old = prepare(old_api, base, medium, old_policy)
            new = prepare(new_api, base, medium, new_policy)
            oracle, report = V.independent_target(base, medium, media, policy)
            digest = compare_prepared(old, new, oracle, report)
            cases.append({'model': item['label'], 'condition': condition.key,
                          'v2_v3_oracle_signature_sha256': digest, 'reports_exactly_equal': True,
                          'boundary_reactions_checked': len(report['boundary_reactions'])})
            reports += 2
        reference = next((c for c in conditions if 'EX_glc__D_e' in c.exchanges), conditions[0])
        ref = V.medium_for(reference)
        _, first_report = V.independent_target(base, ref, media, policy)
        added = first_report['added_exchanges']
        rich = V.Medium('completion_withdrawal_start', 'Synthetic bound-reset stress test',
                        {**ref.uptakes, **{rid: -.001 for rid in added}})
        target = V.Medium('completion_withdrawal_target', 'Remove all newly completed components',
                          {rid: lower for rid, lower in ref.uptakes.items() if rid not in added})
        old_start = prepare(old_api, base, rich, old_policy)
        new_start = prepare(new_api, base, rich, new_policy)
        oracle_start, _ = V.independent_target(base, rich, media, policy)
        old_reuse = prepare(old_api, old_start.model, target, old_policy)
        new_reuse = prepare(new_api, new_start.model, target, new_policy)
        oracle_reuse, reuse_report = V.independent_target(oracle_start, target, media, policy)
        reuse_digest = compare_prepared(old_reuse, new_reuse, oracle_reuse, reuse_report)
        old_fresh = prepare(old_api, base, target, old_policy)
        new_fresh = prepare(new_api, base, target, new_policy)
        oracle_fresh, fresh_report = V.independent_target(base, target, media, policy)
        fresh_digest = compare_prepared(old_fresh, new_fresh, oracle_fresh, fresh_report)
        V.require(reuse_digest == fresh_digest, 'Withdrawal still depends on previous medium')
        closed = {rid: list(new_reuse.model.reactions.get_by_id(rid).bounds) for rid in added}
        V.require(all(bounds == [0., 1000.] for bounds in closed.values()), 'Withdrawn component remains open')
        withdrawals.append({'model': item['label'], 'reference_condition': reference.key,
            'reuse_fresh_v2_v3_oracle_signature_sha256': fresh_digest,
            'withdrawn_exchange_bounds': closed, 'fresh_and_reuse_reports_match_v2_and_oracle': True})
        reports += 4
    V.require(len(inventory) == len(set(inventory)) == plan['case_count_from_source_inventory'] == 156,
              'Case inventory is not the complete declared panel')
    return {'status': 'passed', 'cases_compared': len(cases), 'withdrawal_controls_compared': len(withdrawals),
            'version_reports_checked_against_oracle': reports, 'cases': cases, 'withdrawal_controls': withdrawals,
            'scope': 'Algebra/GPR/chemistry and preparation-report identity for the metadata-isolation revision; no optimization, fitness rescore, or biological validation. Metadata nonaliasing is tested separately in the snapshotted v3 tests.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=ROOT / 'results/medium_preparation_2026_09_06/metadata_isolation_comparison')
    args = parser.parse_args()
    out = args.out
    out.mkdir(parents=True, exist_ok=False)
    started = time.time()
    try:
        snapshots = out / 'source_snapshots'
        snapshots.mkdir()
        input_hashes = {}
        for relative in SOURCES:
            path = V.resolve_input(ROOT, relative)
            destination = snapshots / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(path.read_bytes())
            input_hashes[relative] = V.sha(destination)
        manifest_path = ROOT / 'results/medium_preparation_2026_09_06/runs/main_v2/manifest.json'
        manifest = V.read(manifest_path)
        frozen = V.verify_manifest(ROOT, manifest)
        input_hashes[str(manifest_path.relative_to(ROOT))] = V.sha(manifest_path)
        (out / 'audit_inputs.json').write_text(json.dumps({'created_utc': datetime.now(timezone.utc).isoformat(),
            'additional_input_sha256': input_hashes, 'frozen_input_sha256': frozen,
            'study_fingerprint': manifest['content_fingerprint'],
            'scope': 'Additive structural comparison; all sources fixed before applying either version to real models.'}, indent=2) + '\n')
        logging.getLogger('cobra').setLevel(logging.ERROR)
        warnings.filterwarnings('ignore', category=FutureWarning, module=r'cobra\.medium\.boundary_types')
        def forbidden(*args, **kwargs):
            raise AssertionError('No optimization in metadata-isolation comparison')
        V.cobra.Model.optimize = forbidden
        V.cobra.Model.slim_optimize = forbidden
        plan = manifest['plan']['recipe']
        V.require(plan['study_id'] == 'explicit-medium-preparation-v2' and len(plan['models']) == len(V.LABELS)
                  and {m['label'] for m in plan['models']} == V.LABELS, 'Unexpected panel declaration')
        result = compare(plan, V.load_source())
        V.require(frozen == V.verify_manifest(ROOT, manifest), 'Frozen inputs changed during comparison')
        V.require(input_hashes == {path: V.sha(V.resolve_input(ROOT, path)) for path in input_hashes},
                  'Comparator or implementation inputs changed during comparison')
        result.update(additional_input_sha256=input_hashes, frozen_input_files_verified=len(frozen),
                      study_fingerprint=manifest['content_fingerprint'], elapsed_seconds=time.time() - started)
        (out / 'report.json').write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + '\n')
        print(json.dumps({key: result[key] for key in ('status', 'cases_compared', 'withdrawal_controls_compared',
                                                     'version_reports_checked_against_oracle')}))
    except Exception as error:
        (out / 'failure.json').write_text(json.dumps({'type': type(error).__name__, 'message': str(error),
            'elapsed_seconds': time.time() - started, 'scope': 'Additive comparison attempt retained; frozen source/results unchanged.'}, indent=2) + '\n')
        raise


if __name__ == '__main__':
    main()
