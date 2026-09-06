"""Finish the declared structural comparison separately from numerical solving.

Additive post-failure analysis. No solver calls or numeric fitness reads.
The numerical failures remain failures; this does not certify their solutions.
"""
import argparse
import json
from pathlib import Path
import shutil
import sys

import cobra

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from gembench.cards import sha256_of
from gembench.study import verify_study
from gembench.fitness_browser import base_medium
from gembench.media import Medium, apply_medium
from gembench.protocols.carbon_fitness_generic import complete_medium_transport
from gembench.medium_preparation_v2 import CompartmentPolicy, prepare_medium
from scripts.run_medium_preparation_v2 import inventory_module, signature, digest, medium_for, differences


def write(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.out / 'auditor_source.py')
    def forbidden(*a, **kw):
        raise RuntimeError('Optimization forbidden in the structural audit')
    cobra.Model.optimize = forbidden
    cobra.Model.slim_optimize = forbidden
    manifest_path = args.run / 'manifest.json'
    manifest = json.loads(manifest_path.read_text())
    verification = verify_study(ROOT, manifest)
    if not verification['valid']:
        raise ValueError('Changed frozen inputs')
    plan = manifest['plan']['recipe']
    source = inventory_module()
    cases, stress = [], []
    for item in plan['models']:
        model, _ = source.build_parent(item['label'])
        policy = CompartmentPolicy(**item['compartments'])
        conditions = source.metadata_conditions(item['organism'])
        names = sorted({c.media for c in conditions})
        media = [base_medium(name) for name in names]
        legacy = model.copy()
        added = complete_medium_transport(legacy, names, list(policy.completion_exclude))
        for r in legacy.exchanges:
            r.bounds = (0., 1000.)
        for condition in conditions:
            medium = medium_for(condition)
            old = legacy.copy(); apply_medium(old, medium)
            new = prepare_medium(model, medium, policy=policy, completion_media=media, missing_policy='report')
            a, b = signature(old), signature(new.model)
            diff = differences(a, b)
            if diff:
                write(args.out / 'difference.json', {'model': item['label'], 'condition': condition.key, 'differences': diff})
                raise ValueError('Structural comparison differs')
            cases.append({'model': item['label'], 'condition': condition.key, 'legacy_sha256': digest(a),
                          'new_sha256': digest(b), 'new_boundary_report': new.report})
        reference = next((c for c in conditions if 'EX_glc__D_e' in c.exchanges), conditions[0])
        ref = medium_for(reference)
        supplemented = Medium('completion_withdrawal_start', 'Synthetic bound-reset stress test',
                              {**ref.uptakes, **{rid: -.001 for rid in added}})
        withdrawn = Medium('completion_withdrawal_target', 'Remove all newly completed components',
                           {rid: lb for rid, lb in ref.uptakes.items() if rid not in added})
        old_reuse = legacy.copy(); apply_medium(old_reuse, supplemented); apply_medium(old_reuse, withdrawn)
        old_fresh = legacy.copy(); apply_medium(old_fresh, withdrawn)
        new_start = prepare_medium(model, supplemented, policy=policy, completion_media=media, missing_policy='report')
        new_reuse = prepare_medium(new_start.model, withdrawn, policy=policy, completion_media=media, missing_policy='report')
        new_fresh = prepare_medium(model, withdrawn, policy=policy, completion_media=media, missing_policy='report')
        if signature(new_reuse.model) != signature(new_fresh.model):
            raise ValueError('New preparation depends on medium history')
        stress.append({'label': item['label'], 'legacy_retained_uptakes': {
            rid: list(old_reuse.reactions.get_by_id(rid).bounds) for rid in added if old_reuse.reactions.get_by_id(rid).lower_bound < 0},
            'legacy_reuse_vs_fresh': differences(signature(old_fresh), signature(old_reuse)),
            'strict_reuse_signature_sha256': digest(signature(new_reuse.model)),
            'strict_fresh_signature_sha256': digest(signature(new_fresh.model)),
            'strict_reuse_report': new_reuse.report, 'strict_fresh_report': new_fresh.report})
        print('Structural comparison complete:', item['label'], len(conditions), flush=True)
    if len(cases) != plan['case_count_from_source_inventory']:
        raise ValueError('Incomplete structural panel')
    verification = verify_study(ROOT, manifest)
    if not verification['valid']:
        raise ValueError('Inputs changed during structural audit')
    write(args.out / 'report.json', {'n_cases': len(cases), 'all_exactly_equal': True, 'cases': cases,
        'withdrawal_checks': stress, 'input_verification': verification,
        'input_sha256': {str(Path(__file__).relative_to(ROOT)): sha256_of(__file__),
                         str(manifest_path.resolve().relative_to(ROOT)): sha256_of(manifest_path)},
        'scope': 'Post-failure structural completion using the declared comparison. No optimization; prior numerical failures remain unresolved. This is not the independent artifact checker.'})


if __name__ == '__main__':
    main()
