"""Audit saved Sprint 3 outputs without rerunning or replacing their simulations.

Writes corrected point metrics and paired development comparisons to a separate
result directory. Optional metadata correction never changes numerical outputs.
"""
from __future__ import annotations
import argparse
from dataclasses import asdict
import json
import os
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from gembench import metrics as M
from gembench.cards import BenchmarkCard, carbon_fitness_leakage, sha256_of
from gembench.comparison import load_run, compare_runs


def json_safe(value):
    if isinstance(value, dict):
        return {k: json_safe(v) for k, v in value.items()}
    if isinstance(value, list):
        return [json_safe(v) for v in value]
    if isinstance(value, (float, np.floating)) and not np.isfinite(value):
        return None
    return value


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--results-dir', type=Path, default=ROOT / 'results/carbon_fitness_multi')
    ap.add_argument('--out', type=Path, default=ROOT / 'results/audit_2026_09_06')
    ap.add_argument('--correct-leakage', action='store_true')
    ap.add_argument('--n-boot', type=int, default=500)
    args = ap.parse_args()
    args.results_dir = args.results_dir.resolve()
    args.out = args.out.resolve()
    args.out.mkdir(parents=True, exist_ok=True)
    rows = []
    for path in sorted(args.results_dir.glob('*/*/card.json')):
        card = json.loads(path.read_text())
        r = load_run(path.parent)
        st, ft = r['params']['growth_threshold'], r['params']['fitness_threshold']
        grows = np.isfinite(r['wt_growth']) & (r['wt_growth'] >= st)
        s, f = r['sim_growth'][:, grows], r['fitness'][:, grows]
        org = path.parent.parent.name
        patched = bool(card['protocol'].get('patches'))
        old = card['results']['gene_level_conditions_where_wt_grows']
        metrics = {name: fn(s, f) for name, fn in [
            ('mcc', lambda s,f: M.mcc(s,f,st,ft)),
            ('aucpr_bernstein', lambda s,f: M.aucpr_bernstein(s,f,st)),
            ('aucpr_standard', lambda s,f: M.aucpr_standard(s,f,ft)),
            ('auroc_standard', lambda s,f: M.auroc_standard(s,f,ft))]}
        changes = [k for k,v in metrics.items() if k in old and not np.isclose(v, old[k]['point'], atol=1e-9, rtol=0, equal_nan=True)]
        row = dict(run_dir=os.path.relpath(path.parent, ROOT), matrix_sha256=sha256_of(path.with_name('matrices.npz')),
                   n_nonfinite_fitness=int((~np.isfinite(f)).sum()), n_nonfinite_simulation=int((~np.isfinite(s)).sum()),
                   n_scored_pairs=int((np.isfinite(s)&np.isfinite(f)).sum()), metrics=metrics,
                   changed_point_metrics=changes, evaluation_role='retrospective_development')
        if args.correct_leakage:
            audit = card['protocol'].setdefault('reporting_audit', {
                'date': '2026-09-06', 'base_commit': 'b1997d3',
                'note': 'Leakage metadata corrected after review. Historical simulations and numeric results were not rerun or altered.'})
            card['leakage'] = asdict(carbon_fitness_leakage(org, card['protocol']['variant'], patched=patched,
                                      medium_completion=r['params'].get('complete_medium_transport', False)))
            card['protocol']['evaluation_role'] = 'retrospective_development'
            path.write_text(json.dumps(card, indent=2) + '\n')
            path.with_suffix('.md').write_text(BenchmarkCard(**card).to_markdown())
        rows.append(row)
    if not rows:
        raise ValueError('No saved benchmark cards found')
    paired = []
    for org in ('Btheta', 'Putida', 'MR1', 'Smeli'):
        dirs = list((args.results_dir / org).iterdir())
        base = [d for d in dirs if d.name.endswith('__gapfilled__nodroprich')]
        # The last complete arm at handover; Smeli cycle 6 had not been saved.
        suffix = ('v0.1+model-v0.2+gpr-v0.2+0.3+medium__nodroprich' if org == 'Smeli' else
                  'v0.1+model-v0.3+gpr-v0.2+0.3+0.4+medium__nodroprich')
        final = [d for d in dirs if d.name.endswith('__patched-' + suffix)]
        if len(base) != 1 or len(final) != 1:
            raise ValueError(f'Expected one baseline and completed arm for {org}')
        report = compare_runs(load_run(base[0]), load_run(final[0]), n_boot=args.n_boot)
        report.update(org=org, A=os.path.relpath(base[0], ROOT), B=os.path.relpath(final[0], ROOT))
        paired.append(report)
    audit = dict(base_commit='b1997d3', n_historical_cards=len(rows), runs=rows,
                 note='Recalculation from saved matrices cannot recover solver statuses discarded by the original pipeline; this is not a rerun or independent validation.')
    (args.out / 'historical_card_audit.json').write_text(json.dumps(json_safe(audit), indent=2, allow_nan=False) + '\n')
    (args.out / 'paired_development_comparisons.json').write_text(json.dumps(json_safe(paired), indent=2, allow_nan=False) + '\n')
    print(json.dumps({'cards': len(rows), 'runs_with_changed_point_metrics': sum(bool(r['changed_point_metrics']) for r in rows),
                      'paired': [{'org': r['org'], 'mcc': r['mcc'], 'delta': r['paired_mcc_difference_B_minus_A'],
                                  'wt_grows': r['wt_grows'], 'shared_pairs': r['n_shared_finite_pairs']} for r in paired]}, indent=2))


if __name__ == '__main__':
    main()
