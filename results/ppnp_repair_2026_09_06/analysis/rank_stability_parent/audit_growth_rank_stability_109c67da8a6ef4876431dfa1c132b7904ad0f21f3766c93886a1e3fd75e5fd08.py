"""Post-outcome sensitivity of rank metrics to sub-tolerance growth differences.

This never changes a saved prediction, score, primary analysis, or model. The
grid is a numerical diagnostic, not a fitted precision policy or new benchmark.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.verify_quinone_repair import load_benchmark, metric_points, safe

QUANTA = (1e-12, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6)


def transformed(values, quantum, policy):
    result = np.asarray(values, dtype=float).copy()
    if quantum <= 0 or not np.isfinite(quantum):
        raise ValueError('Quantum must be finite and positive')
    finite = np.isfinite(result)
    if policy == 'zero_band':
        result[finite & (np.abs(result) <= quantum)] = 0
    elif policy == 'round_to_quantum':
        result[finite] = np.round(result[finite] / quantum) * quantum
    else:
        raise ValueError('Unknown precision policy')
    return result


def evaluate(run):
    grows = np.isfinite(run['wt_growth']) & (run['wt_growth'] >= run['st'])
    sim, fit = run['sim_growth'][:, grows], run['fitness'][:, grows]
    finite = np.isfinite(sim) & np.isfinite(fit)
    records = [{'policy': 'raw', 'quantum': None,
                'metrics': metric_points(sim, fit, run['st'], run['ft']),
                'changed_binary_predictions': 0}]
    for policy in ('zero_band', 'round_to_quantum'):
        for quantum in QUANTA:
            changed = transformed(sim, quantum, policy)
            records.append({'policy': policy, 'quantum': quantum,
                'metrics': metric_points(changed, fit, run['st'], run['ft']),
                'changed_binary_predictions': int((finite & ((sim >= run['st']) != (changed >= run['st']))).sum())})
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runs', type=Path, nargs='+', required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError('Preserve earlier sensitivity reports')
    args.out.mkdir(parents=True)
    paths = [Path(__file__).resolve(), ROOT / 'scripts/verify_quinone_repair.py']
    paths.extend(p for d in args.runs for p in (d / 'matrices.npz', d / 'card.json'))
    hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    for p in paths[:2]:
        (args.out / f'{p.stem}_{hashes[str(p)]}.py').write_bytes(p.read_bytes())
    records = [{'run': str(d), 'sensitivity': evaluate(load_benchmark(d))} for d in args.runs]
    if hashes != {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}:
        raise ValueError('Inputs changed during diagnostic')
    report = {'role': 'Post-outcome numerical sensitivity prompted by parent rerun; not a replacement primary analysis.',
        'input_sha256': hashes, 'quanta': QUANTA, 'runs': records,
        'interpretation': 'No precision policy is selected using these scores. Small growth perturbations can reorder tied or nearly tied examples; binary threshold decisions are audited separately. Neither rounding nor near-zero clamping establishes biochemical validity.',
        'bootstrap_intervals_recomputed': False}
    with (args.out / 'report.json').open('x') as stream:
        json.dump(safe(report), stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps({'runs': len(records), 'reports_per_run': len(records[0]['sensitivity'])}))


if __name__ == '__main__':
    main()
