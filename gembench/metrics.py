"""Metrics for growth/no-growth phenotype prediction, with bootstrap confidence intervals.

Two AUC-PR conventions are provided because they answer different questions:

* `aucpr_bernstein` reproduces Bernstein et al. 2023 exactly: the model's binary
  no-growth call is treated as the *label*, experimental fitness (negated) as the
  *score*, and the PR curve is traced over fitness thresholds. It measures how well
  the experimental fitness ranking recovers the model's no-growth calls without
  choosing a fitness cutoff.
* `aucpr_standard` is the conventional direction: the experimental phenotype
  (fitness < fit_thresh => 'important') is the label and the model's growth
  prediction (negated) is the score. The generic protocol passes absolute biomass
  flux, not a wild-type-normalized ratio: pooled rankings therefore depend on
  growth-rate scales across media and on numerical ties. MCC and balanced
  accuracy at fixed thresholds are reported alongside.

Bootstrap CIs resample *genes* (rows) with replacement, because the conditions of
one gene are not independent observations.
"""
from __future__ import annotations

from typing import Callable, Dict, Tuple

import numpy as np
from sklearn.metrics import (auc, average_precision_score, matthews_corrcoef,
                             precision_recall_curve, roc_auc_score)


def _finite_pairs(sim: np.ndarray, fit: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Mask observations before thresholding: NaN is an unknown, not a phenotype."""
    sim = np.asarray(sim, dtype=float)
    fit = np.asarray(fit, dtype=float)
    if sim.shape != fit.shape:
        raise ValueError("simulation and fitness arrays must have identical shapes")
    valid = np.isfinite(sim) & np.isfinite(fit)
    return sim[valid], fit[valid]


def aucpr_bernstein(sim: np.ndarray, fit: np.ndarray, sim_thresh: float = 1e-3) -> float:
    sim, fit = _finite_pairs(sim, fit)
    y = (sim >= sim_thresh).astype(int)  # 1 = growth; threshold matches the protocol
    if y.size == 0 or y.min() == y.max():
        return float("nan")
    pre, rec, _ = precision_recall_curve(y, -fit, pos_label=0)
    return float(auc(rec, pre))


def aucpr_standard(sim: np.ndarray, fit: np.ndarray, fit_thresh: float = -2.0) -> float:
    sim, fit = _finite_pairs(sim, fit)
    label = (fit < fit_thresh).astype(int)
    if label.size == 0 or label.min() == label.max():
        return float("nan")
    return float(average_precision_score(label, -sim))


def auroc_standard(sim: np.ndarray, fit: np.ndarray, fit_thresh: float = -2.0) -> float:
    sim, fit = _finite_pairs(sim, fit)
    label = (fit < fit_thresh).astype(int)
    if label.size == 0 or label.min() == label.max():
        return float("nan")
    return float(roc_auc_score(label, -sim))


def confusion(sim: np.ndarray, fit: np.ndarray, sim_thresh: float = 1e-3, fit_thresh: float = -2.0) -> Dict[str, int]:
    sim, fit = _finite_pairs(sim, fit)
    pred_growth = sim >= sim_thresh
    exp_growth = fit >= fit_thresh
    tp = int(np.sum(pred_growth & exp_growth))      # growth predicted and observed
    tn = int(np.sum(~pred_growth & ~exp_growth))    # no growth predicted and observed
    fp = int(np.sum(pred_growth & ~exp_growth))     # growth predicted, not observed (model too permissive)
    fn = int(np.sum(~pred_growth & exp_growth))     # no growth predicted, growth observed (model too strict)
    return {"tp": tp, "tn": tn, "fp": fp, "fn": fn}


def mcc(sim: np.ndarray, fit: np.ndarray, sim_thresh: float = 1e-3, fit_thresh: float = -2.0) -> float:
    sim, fit = _finite_pairs(sim, fit)
    if sim.size == 0:
        return float("nan")
    return float(matthews_corrcoef(fit >= fit_thresh, sim >= sim_thresh))


def balanced_accuracy(sim: np.ndarray, fit: np.ndarray, sim_thresh: float = 1e-3, fit_thresh: float = -2.0) -> float:
    c = confusion(sim, fit, sim_thresh, fit_thresh)
    tpr = c["tp"] / (c["tp"] + c["fn"]) if (c["tp"] + c["fn"]) else float("nan")
    tnr = c["tn"] / (c["tn"] + c["fp"]) if (c["tn"] + c["fp"]) else float("nan")
    return float((tpr + tnr) / 2)


def accuracy(sim: np.ndarray, fit: np.ndarray, sim_thresh: float = 1e-3, fit_thresh: float = -2.0) -> float:
    c = confusion(sim, fit, sim_thresh, fit_thresh)
    n = sum(c.values())
    return float((c["tp"] + c["tn"]) / n) if n else float("nan")


def bootstrap_ci(metric: Callable[[np.ndarray, np.ndarray], float], sim: np.ndarray, fit: np.ndarray,
                 n_boot: int = 1000, seed: int = 0, alpha: float = 0.05) -> Tuple[float, float, float]:
    """Percentile bootstrap over genes (rows). Returns (point, lower, upper)."""
    sim = np.asarray(sim); fit = np.asarray(fit, dtype=float)
    if sim.shape != fit.shape or sim.ndim not in (1, 2):
        raise ValueError("bootstrap expects matching one- or two-dimensional arrays")
    if sim.ndim == 1:
        sim = sim[:, None]; fit = fit[:, None]
    rng = np.random.default_rng(seed)
    point = metric(sim, fit)
    if sim.shape[0] == 0:
        return point, float("nan"), float("nan")
    vals = []
    n = sim.shape[0]
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        try:
            vals.append(metric(sim[idx], fit[idx]))
        except Exception:
            vals.append(float("nan"))
    vals = np.array(vals, dtype=float)
    vals = vals[np.isfinite(vals)]
    if len(vals) == 0:
        return point, float("nan"), float("nan")
    return point, float(np.percentile(vals, 100 * alpha / 2)), float(np.percentile(vals, 100 * (1 - alpha / 2)))
