"""Metrics for growth/no-growth phenotype prediction, with bootstrap confidence intervals.

Two AUC-PR conventions are provided because they answer different questions:

* `aucpr_bernstein` reproduces Bernstein et al. 2023 exactly: the model's binary
  no-growth call is treated as the *label*, experimental fitness (negated) as the
  *score*, and the PR curve is traced over fitness thresholds. It measures how well
  the experimental fitness ranking recovers the model's no-growth calls without
  choosing a fitness cutoff.
* `aucpr_standard` is the conventional direction: the experimental phenotype
  (fitness < fit_thresh => 'important') is the label and the model's growth ratio
  (negated) is the score. Because FBA growth ratios are close to binary this is
  mostly a precision/recall summary of a binary predictor, which is why MCC and
  balanced accuracy at fixed thresholds are reported alongside.

Bootstrap CIs resample *genes* (rows) with replacement, because the conditions of
one gene are not independent observations.
"""
from __future__ import annotations

from typing import Callable, Dict, Tuple

import numpy as np
from sklearn.metrics import (auc, average_precision_score, matthews_corrcoef,
                             precision_recall_curve, roc_auc_score)


def aucpr_bernstein(sim: np.ndarray, fit: np.ndarray, sim_thresh: float = 1e-3) -> float:
    y = (sim > sim_thresh).astype(int).ravel()       # 1 = predicted growth, 0 = predicted no growth
    score = -np.asarray(fit, dtype=float).ravel()     # more negative fitness -> higher score
    ok = np.isfinite(score)
    if y[ok].min() == y[ok].max():
        return float("nan")
    pre, rec, _ = precision_recall_curve(y[ok], score[ok], pos_label=0)
    return float(auc(rec, pre))


def aucpr_standard(sim: np.ndarray, fit: np.ndarray, fit_thresh: float = -2.0) -> float:
    label = (np.asarray(fit, dtype=float).ravel() < fit_thresh).astype(int)  # 1 = experimentally important
    score = -np.asarray(sim, dtype=float).ravel()                           # lower growth -> higher score
    ok = np.isfinite(label) & np.isfinite(score)
    if label[ok].min() == label[ok].max():
        return float("nan")
    return float(average_precision_score(label[ok], score[ok]))


def auroc_standard(sim: np.ndarray, fit: np.ndarray, fit_thresh: float = -2.0) -> float:
    label = (np.asarray(fit, dtype=float).ravel() < fit_thresh).astype(int)
    score = -np.asarray(sim, dtype=float).ravel()
    ok = np.isfinite(score)
    if label[ok].min() == label[ok].max():
        return float("nan")
    return float(roc_auc_score(label[ok], score[ok]))


def confusion(sim: np.ndarray, fit: np.ndarray, sim_thresh: float = 1e-3, fit_thresh: float = -2.0) -> Dict[str, int]:
    pred_growth = (np.asarray(sim).ravel() > sim_thresh)
    exp_growth = (np.asarray(fit, dtype=float).ravel() >= fit_thresh)
    ok = np.isfinite(np.asarray(fit, dtype=float).ravel())
    pred_growth, exp_growth = pred_growth[ok], exp_growth[ok]
    tp = int(np.sum(pred_growth & exp_growth))      # growth predicted and observed
    tn = int(np.sum(~pred_growth & ~exp_growth))    # no growth predicted and observed
    fp = int(np.sum(pred_growth & ~exp_growth))     # growth predicted, not observed (model too permissive)
    fn = int(np.sum(~pred_growth & exp_growth))     # no growth predicted, growth observed (model too strict)
    return {"tp": tp, "tn": tn, "fp": fp, "fn": fn}


def mcc(sim: np.ndarray, fit: np.ndarray, sim_thresh: float = 1e-3, fit_thresh: float = -2.0) -> float:
    pred = (np.asarray(sim).ravel() > sim_thresh).astype(int)
    lab = (np.asarray(fit, dtype=float).ravel() >= fit_thresh).astype(int)
    ok = np.isfinite(np.asarray(fit, dtype=float).ravel())
    return float(matthews_corrcoef(lab[ok], pred[ok]))


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
    if sim.ndim == 1:
        sim = sim[:, None]; fit = fit[:, None]
    rng = np.random.default_rng(seed)
    point = metric(sim, fit)
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
