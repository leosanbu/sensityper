#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
sensityping_metrics_v2.0.py

Metrics for antibiotic predictions and first-line recommendations
(supports provided combo columns like CRO+AZM).

================================================================================
VERSION HISTORY / CHANGELOG
================================================================================

v2.0 (2026-01-08)
  - Added 95%% confidence intervals (CIs) where applicable:
      * Proportion-type metrics: Wilson score interval
        (accuracy, concordance, coverage_fraction, sensitivity, specificity, PPV, NPV, FDR, ME_rate, VME_rate)
      * Composite metrics: bootstrap percentile CIs (enabled via --ci_method bootstrap/hybrid)
        (f1_score, mcc, kappa, balanced_accuracy, auc, cost_sensitive_error_rate)
  - Added SSD (SampleSizeDiagnostics-style) sample size requirement:
      * Computes required n for a target CI half-width w (approx normal):
          n >= z^2 * p(1-p) / w^2
      * Modes:
          - conservative: p=0.5 (max variance; most conservative)
          - observed:     p = observed proportion (p-hat)
          - fixed:        p = user-provided --ssd_p
      * Reports <metric>_ssd_n_required and <metric>_ssd_met (True/False), plus denominators.
  - Added version metadata and CLI flag:
      * __version__ string
      * --version prints version and exits
  - Python 3.6 compatibility fixes:
      * Avoids PEP604 union types (uses typing.Optional)
      * Escapes percent signs (%%) in argparse help/epilog to avoid argparse formatting crash
      * Uses np.random.RandomState (avoids numpy default_rng dependency)
  - Output quality improvements:
      * Suppresses sklearn RuntimeWarnings locally (no noisy stdout)
      * Adds a short per-step 'note' field flagging single-class/degenerate subsets

v1.x (prior)
  - Baseline metrics computation for:
      * predicted_vs_treatment  : <ABX>_predicted vs <ABX>_treatment
      * first_line_vs_treatment : sequential <TOKEN>_recommend vs <TOKEN>_treatment using --order
  - Radar chart output (Plotly)
  - ID extraction tables for TP/TN/FP/FN

================================================================================
EXAMPLE COMMANDS
================================================================================

1) Per-antibiotic predictions vs treatment (no combos here) + CIs + SSD:
   python sensityping_metrics_v2.0.py -i results.tsv -o out.txt -d ./out \
      --analysis_type predicted_vs_treatment \
      --ci_flag --ci_method hybrid --ci_level 0.95 --n_boot 2000 --seed 1 \
      --ssd_flag --ssd_width 0.05 --ssd_mode conservative \
      --radar_flag --radar_metrics PPV,one_minus_FDR,coverage_fraction

2) First-line (rule-based) vs treatment with order including a combo token present in the table:
   python sensityping_metrics_v2.0.py -i treatment_output.tsv -o firstline_metrics.txt -d ./out \
      --analysis_type first_line_vs_treatment \
      --order CIP,CRO+AZM,CRO,SPC \
      --ci_flag --ci_method hybrid \
      --ssd_flag --ssd_width 0.05 \
      --radar_flag --radar_metrics PPV,one_minus_FDR,coverage_fraction \
      --id_extraction CIP,CRO+AZM

================================================================================
NOTES
================================================================================
- 'Yes' parsing is strict: only literal 'Yes' (after str().strip()) counts as Yes.
- For first_line_vs_treatment:
    After evaluating a token, isolates predicted 'Yes' for that token are removed before the next step.
- TOKEN can be a combo (e.g., CRO+AZM) if such columns exist in the input table.
- CI:
    Proportion metrics use Wilson CIs.
    Composite metrics use bootstrap CIs if --ci_method is bootstrap/hybrid.
- SSD:
    --ssd_width is the CI half-width w (0.05 => ±0.05; total width 0.10).
    --ssd_mode conservative uses p=0.5; observed uses p-hat; fixed uses --ssd_p.
"""

__version__ = "2.0"

import os
import argparse
import math
import warnings
from typing import Optional

import numpy as np
import pandas as pd
from sklearn.metrics import (
    confusion_matrix,
    f1_score,
    matthews_corrcoef,
    roc_auc_score,
    cohen_kappa_score
)
import plotly.graph_objects as go

# ---------------------------
# Utilities
# ---------------------------

def normalize_yes_no(seq):
    """Return list with exact 'Yes' or 'No' strings. Strict: only literal 'Yes' counts as Yes."""
    return ['Yes' if str(x).strip() == 'Yes' else 'No' for x in seq]

def classify_ids(id_series, true_labels, pred_labels):
    """
    Given parallel sequences of ids, y_true ('Yes'/'No'), y_pred ('Yes'/'No'),
    return dict with lists of IDs per confusion quadrant: TP, TN, FP, FN.
    """
    ids = list(id_series)
    y_true = normalize_yes_no(true_labels)
    y_pred = normalize_yes_no(pred_labels)

    TP, TN, FP, FN = [], [], [], []
    for i in range(len(ids)):
        t, p = y_true[i], y_pred[i]
        if t == 'Yes' and p == 'Yes':
            TP.append(ids[i])
        elif t == 'Yes' and p == 'No':
            FN.append(ids[i])
        elif t == 'No' and p == 'Yes':
            FP.append(ids[i])
        else:  # t == 'No' and p == 'No'
            TN.append(ids[i])
    return {'TP': TP, 'TN': TN, 'FP': FP, 'FN': FN}

def write_id_extraction_table(path, id_lists):
    """
    Write a .tab file with columns TP, TN, FP, FN.
    Each row aligns the kth element of each list (padded with blank if shorter).
    """
    out_dir = os.path.dirname(path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    max_len = max((len(id_lists[k]) for k in ['TP', 'TN', 'FP', 'FN']), default=0)
    with open(path, 'w') as f:
        f.write("TP\tTN\tFP\tFN\n")
        for i in range(max_len):
            row = [
                id_lists['TP'][i] if i < len(id_lists['TP']) else '',
                id_lists['TN'][i] if i < len(id_lists['TN']) else '',
                id_lists['FP'][i] if i < len(id_lists['FP']) else '',
                id_lists['FN'][i] if i < len(id_lists['FN']) else '',
            ]
            f.write('\t'.join(row) + '\n')

def parse_order_token(token):
    """Return token string as-is (no AND logic)."""
    return token.strip()

def needed_columns_for_token(token):
    """For token, require <token>_treatment and <token>_recommend."""
    return ["%s_treatment" % token, "%s_recommend" % token]

# ---------------------------
# CI + SSD helpers
# ---------------------------

def z_from_confidence(conf_level):
    """
    Convert confidence level to z-score using an approximation to inverse normal CDF.
    Avoids external deps (scipy).
    """
    # For 95%% CI, return exact common constant.
    if abs(conf_level - 0.95) < 1e-12:
        return 1.959963984540054

    p = 0.5 + conf_level / 2.0
    if p <= 0.0 or p >= 1.0:
        raise ValueError("conf_level must be between 0 and 1 (exclusive).")

    # Acklam inverse normal CDF approximation
    a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
         1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00]
    b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
         6.680131188771972e+01, -1.328068155288572e+01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
         -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00]
    d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
         3.754408661907416e+00]

    plow = 0.02425
    phigh = 1 - plow

    if p < plow:
        q = math.sqrt(-2 * math.log(p))
        num = (((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5])
        den = ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1)
        return num / den
    if p > phigh:
        q = math.sqrt(-2 * math.log(1 - p))
        num = -(((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5])
        den = ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1)
        return num / den

    q = p - 0.5
    r = q * q
    num = (((((a[0]*r + a[1])*r + a[2])*r + a[3])*r + a[4])*r + a[5]) * q
    den = (((((b[0]*r + b[1])*r + b[2])*r + b[3])*r + b[4])*r + 1)
    return num / den

def wilson_ci(k, n, conf_level=0.95):
    """Wilson score interval for a binomial proportion."""
    if n <= 0:
        return None
    z = z_from_confidence(conf_level)
    p = float(k) / float(n)
    denom = 1.0 + (z*z)/float(n)
    center = (p + (z*z)/(2.0*float(n))) / denom
    half = (z / denom) * math.sqrt((p*(1.0-p)/float(n)) + (z*z)/(4.0*float(n)*float(n)))
    lo = max(0.0, center - half)
    hi = min(1.0, center + half)
    return (lo, hi)

def clamp01(x):
    return max(0.0, min(1.0, x))

def bootstrap_ci(y_true_bin,
                 y_pred_bin,
                 metric_fn,
                 conf_level=0.95,
                 n_boot=2000,
                 seed=1,
                 min_valid_fraction=0.8):
    """
    Nonparametric bootstrap CI for metrics.
    - y_true_bin, y_pred_bin are 0/1 arrays (1 = 'Yes').
    - metric_fn returns float or raises.
    - If too many undefined replicates, return None.
    """
    n = len(y_true_bin)
    if n == 0:
        return None

    rng = np.random.RandomState(seed)
    vals = []
    for _ in range(int(n_boot)):
        idx = rng.randint(0, n, size=n)
        try:
            v = metric_fn(y_true_bin[idx], y_pred_bin[idx])
            if v is None:
                continue
            v = float(v)
            if math.isnan(v) or math.isinf(v):
                continue
            vals.append(v)
        except Exception:
            continue

    if len(vals) < max(1, int(min_valid_fraction * n_boot)):
        return None

    vals = np.array(vals, dtype=float)
    alpha = 1.0 - conf_level
    lo = float(np.quantile(vals, alpha/2.0))
    hi = float(np.quantile(vals, 1.0 - alpha/2.0))
    return (lo, hi)

def ssd_required_n(half_width, conf_level=0.95, p=0.5):
    """
    Sample size required for a binomial proportion to achieve a target CI half-width (approx, normal):
      n >= z^2 * p(1-p) / w^2
    """
    if half_width <= 0:
        raise ValueError("SSD half-width (w) must be > 0.")
    z = z_from_confidence(conf_level)
    p = clamp01(float(p))
    n = (z*z) * p * (1.0 - p) / (half_width * half_width)
    return int(math.ceil(n))

def add_prop_ci_and_ssd(out,
                        name,
                        k,
                        n,
                        conf_level,
                        ssd_flag,
                        ssd_half_width,
                        ssd_mode,
                        ssd_p_override):
    """
    Adds:
      - <name>_ci_low / <name>_ci_high (Wilson)
      - <name>_denom (n)
      - <name>_ssd_n_required and <name>_ssd_met (if ssd_flag)
    """
    if n is None or int(n) <= 0:
        out["%s_ci_low" % name] = 'N/A'
        out["%s_ci_high" % name] = 'N/A'
        out["%s_denom" % name] = 'N/A'
        if ssd_flag:
            out["%s_ssd_n_required" % name] = 'N/A'
            out["%s_ssd_met" % name] = 'N/A'
        return

    ci = wilson_ci(int(k), int(n), conf_level=conf_level)
    out["%s_ci_low" % name] = round(ci[0], 5) if ci else 'N/A'
    out["%s_ci_high" % name] = round(ci[1], 5) if ci else 'N/A'
    out["%s_denom" % name] = int(n)

    if ssd_flag:
        if ssd_mode == 'conservative':
            p_for_ssd = 0.5
        elif ssd_mode == 'observed':
            p_for_ssd = (float(k) / float(n)) if int(n) > 0 else 0.5
        elif ssd_mode == 'fixed':
            if ssd_p_override is None:
                raise ValueError("--ssd_mode fixed requires --ssd_p")
            p_for_ssd = float(ssd_p_override)
        else:
            raise ValueError("Unknown ssd_mode: %s" % ssd_mode)

        n_req = ssd_required_n(ssd_half_width, conf_level=conf_level, p=p_for_ssd)
        out["%s_ssd_n_required" % name] = int(n_req)
        out["%s_ssd_met" % name] = bool(int(n) >= int(n_req))

# ---------------------------
# Metrics
# ---------------------------

def calculate_metrics(true_labels, predicted_labels,
                      fp_weight=5, fn_weight=1, yes_means='susceptible',
                      ci_flag=False, ci_method='hybrid', ci_level=0.95,
                      n_boot=2000, seed=1,
                      ssd_flag=False, ssd_half_width=0.05, ssd_mode='conservative', ssd_p=None):
    """
    yes_means: 'susceptible' or 'resistant' (default: 'susceptible')

    CI strategy:
      - Proportion metrics: Wilson
      - Composite metrics: bootstrap if ci_method in {'bootstrap','hybrid'}
    """
    true_labels = normalize_yes_no(true_labels)
    predicted_labels = normalize_yes_no(predicted_labels)

    # Better fix: silence sklearn RuntimeWarnings locally (division-by-zero etc.)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        cm = confusion_matrix(true_labels, predicted_labels, labels=['Yes', 'No'])
        tp = int(cm[0, 0])
        fn = int(cm[0, 1])
        fp = int(cm[1, 0])
        tn = int(cm[1, 1])

    # Short note for degenerate / single-class subsets
    note = None
    n_true_classes = len(set(true_labels))
    n_pred_classes = len(set(predicted_labels))

    if n_true_classes == 1 and n_pred_classes == 1:
        only = true_labels[0] if len(true_labels) else "N/A"
        note = "Single-class subset: all isolates are '%s'; sensitivity/PPV/NPV/AUC may be not estimable." % only
    elif n_true_classes == 1:
        note = "Single-class truth: some rate metrics (e.g., sensitivity/specificity) not estimable."
    elif n_pred_classes == 1:
        note = "Single-class prediction: some predictive value metrics (PPV/NPV) not estimable."
    else:
        note = "OK"

    # Print counts for sanity-check
    print("TP: %d, TN: %d, FP: %d, FN: %d" % (tp, tn, fp, fn))

    total = tp + tn + fp + fn
    assigned = tp + fp  # predicted 'Yes'

    accuracy = (tp + tn) / float(total) if total else 0.0
    sensitivity = tp / float(tp + fn) if (tp + fn) else 0.0
    specificity = tn / float(tn + fp) if (tn + fp) else 0.0
    ppv = tp / float(tp + fp) if (tp + fp) else 0.0
    npv = tn / float(tn + fn) if (tn + fn) else 0.0

    fdr = fp / float(tp + fp) if (tp + fp) else 'N/A'
    one_minus_fdr = 1.0 - fdr if isinstance(fdr, float) else 'N/A'
    coverage_fraction = assigned / float(total) if total else 0.0

    # AUC: needs both classes in truth
    auc = 'N/A'
    if len(set(true_labels)) > 1:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            try:
                auc = roc_auc_score(
                    [1 if x == 'Yes' else 0 for x in true_labels],
                    [1 if x == 'Yes' else 0 for x in predicted_labels]
                )
            except Exception:
                auc = 'N/A'

    balanced_accuracy = (sensitivity + specificity) / 2.0 if ((tp + fn) and (tn + fp)) else 'N/A'
    cost_sensitive_error_rate = (fp * fp_weight + fn * fn_weight) / float(total) if total else 0.0

    # Composite metrics (sklearn)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        if len(set(true_labels)) > 1 and len(set(predicted_labels)) > 1:
            try:
                f1 = f1_score(true_labels, predicted_labels, pos_label='Yes', average='binary', zero_division=0)
            except Exception:
                f1 = 0.0
            try:
                mcc = matthews_corrcoef(true_labels, predicted_labels)
            except Exception:
                mcc = 0.0
            try:
                kappa = cohen_kappa_score(true_labels, predicted_labels, labels=['Yes', 'No'])
            except Exception:
                kappa = 0.0
        else:
            f1 = 0.0
            mcc = 0.0
            kappa = 0.0

    me_rate = 'N/A'
    vme_rate = 'N/A'
    me_denom = None
    vme_denom = None
    me_num = None
    vme_num = None

    # Preserve your v1 logic for ME/VME semantics driven by yes_means
    if yes_means == 'resistant':
        me_denom = (tn + fp)   # true susceptible
        vme_denom = (tp + fn)  # true resistant
        me_num = fp
        vme_num = fn
        if me_denom > 0:
            me_rate = me_num / float(me_denom)
        if vme_denom > 0:
            vme_rate = vme_num / float(vme_denom)
    else:
        me_denom = (tp + fn)   # true susceptible (Yes=susceptible)
        vme_denom = (tn + fp)  # true resistant
        me_num = fn
        vme_num = fp
        if me_denom > 0:
            me_rate = me_num / float(me_denom)
        if vme_denom > 0:
            vme_rate = vme_num / float(vme_denom)

    one_minus_me = 1.0 - me_rate if isinstance(me_rate, float) else 'N/A'
    one_minus_vme = 1.0 - vme_rate if isinstance(vme_rate, float) else 'N/A'

    out = {
        'note': note,

        'tp': tp, 'tn': tn, 'fp': fp, 'fn': fn,
        'assigned_count': assigned,

        'accuracy': round(accuracy, 5),
        'concordance': round(accuracy, 5),
        'sensitivity': round(sensitivity, 5),
        'specificity': round(specificity, 5) if (tn + fp) else 'N/A',
        'PPV': round(ppv, 5) if (tp + fp) else 'N/A',
        'NPV': round(npv, 5) if (tn + fn) else 'N/A',
        'f1_score': round(float(f1), 5),
        'mcc': round(float(mcc), 5),
        'auc': round(float(auc), 5) if auc != 'N/A' else 'N/A',
        'balanced_accuracy': round(float(balanced_accuracy), 5) if balanced_accuracy != 'N/A' else 'N/A',
        'cost_sensitive_error_rate': round(float(cost_sensitive_error_rate), 5),
        'kappa': round(float(kappa), 5),

        'ME_rate': round(me_rate, 5) if me_rate != 'N/A' else 'N/A',
        'VME_rate': round(vme_rate, 5) if vme_rate != 'N/A' else 'N/A',
        'ME_denom': me_denom if me_denom is not None else 'N/A',
        'VME_denom': vme_denom if vme_denom is not None else 'N/A',
        'one_minus_ME': round(one_minus_me, 5) if one_minus_me != 'N/A' else 'N/A',
        'one_minus_VME': round(one_minus_vme, 5) if one_minus_vme != 'N/A' else 'N/A',

        'FDR': round(fdr, 5) if fdr != 'N/A' else 'N/A',
        'one_minus_FDR': round(one_minus_fdr, 5) if one_minus_fdr != 'N/A' else 'N/A',
        'coverage_fraction': round(float(coverage_fraction), 5),
    }

    # ---- Wilson CIs + SSD for proportion metrics ----
    if ci_flag and ci_method in ('wilson', 'hybrid', 'bootstrap'):
        add_prop_ci_and_ssd(out, "accuracy", tp + tn, total, ci_level, ssd_flag, ssd_half_width, ssd_mode, ssd_p)
        add_prop_ci_and_ssd(out, "concordance", tp + tn, total, ci_level, ssd_flag, ssd_half_width, ssd_mode, ssd_p)
        add_prop_ci_and_ssd(out, "coverage_fraction", assigned, total, ci_level, ssd_flag, ssd_half_width, ssd_mode, ssd_p)
        add_prop_ci_and_ssd(out, "sensitivity", tp, tp + fn, ci_level, ssd_flag, ssd_half_width, ssd_mode, ssd_p)
        add_prop_ci_and_ssd(out, "specificity", tn, tn + fp, ci_level, ssd_flag, ssd_half_width, ssd_mode, ssd_p)
        add_prop_ci_and_ssd(out, "PPV", tp, tp + fp, ci_level, ssd_flag, ssd_half_width, ssd_mode, ssd_p)
        add_prop_ci_and_ssd(out, "NPV", tn, tn + fn, ci_level, ssd_flag, ssd_half_width, ssd_mode, ssd_p)

        # FDR and transform
        if isinstance(fdr, float):
            add_prop_ci_and_ssd(out, "FDR", fp, tp + fp, ci_level, ssd_flag, ssd_half_width, ssd_mode, ssd_p)
            lo = out.get("FDR_ci_low")
            hi = out.get("FDR_ci_high")
            if isinstance(lo, float) and isinstance(hi, float):
                out["one_minus_FDR_ci_low"] = round(1.0 - hi, 5)
                out["one_minus_FDR_ci_high"] = round(1.0 - lo, 5)
            else:
                out["one_minus_FDR_ci_low"] = 'N/A'
                out["one_minus_FDR_ci_high"] = 'N/A'
        else:
            out["FDR_ci_low"] = 'N/A'
            out["FDR_ci_high"] = 'N/A'
            out["FDR_denom"] = 'N/A'
            out["one_minus_FDR_ci_low"] = 'N/A'
            out["one_minus_FDR_ci_high"] = 'N/A'
            if ssd_flag:
                out["FDR_ssd_n_required"] = 'N/A'
                out["FDR_ssd_met"] = 'N/A'

        # ME and transform
        if isinstance(me_rate, float):
            add_prop_ci_and_ssd(out, "ME_rate", me_num, me_denom, ci_level, ssd_flag, ssd_half_width, ssd_mode, ssd_p)
            lo = out.get("ME_rate_ci_low")
            hi = out.get("ME_rate_ci_high")
            if isinstance(lo, float) and isinstance(hi, float):
                out["one_minus_ME_ci_low"] = round(1.0 - hi, 5)
                out["one_minus_ME_ci_high"] = round(1.0 - lo, 5)
            else:
                out["one_minus_ME_ci_low"] = 'N/A'
                out["one_minus_ME_ci_high"] = 'N/A'
        else:
            out["ME_rate_ci_low"] = 'N/A'
            out["ME_rate_ci_high"] = 'N/A'
            out["ME_rate_denom"] = 'N/A'
            out["one_minus_ME_ci_low"] = 'N/A'
            out["one_minus_ME_ci_high"] = 'N/A'
            if ssd_flag:
                out["ME_rate_ssd_n_required"] = 'N/A'
                out["ME_rate_ssd_met"] = 'N/A'

        # VME and transform
        if isinstance(vme_rate, float):
            add_prop_ci_and_ssd(out, "VME_rate", vme_num, vme_denom, ci_level, ssd_flag, ssd_half_width, ssd_mode, ssd_p)
            lo = out.get("VME_rate_ci_low")
            hi = out.get("VME_rate_ci_high")
            if isinstance(lo, float) and isinstance(hi, float):
                out["one_minus_VME_ci_low"] = round(1.0 - hi, 5)
                out["one_minus_VME_ci_high"] = round(1.0 - lo, 5)
            else:
                out["one_minus_VME_ci_low"] = 'N/A'
                out["one_minus_VME_ci_high"] = 'N/A'
        else:
            out["VME_rate_ci_low"] = 'N/A'
            out["VME_rate_ci_high"] = 'N/A'
            out["VME_rate_denom"] = 'N/A'
            out["one_minus_VME_ci_low"] = 'N/A'
            out["one_minus_VME_ci_high"] = 'N/A'
            if ssd_flag:
                out["VME_rate_ssd_n_required"] = 'N/A'
                out["VME_rate_ssd_met"] = 'N/A'

    # ---- Bootstrap CIs for composite metrics ----
    if ci_flag and ci_method in ('bootstrap', 'hybrid'):
        y_true_bin = np.array([1 if x == 'Yes' else 0 for x in true_labels], dtype=int)
        y_pred_bin = np.array([1 if x == 'Yes' else 0 for x in predicted_labels], dtype=int)

        def _f1(t, p):
            # handle single-class within replicate safely
            if len(np.unique(t)) < 2 and len(np.unique(p)) < 2:
                return None
            return f1_score(t, p, pos_label=1, average='binary', zero_division=0)

        def _mcc(t, p):
            if len(np.unique(t)) < 2 or len(np.unique(p)) < 2:
                return None
            return float(matthews_corrcoef(t, p))

        def _kappa(t, p):
            if len(np.unique(t)) < 2 or len(np.unique(p)) < 2:
                return None
            return float(cohen_kappa_score(t, p, labels=[1, 0]))

        def _ba(t, p):
            cm2 = confusion_matrix(t, p, labels=[1, 0])
            tp2 = cm2[0, 0]
            fn2 = cm2[0, 1]
            fp2 = cm2[1, 0]
            tn2 = cm2[1, 1]
            if (tp2 + fn2) == 0 or (tn2 + fp2) == 0:
                return None
            sens2 = tp2 / float(tp2 + fn2)
            spec2 = tn2 / float(tn2 + fp2)
            return (sens2 + spec2) / 2.0

        def _auc(t, p):
            if len(np.unique(t)) < 2:
                return None
            return float(roc_auc_score(t, p))

        def _cserr(t, p):
            cm2 = confusion_matrix(t, p, labels=[1, 0])
            tp2 = cm2[0, 0]
            fn2 = cm2[0, 1]
            fp2 = cm2[1, 0]
            tn2 = cm2[1, 1]
            tot2 = tp2 + tn2 + fp2 + fn2
            if tot2 == 0:
                return None
            return float((fp2 * fp_weight + fn2 * fn_weight) / float(tot2))

        composite = [
            ("f1_score", _f1),
            ("mcc", _mcc),
            ("kappa", _kappa),
            ("balanced_accuracy", _ba),
            ("auc", _auc),
            ("cost_sensitive_error_rate", _cserr),
        ]
        for name, fn in composite:
            ci = bootstrap_ci(y_true_bin, y_pred_bin, fn, conf_level=ci_level, n_boot=n_boot, seed=seed)
            out["%s_ci_low" % name] = round(ci[0], 5) if ci else 'N/A'
            out["%s_ci_high" % name] = round(ci[1], 5) if ci else 'N/A'

    return out

# ---------------------------
# Plotting
# ---------------------------

def create_combined_radar_chart(tokens, metrics_dict, output_dir, selected_metrics, suffix=''):
    fig = go.Figure()
    for token in tokens:
        if token in metrics_dict:
            vals = []
            for metric in selected_metrics:
                val = metrics_dict[token].get(metric, 'N/A')
                if val == 'N/A':
                    val = 0
                vals.append(val)
            fig.add_trace(go.Scatterpolar(
                r=vals + [vals[0]],
                theta=selected_metrics + [selected_metrics[0]],
                fill='toself',
                name=str(token)
            ))
    fig.update_layout(
        polar=dict(radialaxis=dict(visible=True, range=[0, 1])),
        showlegend=True,
        title="Radar Chart for %s%s" % (", ".join([str(t) for t in tokens]), suffix)
    )
    os.makedirs(output_dir, exist_ok=True)
    fig.write_html(os.path.join(output_dir, "combined_radar%s.html" % suffix))

# ---------------------------
# Main analysis flows
# ---------------------------

def main(input_file, output_file, output_dir, analysis_type, order,
         radar_flag, radar_antibiotics, radar_metrics,
         fp_weight, fn_weight, yes_means,
         id_extraction, id_column,
         ci_flag, ci_method, ci_level, n_boot, seed,
         ssd_flag, ssd_width, ssd_mode, ssd_p):

    os.makedirs(output_dir, exist_ok=True)
    df = pd.read_csv(input_file, sep='\t')

    # Determine ID column
    if id_column and id_column in df.columns:
        id_col = id_column
    else:
        id_col = df.columns[0]

    tokens_order = order.split(',') if order else None
    selected_antibiotics = [x.strip() for x in radar_antibiotics.split(',') if x.strip()]
    selected_metrics = [x.strip() for x in radar_metrics.split(',') if x.strip()]

    if id_extraction:
        requested_tokens = [t.strip() for t in id_extraction.split(',') if t.strip()]
    else:
        requested_tokens = []

    metrics_dict = {}

    with open(output_file, 'w') as f_out:
        f_out.write("# sensityping_metrics_v2.0\n")
        f_out.write("# version: %s\n" % __version__)
        if ci_flag:
            f_out.write("# CI enabled: method=%s, level=%s, n_boot=%s, seed=%s\n" % (ci_method, ci_level, n_boot, seed))
        if ssd_flag:
            f_out.write("# SSD enabled: half-width(w)=%s, mode=%s, p=%s, level=%s\n" % (ssd_width, ssd_mode, str(ssd_p), ci_level))

        if analysis_type == 'predicted_vs_treatment':
            # Per-antibiotic: <ABX>_predicted vs <ABX>_treatment
            for abx in selected_antibiotics:
                t_col = "%s_treatment" % abx
                p_col = "%s_predicted" % abx
                if t_col in df.columns and p_col in df.columns:
                    filtered_df = df.dropna(subset=[t_col, p_col])
                    if filtered_df.empty:
                        continue

                    y_true = filtered_df[t_col].values
                    y_pred = filtered_df[p_col].values

                    metrics = calculate_metrics(
                        y_true, y_pred,
                        fp_weight=fp_weight,
                        fn_weight=fn_weight,
                        yes_means=yes_means,
                        ci_flag=ci_flag,
                        ci_method=ci_method,
                        ci_level=ci_level,
                        n_boot=n_boot,
                        seed=seed,
                        ssd_flag=ssd_flag,
                        ssd_half_width=ssd_width,
                        ssd_mode=ssd_mode,
                        ssd_p=ssd_p
                    )
                    metrics_dict[abx] = metrics

                    f_out.write("\nMetrics for %s (Prediction vs. Treatment):\n" % abx)
                    f_out.write("Total Isolates Analyzed: %d\n" % len(filtered_df))
                    for k, v in metrics.items():
                        f_out.write("%s: %s\n" % (k, v))

                    if abx in requested_tokens:
                        id_lists = classify_ids(filtered_df[id_col], y_true, y_pred)
                        out_path = os.path.join(output_dir, "%s_id_extracted.tab" % abx)
                        write_id_extraction_table(out_path, id_lists)
                        print("[ID extraction] Saved: %s" % out_path)

        elif analysis_type == 'first_line_vs_treatment':
            if tokens_order is None:
                raise ValueError("The --order flag must be specified for first_line_vs_treatment.")
            remaining_df = df.copy()

            for raw_token in tokens_order:
                token = parse_order_token(raw_token)
                needed_cols = needed_columns_for_token(token)

                if not all(c in remaining_df.columns for c in needed_cols):
                    continue

                step_df = remaining_df.dropna(subset=needed_cols)
                if step_df.empty:
                    continue

                y_true = step_df["%s_treatment" % token].values
                y_pred = step_df["%s_recommend" % token].values

                metrics = calculate_metrics(
                    y_true, y_pred,
                    fp_weight=fp_weight,
                    fn_weight=fn_weight,
                    yes_means=yes_means,
                    ci_flag=ci_flag,
                    ci_method=ci_method,
                    ci_level=ci_level,
                    n_boot=n_boot,
                    seed=seed,
                    ssd_flag=ssd_flag,
                    ssd_half_width=ssd_width,
                    ssd_mode=ssd_mode,
                    ssd_p=ssd_p
                )
                metrics_dict[token] = metrics

                f_out.write("\nMetrics for %s (First Line vs. Treatment):\n" % token)
                f_out.write("Total Isolates Analyzed: %d\n" % len(step_df))
                for k, v in metrics.items():
                    f_out.write("%s: %s\n" % (k, v))

                if token in requested_tokens:
                    id_lists = classify_ids(step_df[id_col], y_true, y_pred)
                    safe_token = token.replace('/', '_').replace('\\', '_')
                    out_path = os.path.join(output_dir, "%s_id_extracted.tab" % safe_token)
                    write_id_extraction_table(out_path, id_lists)
                    print("[ID extraction] Saved: %s" % out_path)

                # Remove predicted Yes before next step
                mask_pred_yes = pd.Series([x == 'Yes' for x in normalize_yes_no(y_pred)], index=step_df.index)
                remaining_df = remaining_df.drop(index=step_df.index[mask_pred_yes])

            leftovers = remaining_df.copy()
            if not leftovers.empty:
                f_out.write("\nUnassigned after order: %d isolates were not predicted 'Yes' for any token (or had missing columns).\n" % len(leftovers))
                unassigned_path = os.path.join(output_dir, "UNASSIGNED_after_order_ids.tab")
                leftovers[[id_col]].to_csv(unassigned_path, sep='\t', index=False)
                print("[Unassigned] Saved IDs: %s" % unassigned_path)

        else:
            raise ValueError("Unsupported analysis_type. Use 'predicted_vs_treatment' or 'first_line_vs_treatment'.")

        if radar_flag:
            radar_tokens = tokens_order if tokens_order else selected_antibiotics
            create_combined_radar_chart(
                radar_tokens,
                metrics_dict,
                output_dir,
                selected_metrics,
                suffix='_%s' % analysis_type
            )

# ---------------------------
# CLI
# ---------------------------

if __name__ == "__main__":

    # IMPORTANT: argparse uses percent-formatting internally on help strings in Python 3.6.
    # Any literal '%' must be escaped as '%%' to avoid crashes when rendering -h.

    epilog = """
Examples
--------
1) Per-antibiotic predictions vs treatment (no combos here) + CIs + SSD:
   python sensityping_metrics_v2.0.py -i results.tsv -o out.txt -d ./out \\
      --analysis_type predicted_vs_treatment \\
      --ci_flag --ci_method hybrid --ci_level 0.95 --n_boot 2000 --seed 1 \\
      --ssd_flag --ssd_width 0.05 --ssd_mode conservative \\
      --radar_flag --radar_metrics PPV,one_minus_FDR,coverage_fraction

2) First-line (rule-based) vs treatment with order including a combo token present in the table:
   python sensityping_metrics_v2.0.py -i treatment_output.tsv -o firstline_metrics.txt -d ./out \\
      --analysis_type first_line_vs_treatment \\
      --order CIP,CRO+AZM,CRO,SPC \\
      --ci_flag --ci_method hybrid \\
      --ssd_flag --ssd_width 0.05 \\
      --radar_flag --radar_metrics PPV,one_minus_FDR,coverage_fraction \\
      --id_extraction CIP,CRO+AZM

Notes
-----
- CI: proportion metrics use Wilson intervals; composite metrics use bootstrap if --ci_method is bootstrap/hybrid.
- SSD width (--ssd_width) is the *half-width* w of the desired CI (e.g., 0.05 means ±0.05, total width 0.10).
- SSD mode:
    conservative: uses p=0.5 (max variance; most conservative)
    observed:     uses p-hat
    fixed:        uses --ssd_p
"""

    parser = argparse.ArgumentParser(
        description="Metrics for antibiotic predictions and first-line recommendations (supports provided combo columns like CRO+AZM).",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=epilog
    )

    parser.add_argument("--version", action="store_true",
                        help="Print version and exit.")

    # Files / dirs
    parser.add_argument("-i", "--input", required=False,
                        help=("Input TSV. For 'predicted_vs_treatment' needs <ABX>_predicted and <ABX>_treatment. "
                              "For 'first_line_vs_treatment' needs <TOKEN>_recommend and <TOKEN>_treatment "
                              "(TOKEN may be 'CRO+AZM')."))
    parser.add_argument("-o", "--output", required=False,
                        help="Output text file for metrics.")
    parser.add_argument("-d", "--output_dir", required=False,
                        help="Directory to write radar HTML and ID tables.")

    # Analysis mode
    parser.add_argument("--analysis_type", required=False,
                        choices=['predicted_vs_treatment', 'first_line_vs_treatment'],
                        help=("Type of analysis:\n"
                              "  predicted_vs_treatment : compare <ABX>_predicted vs <ABX>_treatment per antibiotic\n"
                              "  first_line_vs_treatment : sequentially compare <TOKEN>_recommend vs <TOKEN>_treatment "
                              "using the user-supplied --order; TOKEN may be a combo like CRO+AZM"))

    parser.add_argument("--order",
                        help=("Comma-separated order of tokens (REQUIRED for first_line_vs_treatment). "
                              "Tokens can be single ABX (CIP) or combos present as columns (CRO+AZM). "
                              "Example: CIP,CRO+AZM,CRO,SPC"))

    # Radar settings
    parser.add_argument("--radar_flag", action="store_true",
                        help="If set, writes a radar chart HTML to --output_dir.")
    parser.add_argument("--radar_antibiotics", default="CRO,AZM,CIP,TET,PCG,SPC,ZOL",
                        help="Comma-separated list for predicted_vs_treatment radar (single drugs only).")
    parser.add_argument("--radar_metrics",
                        default="PPV,one_minus_FDR,coverage_fraction",
                        help=("Comma-separated metrics for radar (0..1 scale). Useful options: "
                              "PPV, one_minus_FDR, coverage_fraction, accuracy, one_minus_ME, one_minus_VME"))

    # Cost-sensitive weights
    parser.add_argument("--fp_weight", type=int, default=5,
                        help="Weight for FP in cost-sensitive error (default 5).")
    parser.add_argument("--fn_weight", type=int, default=1,
                        help="Weight for FN in cost-sensitive error (default 1).")

    # Semantics of 'Yes'
    parser.add_argument("--yes_means", choices=['susceptible', 'resistant'], default='susceptible',
                        help="Meaning of 'Yes' in your labels. Default 'susceptible' for *_treatment.")

    # ID extraction
    parser.add_argument("--id_extraction",
                        help="Comma-separated tokens (e.g., CIP,SPC or CRO+AZM) for which to write <TOKEN>_id_extracted.tab with TP/TN/FP/FN IDs.")
    parser.add_argument("--id_column",
                        help="Name of the ID column; if omitted, the first column of the input is used.")

    # CI controls (escape %% in help!)
    parser.add_argument("--ci_flag", action="store_true",
                        help="If set, adds 95%% CIs (proportions via Wilson; composites via bootstrap if enabled).")
    parser.add_argument("--ci_method", choices=["wilson", "bootstrap", "hybrid"], default="hybrid",
                        help="CI method: wilson (proportions only), bootstrap (bootstrap for composites), hybrid (recommended).")
    parser.add_argument("--ci_level", type=float, default=0.95,
                        help="Confidence level (default 0.95).")
    parser.add_argument("--n_boot", type=int, default=2000,
                        help="Bootstrap replicates for composite-metric CIs (default 2000).")
    parser.add_argument("--seed", type=int, default=1,
                        help="Random seed for bootstrap (default 1).")

    # SSD controls
    parser.add_argument("--ssd_flag", action="store_true",
                        help="If set, computes SSD required sample size for specified CI half-width (approx, normal).")
    parser.add_argument("--ssd_width", type=float, default=0.05,
                        help="SSD CI half-width w (default 0.05 => ±0.05; total width 0.10).")
    parser.add_argument("--ssd_mode", choices=["conservative", "observed", "fixed"], default="conservative",
                        help="SSD p choice: conservative (p=0.5), observed (p-hat), or fixed (use --ssd_p).")
    parser.add_argument("--ssd_p", type=float, default=None,
                        help="If --ssd_mode fixed, set p used in SSD formula (0..1).")

    args = parser.parse_args()

    if args.version:
        print("sensityping_metrics_v2.0.py version %s" % __version__)
        raise SystemExit(0)

    # Enforce required args only when actually running analysis
    missing = []
    for req in ["input", "output", "output_dir", "analysis_type"]:
        if getattr(args, req) is None:
            missing.append(req)
    if missing:
        parser.error("Missing required arguments: %s" % ", ".join(missing))

    main(args.input, args.output, args.output_dir,
         args.analysis_type, args.order,
         args.radar_flag, args.radar_antibiotics, args.radar_metrics,
         args.fp_weight, args.fn_weight, args.yes_means,
         args.id_extraction, args.id_column,
         args.ci_flag, args.ci_method, args.ci_level, args.n_boot, args.seed,
         args.ssd_flag, args.ssd_width, args.ssd_mode, args.ssd_p)
