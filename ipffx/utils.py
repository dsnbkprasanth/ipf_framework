"""Shared statistics helpers."""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy import stats


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR. Uses scipy if available, else manual."""
    pvals = np.asarray(pvals, dtype=float)
    try:
        return stats.false_discovery_control(pvals, method="bh")
    except Exception:
        n = len(pvals)
        order = np.argsort(pvals)
        ranked = pvals[order] * n / (np.arange(n) + 1)
        ranked = np.minimum.accumulate(ranked[::-1])[::-1]
        out = np.empty(n)
        out[order] = np.clip(ranked, 0, 1)
        return out


def differential_expression(expr: pd.DataFrame, groups: pd.Series,
                            case_label, control_label,
                            log_input: bool = True) -> pd.DataFrame:
    """Per-gene Welch t-test between case and control.

    expr : genes (rows) x samples (columns), log-space (e.g. log2 CPM / microarray).
    groups : sample -> group label, indexed like expr.columns.
    Returns DataFrame indexed by gene with columns:
      log2FC, t, pvalue, fdr, mean_case, mean_control.
    log2FC = mean_case - mean_control (assumes log-space input).
    """
    groups = groups.reindex(expr.columns)
    case = expr.loc[:, groups == case_label]
    ctrl = expr.loc[:, groups == control_label]
    if case.shape[1] < 2 or ctrl.shape[1] < 2:
        raise ValueError("Need >=2 samples per group for DE.")

    a = case.to_numpy(dtype=float)
    b = ctrl.to_numpy(dtype=float)
    mean_case = np.nanmean(a, axis=1)
    mean_ctrl = np.nanmean(b, axis=1)
    if log_input:
        log2fc = mean_case - mean_ctrl
    else:
        eps = 1e-9
        log2fc = np.log2((mean_case + eps) / (mean_ctrl + eps))
    t, p = stats.ttest_ind(a, b, axis=1, equal_var=False, nan_policy="omit")
    p = np.asarray(p, dtype=float)
    p[np.isnan(p)] = 1.0
    fdr = bh_fdr(p)
    return pd.DataFrame(
        {"log2FC": log2fc, "t": np.asarray(t, dtype=float),
         "pvalue": p, "fdr": fdr,
         "mean_case": mean_case, "mean_control": mean_ctrl},
        index=expr.index,
    )


def zscore_rows(df: pd.DataFrame) -> pd.DataFrame:
    return df.sub(df.mean(axis=1), axis=0).div(df.std(axis=1).replace(0, np.nan), axis=0)
