#!/usr/bin/env python3
"""Count-based DESeq2 sensitivity analysis for the RNA-seq discovery cohort (GSE150910),
using pydeseq2 (pure Python — no R required). Runs an optional negative-binomial DESeq2 sensitivity analysis
for the RNA-seq discovery cohort and quantifies overlap with the primary
cross-cohort signature.

Install (once), in the ipffx env:
    pip install pydeseq2

Run from the project root:
    python scripts/de_sensitivity_deseq2.py \
        --counts data/GSE150910/GSE150910_gene-level_count_file.csv.gz \
        --primary results_real/tables/S_stage1_replicated_signature.csv \
        --out results_real/tables/S2b_deseq2.csv

Outputs a one-row-per-method summary CSV and prints the overlap of the
DESeq2-based cross-cohort replicated set with the primary (Welch) signature.
"""
from __future__ import annotations
import argparse
import re
from pathlib import Path
import numpy as np
import pandas as pd


def label(col: str):
    c = col.lower()
    if c.startswith("ipf"):
        return "IPF"
    if c.startswith("control"):
        return "control"
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", required=True,
                    help="RNA-seq raw counts, genes x samples (csv/tsv, .gz ok)")
    ap.add_argument("--primary", default="results_real/tables/S_stage1_replicated_signature.csv",
                    help="primary replicated-signature table from the pipeline")
    ap.add_argument("--out", default="results_real/tables/S2b_deseq2.csv")
    ap.add_argument("--alpha", type=float, default=0.05)
    args = ap.parse_args()

    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
    except ImportError:
        raise SystemExit("pydeseq2 not installed. Run:  pip install pydeseq2")

    sep = "\t" if re.search(r"\.(txt|tsv)", args.counts) else ","
    counts = pd.read_csv(args.counts, sep=sep, index_col=0, compression="infer")
    counts = counts.groupby(counts.index).sum()          # collapse duplicate symbols
    grp = pd.Series({c: label(c) for c in counts.columns}).dropna()
    grp = grp[grp.isin(["IPF", "control"])]
    counts = counts[grp.index]
    # DESeq2 needs integer counts, samples x genes
    X = counts.T.round().astype(int)
    X = X.loc[:, X.sum(axis=0) > 0]                       # drop all-zero genes
    meta = pd.DataFrame({"condition": grp.reindex(X.index).values}, index=X.index)
    print(f"DESeq2 on {X.shape[0]} samples x {X.shape[1]} genes "
          f"({(meta.condition=='IPF').sum()} IPF / {(meta.condition=='control').sum()} control)")

    dds = DeseqDataSet(counts=X, metadata=meta, design_factors="condition", quiet=True)
    dds.deseq2()
    st = DeseqStats(dds, contrast=["condition", "IPF", "control"], quiet=True)
    st.summary()
    res = st.results_df.rename(columns={"log2FoldChange": "lfc", "padj": "fdr"})
    res = res.dropna(subset=["fdr"])

    ds_up = set(res.index[(res.fdr < args.alpha) & (res.lfc > 0)])
    ds_dn = set(res.index[(res.fdr < args.alpha) & (res.lfc < 0)])
    print(f"DESeq2 (GSE150910) significant: up={len(ds_up)} down={len(ds_dn)}")

    # cross with the microarray cohort's Welch significance (from the primary table)
    prim = pd.read_csv(args.primary, index_col=0)
    fdr47 = [c for c in prim.columns if c.startswith("fdr_") and "47460" in c]
    lfc47 = [c for c in prim.columns if c.startswith("log2FC_") and "47460" in c]
    pu = set(prim.index[prim["replicated_up"]]); pd_ = set(prim.index[prim["replicated_down"]])
    summary = {"method": "DESeq2 (pydeseq2), RNA-seq cohort + microarray",
               "primary_up": len(pu), "primary_down": len(pd_),
               "deseq2_sig_up_RNAseq": len(ds_up), "deseq2_sig_down_RNAseq": len(ds_dn)}

    if fdr47 and lfc47:
        common = res.index.intersection(prim.index)
        p = prim.loc[common]
        a_sig_up = set(common[(p[fdr47[0]] < args.alpha) & (p[lfc47[0]] > 0)])
        a_sig_dn = set(common[(p[fdr47[0]] < args.alpha) & (p[lfc47[0]] < 0)])
        hy_up = (ds_up & a_sig_up)
        hy_dn = (ds_dn & a_sig_dn)
        ov_up = 100 * len(hy_up & pu) / max(len(hy_up), 1)
        ov_dn = 100 * len(hy_dn & pd_) / max(len(hy_dn), 1)
        # log2FC correlation vs primary mean
        if "mean_log2FC" in prim.columns:
            cc = res.index.intersection(prim.index)
            r = np.corrcoef(res.loc[cc, "lfc"], prim.loc[cc, "mean_log2FC"])[0, 1]
        else:
            r = np.nan
        summary.update({
            "deseq2_hybrid_replicated_up": len(hy_up),
            "deseq2_hybrid_replicated_down": len(hy_dn),
            "pct_hybrid_up_within_primary": round(ov_up, 1),
            "pct_hybrid_down_within_primary": round(ov_dn, 1),
            "log2FC_corr_vs_primary": round(float(r), 3)})
        print(f"DESeq2-hybrid replicated (DESeq2 RNA-seq + Welch microarray): "
              f"up={len(hy_up)} down={len(hy_dn)}")
        print(f"  % within primary signature: up={ov_up:.1f}%  down={ov_dn:.1f}%")
        print(f"  log2FC correlation vs primary: r={r:.3f}")

    outp = Path(args.out); outp.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([summary]).to_csv(outp, index=False)
    res.to_csv(outp.with_name("S2b_deseq2_full_GSE150910.csv"))
    print(f"\nWrote {outp}\nWrote {outp.with_name('S2b_deseq2_full_GSE150910.csv')} (full per-gene DESeq2)")



if __name__ == "__main__":
    main()
