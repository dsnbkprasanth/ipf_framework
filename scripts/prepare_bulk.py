#!/usr/bin/env python3
"""Turn downloaded GEO bulk cohorts into expr.csv + meta.csv the pipeline expects.

Runs after scripts/download_data.sh. For each cohort it:
  * loads phenotype from the series matrix (via GEOparse),
  * auto-detects the diagnosis column and IPF / control labels (override with flags),
  * builds an expression matrix (genes x samples):
      - microarray  -> series-matrix VALUE table, probes mapped to gene symbols,
      - RNA-seq      -> the largest supplementary count table in data/<GSE>/,
  * writes data/<GSE>/expr.csv and data/<GSE>/meta.csv.

Usage:
    python scripts/prepare_bulk.py                       # all default cohorts
    python scripts/prepare_bulk.py --gse GSE150910       # one cohort
    python scripts/prepare_bulk.py --gse GSE134692 \
        --group-col "characteristics_ch1.1.diagnosis" \
        --case IPF --control control                     # manual override
"""
from __future__ import annotations
import argparse
import glob
import gzip
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"

IPF_KEYS = ["ipf", "idiopathic pulmonary fibrosis", "uip", "usual interstitial"]
CTRL_KEYS = ["control", "normal", "donor", "non-fibrotic", "healthy", "unused donor",
             "ctrl", "unaffected"]

DEFAULT_COHORTS = {
    "GSE47460":  "microarray",
    "GSE150910": "rnaseq",
    "GSE134692": "rnaseq",
}


def _open_maybe_gz(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def load_phenotype(gse: str):
    import GEOparse
    dest = DATA / gse
    g = GEOparse.get_GEO(gse, destdir=str(dest), silent=True)
    return g


def detect_group_column(pheno: pd.DataFrame):
    """Find the characteristics column that best separates IPF vs control."""
    best, best_score = None, -1
    for col in pheno.columns:
        vals = pheno[col].astype(str).str.lower()
        joined = " ".join(vals.unique())
        has_ipf = any(k in joined for k in IPF_KEYS)
        has_ctrl = any(k in joined for k in CTRL_KEYS)
        if has_ipf and has_ctrl:
            score = vals.nunique()
            # prefer few categories (a clean diagnosis field)
            score = 100 - min(vals.nunique(), 50)
            if score > best_score:
                best, best_score = col, score
    return best


def map_labels(series: pd.Series):
    out = {}
    for s, v in series.items():
        lv = str(v).lower()
        if any(k in lv for k in IPF_KEYS):
            out[s] = "IPF"
        elif any(k in lv for k in CTRL_KEYS):
            out[s] = "control"
    return pd.Series(out)


def build_meta(gse, group_col=None, case=None, control=None):
    g = load_phenotype(gse)
    pheno = g.phenotype_data
    col = group_col or detect_group_column(pheno)
    if col is None:
        raise SystemExit(f"{gse}: could not auto-detect a diagnosis column. "
                         f"Columns:\n" + "\n".join(pheno.columns) +
                         "\nRe-run with --group-col.")
    if case and control:
        raw = pheno[col].astype(str)
        groups = pd.Series({s: ("IPF" if case.lower() in v.lower()
                                 else "control" if control.lower() in v.lower()
                                 else np.nan)
                            for s, v in raw.items()}).dropna()
    else:
        groups = map_labels(pheno[col])
    groups = groups[groups.isin(["IPF", "control"])]
    print(f"  {gse}: group column '{col}'  ->  "
          f"{(groups=='IPF').sum()} IPF / {(groups=='control').sum()} control")
    meta = pd.DataFrame({"diagnosis": groups})
    return g, meta


def expr_microarray(g, meta):
    """Series-matrix VALUE table; map probes to gene symbols via the platform."""
    expr = g.pivot_samples("VALUE")            # probes x samples
    # platform annotation for probe -> symbol
    gpl = list(g.gpls.values())[0]
    tab = gpl.table
    sym_col = next((c for c in tab.columns
                    if c.upper() in ("GENE_SYMBOL", "GENE SYMBOL", "SYMBOL",
                                     "GENESYMBOL", "ILMN_GENE")), None)
    id_col = tab.columns[0]
    if sym_col:
        mapping = tab.set_index(id_col)[sym_col].astype(str)
        expr = expr.copy()
        expr.index = expr.index.map(mapping)
        expr = expr[~expr.index.isin(["", "nan", "NaN", "---"])]
        expr = expr.groupby(expr.index).mean()   # collapse probes per symbol
    expr = expr.reindex(columns=meta.index)
    return expr.dropna(how="all")


def expr_rnaseq(gse, meta):
    """Load the largest supplementary count table and coerce to genes x samples."""
    files = [f for f in glob.glob(str(DATA / gse / "*"))
             if re.search(r"\.(txt|tsv|csv)(\.gz)?$", f)
             and "series_matrix" not in f]
    if not files:
        raise SystemExit(f"{gse}: no supplementary count table found in data/{gse}/. "
                         f"Check the download or pass one manually.")
    files.sort(key=lambda f: os.path.getsize(f), reverse=True)
    path = files[0]
    print(f"  {gse}: using count table {os.path.basename(path)}")
    sep = "\t" if re.search(r"\.(txt|tsv)", path) else ","
    df = pd.read_csv(path, sep=sep, index_col=0, compression="infer")
    # drop non-numeric annotation columns
    df = df.apply(pd.to_numeric, errors="coerce").dropna(how="all")
    # align samples to metadata (match by substring if names differ)
    cols = {}
    for s in meta.index:
        if s in df.columns:
            cols[s] = s
        else:
            hit = [c for c in df.columns if s in c or c in s]
            if hit:
                cols[s] = hit[0]
    if len(cols) < 4:
        print(f"    WARNING: only {len(cols)} samples matched by name. You may need "
              f"to reconcile sample IDs between meta and the count table.")
    df = df[[cols[s] for s in meta.index if s in cols]]
    df.columns = [s for s in meta.index if s in cols]
    # counts -> log2 CPM
    lib = df.sum(axis=0).replace(0, np.nan)
    df = np.log2(df.div(lib, axis=1) * 1e6 + 1.0)
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gse", action="append", help="cohort(s); default = all")
    ap.add_argument("--kind", choices=["microarray", "rnaseq"],
                    help="override platform type for a single --gse")
    ap.add_argument("--group-col")
    ap.add_argument("--case")
    ap.add_argument("--control")
    args = ap.parse_args()

    cohorts = {k: DEFAULT_COHORTS[k] for k in (args.gse or DEFAULT_COHORTS)}
    if args.kind and args.gse and len(args.gse) == 1:
        cohorts[args.gse[0]] = args.kind

    for gse, kind in cohorts.items():
        print(f"\n=== {gse} ({kind}) ===")
        g, meta = build_meta(gse, args.group_col, args.case, args.control)
        expr = (expr_microarray(g, meta) if kind == "microarray"
                else expr_rnaseq(gse, meta))
        meta = meta.reindex(expr.columns)
        out = DATA / gse
        expr.to_csv(out / "expr.csv")
        meta.to_csv(out / "meta.csv")
        print(f"  wrote {out/'expr.csv'}  ({expr.shape[0]} genes x {expr.shape[1]} samples)")
        print(f"  wrote {out/'meta.csv'}")

    print("\nDone. Group column values must contain IPF and control. "
          "Now run: python scripts/build_atlas_h5ad.py  (or place ipf_atlas.h5ad), "
          "then bash scripts/run_real.sh")


if __name__ == "__main__":
    main()
