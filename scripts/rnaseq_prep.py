#!/usr/bin/env python3
"""Robust RNA-seq cohort prep: counts (+ optional Ensembl->symbol annotation and a
design file for grouping) -> expr.csv (gene SYMBOLS x samples) + meta.csv.

Handles the real GEO layouts of GSE134692 and GSE150910.

Examples
--------
GSE134692 (raw counts, Ensembl index, design file with DiseaseStatus):
  python scripts/rnaseq_prep.py --gse GSE134692 \
    --counts data/GSE134692/GSE134692_raw_counts.txt.gz --sep tab --index-col ensembl \
    --annot data/GSE134692/GSE134692_gene_annotation.txt.gz --annot-id ensembl --annot-sym GeneName \
    --design data/GSE134692/GSE134692_design.txt.gz --design-sample sample_id --design-group DiseaseStatus \
    --case IPF --control Normal

GSE150910 (gene-level csv; grouping from GEO series-matrix phenotype):
  python scripts/rnaseq_prep.py --gse GSE150910 \
    --counts data/GSE150910/GSE150910_gene-level_count_file.csv.gz --sep comma \
    --geo-group          # auto-match columns to GEO diagnosis
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
IPF_KEYS = ["ipf", "idiopathic pulmonary fibrosis", "uip"]
CTRL_KEYS = ["control", "normal", "donor", "non-fibrotic", "healthy", "unused"]
SEPS = {"tab": "\t", "comma": ",", "auto": None}


def read_any(path, sep, index_col=0):
    s = SEPS.get(sep, sep)
    if s is None:
        s = "\t" if str(path).endswith((".txt.gz", ".tsv.gz", ".txt", ".tsv")) else ","
    return pd.read_csv(path, sep=s, index_col=index_col, compression="infer",
                       low_memory=False)


def to_symbols(counts, annot_path, id_col, sym_col, sep):
    annot = read_any(annot_path, sep, index_col=None)
    m = annot.set_index(id_col)[sym_col].astype(str)
    counts = counts.copy()
    counts.index = counts.index.map(m)
    counts = counts[~counts.index.isin(["", "nan", "NaN", "None"])]
    counts = counts.groupby(counts.index).sum()   # collapse to symbol level
    return counts


def classify(value, case, control):
    lv = str(value).lower()
    if case and case.lower() in lv:
        return "IPF"
    if control and control.lower() in lv:
        return "control"
    if any(k in lv for k in IPF_KEYS):
        return "IPF"
    if any(k in lv for k in CTRL_KEYS):
        return "control"
    return np.nan


def group_from_design(path, sample_col, group_col, case, control, sep):
    d = read_any(path, sep, index_col=None)
    lab = d[group_col].map(lambda v: classify(v, case, control))
    g = pd.Series(lab.values, index=d[sample_col].astype(str))
    return g.dropna()


def group_from_geo(gse, columns, case, control):
    import GEOparse
    g = GEOparse.get_GEO(gse, destdir=str(ROOT / "data" / gse), silent=True)
    p = g.phenotype_data
    ch = [c for c in p.columns if "characteristics" in c or "title" in c
          or "description" in c or "source" in c]
    lab = p[ch].apply(lambda r: classify(" ".join(map(str, r)), case, control), axis=1).dropna()
    # match GEO sample identifiers to count columns via any phenotype text field
    colset = list(columns)
    out = {}
    for gsm, grp in lab.items():
        # candidate identifiers for this sample
        cands = [gsm] + [str(p.loc[gsm, c]) for c in p.columns
                         if any(k in c for k in ["title", "description", "supplementary", "geo"])]
        for col in colset:
            if any(col == c or col in c or c in col for c in cands):
                out[col] = grp
                break
    return pd.Series(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gse", required=True)
    ap.add_argument("--counts", required=True)
    ap.add_argument("--sep", default="auto", help="tab|comma|auto or a literal char")
    ap.add_argument("--index-col", default=0, help="counts index column (name or 0)")
    ap.add_argument("--annot"); ap.add_argument("--annot-id"); ap.add_argument("--annot-sym")
    ap.add_argument("--design"); ap.add_argument("--design-sample"); ap.add_argument("--design-group")
    ap.add_argument("--geo-group", action="store_true",
                    help="derive grouping from GEO series-matrix phenotype")
    ap.add_argument("--case", default="IPF")
    ap.add_argument("--control", default="control")
    ap.add_argument("--normalized", action="store_true",
                    help="counts are already log/normalized; do NOT CPM-log2 again")
    args = ap.parse_args()

    idx = args.index_col
    try:
        idx = int(idx)
    except (ValueError, TypeError):
        pass
    print(f"=== {args.gse} rnaseq prep ===")
    counts = read_any(args.counts, args.sep, index_col=idx)
    counts = counts.apply(pd.to_numeric, errors="coerce").dropna(how="all")
    print(f"  counts: {counts.shape[0]} features x {counts.shape[1]} samples")

    if args.annot and args.annot_id and args.annot_sym:
        counts = to_symbols(counts, args.annot, args.annot_id, args.annot_sym, args.sep)
        print(f"  mapped to {counts.shape[0]} gene symbols")

    # grouping
    if args.design:
        groups = group_from_design(args.design, args.design_sample, args.design_group,
                                   args.case, args.control, args.sep)
    elif args.geo_group:
        groups = group_from_geo(args.gse, counts.columns, args.case, args.control)
    else:
        raise SystemExit("Provide --design ... or --geo-group for sample grouping.")
    groups = groups[groups.isin(["IPF", "control"])]
    # keep only samples present in both
    common = [s for s in counts.columns if s in groups.index]
    if len(common) < 4:
        raise SystemExit(f"Only {len(common)} samples matched between counts and grouping. "
                         f"counts cols e.g. {list(counts.columns)[:3]}; "
                         f"group ids e.g. {list(groups.index)[:3]}")
    counts = counts[common]
    groups = groups.reindex(common)
    print(f"  matched {len(common)} samples: "
          f"{(groups=='IPF').sum()} IPF / {(groups=='control').sum()} control")

    if args.normalized:
        expr = counts
    else:
        lib = counts.sum(axis=0).replace(0, np.nan)
        expr = np.log2(counts.div(lib, axis=1) * 1e6 + 1.0)

    out = ROOT / "data" / args.gse
    expr.to_csv(out / "expr.csv")
    pd.DataFrame({"diagnosis": groups}).to_csv(out / "meta.csv")
    print(f"  wrote {out/'expr.csv'} ({expr.shape[0]} genes x {expr.shape[1]} samples)")
    print(f"  wrote {out/'meta.csv'}")


if __name__ == "__main__":
    main()
