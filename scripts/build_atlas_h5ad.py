#!/usr/bin/env python3
"""Assemble a pseudobulk-ready .h5ad from the GSE135893 (Habermann 2020) atlas
supplementary files, OR validate an atlas you downloaded from the IPF Cell Atlas.

Two paths:

  A) You already have an .h5ad (e.g. from https://www.ipfcellatlas.com/):
       python scripts/build_atlas_h5ad.py --h5ad /path/to/atlas.h5ad
     -> checks it has donor / cell-type / condition columns and normalizes their
        names, writing data/GSE135893/ipf_atlas.h5ad.

  B) Build from GEO supplementary (matrix.mtx + genes + barcodes + metadata):
       python scripts/build_atlas_h5ad.py          # looks in data/GSE135893/

The output has raw counts in .X and .obs columns: donor, cell_type, disease
(disease in {IPF, control}), which is what pseudobulk_from_h5ad expects.
"""
from __future__ import annotations
import argparse
import glob
import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
ATLAS_DIR = ROOT / "data" / "GSE135893"
OUT = ATLAS_DIR / "ipf_atlas.h5ad"

DONOR_KEYS = ["donor", "sample", "orig.ident", "patient", "subject", "sample_name"]
CT_KEYS = ["celltype", "cell_type", "cell.type", "cluster", "annotation",
           "population", "manual_annotation"]
COND_KEYS = ["diagnosis", "disease", "condition", "status", "group"]
IPF_VALS = ["ipf", "idiopathic pulmonary fibrosis", "uip"]
CTRL_VALS = ["control", "normal", "donor", "non-fibrotic", "healthy", "ctrl"]


def _find(cols, keys):
    low = {c.lower(): c for c in cols}
    for k in keys:
        for lc, orig in low.items():
            if k in lc:
                return orig
    return None


def normalize_condition(series: pd.Series) -> pd.Series:
    def m(v):
        lv = str(v).lower()
        if any(k in lv for k in IPF_VALS):
            return "IPF"
        if any(k in lv for k in CTRL_VALS):
            return "control"
        return np.nan
    return series.map(m)


def finalize(adata, donor_col=None, cell_col=None, cond_col=None):
    import anndata as ad
    obs = adata.obs
    dcol = donor_col or _find(obs.columns, DONOR_KEYS)
    ccol = cell_col or _find(obs.columns, CT_KEYS)
    ncol = cond_col or _find(obs.columns, COND_KEYS)
    missing = [n for n, c in [("donor", dcol), ("cell_type", ccol),
                              ("condition", ncol)] if c is None]
    if missing:
        raise SystemExit("Atlas obs is missing " + ", ".join(missing) +
                         f".\nAvailable obs columns:\n" + "\n".join(obs.columns))
    obs2 = pd.DataFrame(index=obs.index)
    obs2["donor"] = obs[dcol].astype(str)
    obs2["cell_type"] = obs[ccol].astype(str)
    obs2["disease"] = normalize_condition(obs[ncol])
    keep = obs2["disease"].notna()
    adata = adata[keep.values].copy()
    adata.obs = obs2[keep]
    print(f"  atlas: {adata.n_obs} cells, {adata.obs['donor'].nunique()} donors, "
          f"{adata.obs['cell_type'].nunique()} cell types, "
          f"IPF/control = "
          f"{(adata.obs['disease']=='IPF').sum()}/{(adata.obs['disease']=='control').sum()}")
    ATLAS_DIR.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(OUT)
    print(f"  wrote {OUT}")


def from_existing(path, **cols):
    import anndata as ad
    print(f"Validating existing atlas {path}")
    finalize(ad.read_h5ad(path), **cols)


def from_geo(**cols):
    import scanpy as sc
    import anndata as ad
    from scipy.io import mmread
    import scipy.sparse as sp
    print(f"Assembling atlas from GEO suppl in {ATLAS_DIR}")
    mtx = glob.glob(str(ATLAS_DIR / "*matrix.mtx*"))
    genes = glob.glob(str(ATLAS_DIR / "*genes.tsv*")) + glob.glob(str(ATLAS_DIR / "*features.tsv*"))
    bcs = glob.glob(str(ATLAS_DIR / "*barcodes.tsv*"))
    meta = [f for f in glob.glob(str(ATLAS_DIR / "*metadata*")) + glob.glob(str(ATLAS_DIR / "*.csv*"))
            if "metadata" in f.lower() or "annot" in f.lower()]
    if not (mtx and genes and bcs):
        raise SystemExit(
            "Could not find matrix.mtx/genes.tsv/barcodes.tsv in data/GSE135893/.\n"
            "Easiest fix: download the processed .h5ad from https://www.ipfcellatlas.com/ "
            "and run:  python scripts/build_atlas_h5ad.py --h5ad <file>")
    X = sp.csr_matrix(mmread(mtx[0]).T)  # cells x genes
    gnames = pd.read_csv(genes[0], sep="\t", header=None)
    var = pd.DataFrame(index=gnames.iloc[:, -1].astype(str).values)
    barcodes = pd.read_csv(bcs[0], sep="\t", header=None).iloc[:, 0].astype(str).values
    adata = ad.AnnData(X=X, obs=pd.DataFrame(index=barcodes), var=var)
    if meta:
        md = pd.read_csv(meta[0], index_col=0, compression="infer")
        common = adata.obs_names.intersection(md.index)
        adata = adata[common].copy()
        adata.obs = md.loc[common]
    else:
        raise SystemExit("No metadata CSV found; need per-cell donor/celltype/diagnosis.")
    finalize(adata, **cols)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", help="validate/normalize an existing atlas .h5ad")
    ap.add_argument("--donor-col", help="obs column to use as donor")
    ap.add_argument("--cell-col", help="obs column to use as cell_type")
    ap.add_argument("--cond-col", help="obs column to use as condition/diagnosis")
    args = ap.parse_args()
    cols = {"donor_col": args.donor_col, "cell_col": args.cell_col,
            "cond_col": args.cond_col}
    if args.h5ad:
        from_existing(args.h5ad, **cols)
    else:
        from_geo(**cols)


if __name__ == "__main__":
    main()
