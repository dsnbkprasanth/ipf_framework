"""Data loading and a synthetic IPF-like demo generator.

Two ways to feed the pipeline:

* REAL data — use the loader functions (load_bulk_cohort, pseudobulk_from_h5ad,
  load_lincs_matrix). Point them at files you downloaded (see
  scripts/DOWNLOAD_DATA.md). Nothing is bundled with this package.

* DEMO data — make_demo_bundle() fabricates a small, biologically-styled IPF-like
  dataset for every layer so the whole 7-stage workflow runs with no downloads.
  Demo results are clearly labelled and must never be interpreted as real
  findings.

All layers are returned in one dict ("bundle") consumed by ipffx.pipeline.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd

from . import config
from .utils import differential_expression

CELL_TYPES = ["Epithelial", "Fibroblast", "Myeloid", "Endothelial", "Lymphoid"]


# ----------------------------------------------------------------------------- #
# Real-data loaders
# ----------------------------------------------------------------------------- #
def _read_table(path: str | Path) -> pd.DataFrame:
    path = str(path)
    sep = "\t" if path.endswith((".tsv", ".txt", ".tsv.gz", ".txt.gz")) else ","
    return pd.read_csv(path, sep=sep, index_col=0)


def load_bulk_cohort(cohort_id: str, expr_path: str | Path, meta_path: str | Path,
                     group_col: str, case: str, control: str) -> dict:
    """Load one bulk cohort from local files.

    expr_path : genes x samples table (csv/tsv). Values may be counts or already
                log/normalized; counts are auto log2(x+1) transformed.
    meta_path : samples x metadata table containing `group_col`.
    """
    expr = _read_table(expr_path)
    meta = _read_table(meta_path)
    common = [s for s in expr.columns if s in meta.index]
    if not common:
        raise ValueError(f"{cohort_id}: no overlapping samples between expr and meta.")
    expr = expr[common]
    groups = meta.loc[common, group_col].astype(str)
    groups = groups.replace({case: "IPF", control: "control"})
    groups = groups[groups.isin(["IPF", "control"])]
    expr = expr[groups.index]
    # auto log-transform if values look like raw counts
    if np.nanmax(expr.to_numpy()) > 100 and (expr.to_numpy() >= 0).all():
        expr = np.log2(expr + 1.0)
    return {"id": cohort_id, "expr": expr, "meta": groups}


def compute_bulk_de(cohort: dict) -> pd.DataFrame:
    return differential_expression(cohort["expr"], cohort["meta"], "IPF", "control")


def pseudobulk_from_h5ad(path: str | Path, donor_col: str, celltype_col: str,
                         condition_col: str, case: str, control: str,
                         min_cells: int = 10) -> dict:
    """Donor x cell-type pseudobulk, then per-cell-type disease DE.

    Requires `anndata`. Aggregates raw counts to donor-by-cell-type pseudobulk
    samples (sum), CPM-normalizes and log2-transforms, then Welch t-test of
    IPF vs control pseudobulk samples within each cell type. This follows
    multi-sample pseudobulk recommendations and avoids cell-level
    pseudoreplication.
    """
    import anndata as ad
    adata = ad.read_h5ad(path)
    obs = adata.obs
    X = adata.X
    genes = adata.var_names.to_numpy()
    # build donor x celltype pseudobulk
    keys = list(zip(obs[donor_col].astype(str), obs[celltype_col].astype(str)))
    key_index: dict[tuple, int] = {}
    rows, sample_meta = [], []
    import scipy.sparse as sp
    is_sparse = sp.issparse(X)
    for k in dict.fromkeys(keys):
        mask = np.array([kk == k for kk in keys])
        if mask.sum() < min_cells:
            continue
        sub = X[mask]
        summed = np.asarray(sub.sum(axis=0)).ravel() if is_sparse else sub.sum(axis=0)
        rows.append(summed)
        cond = obs.loc[mask, condition_col].astype(str).iloc[0]
        sample_meta.append({"donor": k[0], "cell_type": k[1], "condition": cond})
    pb = pd.DataFrame(np.vstack(rows), columns=genes)
    smeta = pd.DataFrame(sample_meta)
    # CPM + log2
    lib = pb.sum(axis=1).replace(0, np.nan)
    pb = np.log2(pb.div(lib, axis=0) * 1e6 + 1.0)
    de, deg_counts = {}, []
    for ct in smeta["cell_type"].unique():
        idx = smeta.index[smeta["cell_type"] == ct]
        cond = smeta.loc[idx, "condition"].replace({case: "IPF", control: "control"})
        expr_ct = pb.loc[idx].T
        expr_ct.columns = idx
        cond = cond[cond.isin(["IPF", "control"])]
        if (cond == "IPF").sum() < 2 or (cond == "control").sum() < 2:
            continue
        res = differential_expression(expr_ct[cond.index], cond, "IPF", "control")
        de[ct] = res
        sig = res[res["fdr"] < 0.05]
        deg_counts.append({"cell_type": ct,
                           "up": int((sig["log2FC"] > 0).sum()),
                           "down": int((sig["log2FC"] < 0).sum())})
    return {"de": de, "deg_counts": pd.DataFrame(deg_counts), "source": str(path)}


def load_lincs_matrix(path: str | Path) -> pd.DataFrame:
    """Load a LINCS/L1000 signature matrix: compounds (rows) x genes (cols),
    values = signed differential expression / z-scores induced by each compound."""
    return _read_table(path)


# ----------------------------------------------------------------------------- #
# Synthetic demo generator
# ----------------------------------------------------------------------------- #
def _module_gene_universe() -> list[str]:
    mods = config.load_modules()
    genes = set()
    for m in mods.values():
        genes.update(m["genes"])
    return sorted(genes)


def make_demo_bundle(seed: int = 7, n_background: int = 1500) -> dict:
    """Fabricate a coherent IPF-like dataset across all layers.

    The synthetic ground truth: module genes (fibrotic/epithelial/myeloid/TGFb)
    are up in IPF; a curated 'down' set is down. Two discovery cohorts share the
    signal (with noise) so cross-cohort replication is real. Pseudobulk localizes
    the fibrotic program to Fibroblast/Epithelial/Myeloid. LINCS positive controls
    (nintedanib, pirfenidone) reverse the disease signature.
    """
    rng = np.random.default_rng(seed)
    mods = config.load_modules()
    up_genes = _module_gene_universe()
    down_genes = [f"DN{i:03d}" for i in range(120)]  # synthetic down-regulated set
    bg_genes = [f"BG{i:04d}" for i in range(n_background)]
    all_genes = up_genes + down_genes + bg_genes

    # "true" IPF effect per gene (log2FC)
    true_fc = pd.Series(0.0, index=all_genes)
    true_fc[up_genes] = rng.normal(1.4, 0.5, len(up_genes)).clip(0.3, None)
    true_fc[down_genes] = rng.normal(-1.2, 0.4, len(down_genes)).clip(None, -0.3)
    true_fc[bg_genes] = rng.normal(0.0, 0.15, len(bg_genes))

    def make_cohort(cid, n_case, n_ctrl, fc_scale, noise):
        base = pd.Series(rng.uniform(4, 10, len(all_genes)), index=all_genes)
        cols, meta = [], {}
        for i in range(n_case):
            s = f"{cid}_IPF_{i:02d}"
            cols.append(base + true_fc * fc_scale + rng.normal(0, noise, len(all_genes)))
            meta[s] = "IPF"
        for i in range(n_ctrl):
            s = f"{cid}_CTL_{i:02d}"
            cols.append(base + rng.normal(0, noise, len(all_genes)))
            meta[s] = "control"
        expr = pd.DataFrame({c: v for c, v in zip(meta.keys(), cols)}, index=all_genes)
        return {"id": cid, "expr": expr, "meta": pd.Series(meta)}

    bulk = {
        "GSE47460": make_cohort("GSE47460", 30, 20, 1.0, 0.8),
        "GSE150910": make_cohort("GSE150910", 25, 18, 0.9, 0.9),
    }
    bulk_external = make_cohort("GSE134692", 15, 12, 0.7, 1.0)  # weaker => etiology-ish

    # --- pseudobulk: per-cell-type disease log2FC ---
    # localization weights: fibrotic program strongest in Fibroblast/Epithelial/Myeloid
    ct_weight = {
        "Fibroblast":  {"myofibroblast_ecm": 1.3, "fibroinflammatory_core": 1.1,
                        "tgfb_signaling": 1.0, "profibrotic_macrophage": 0.2,
                        "aberrant_epithelial": 0.2, "epithelial_stress": 0.3},
        "Epithelial":  {"aberrant_epithelial": 1.4, "epithelial_stress": 1.1,
                        "tgfb_signaling": 0.7, "myofibroblast_ecm": 0.3,
                        "fibroinflammatory_core": 0.4, "profibrotic_macrophage": 0.1},
        "Myeloid":     {"profibrotic_macrophage": 1.4, "fibroinflammatory_core": 1.0,
                        "tgfb_signaling": 0.6, "myofibroblast_ecm": 0.2,
                        "aberrant_epithelial": 0.1, "epithelial_stress": 0.3},
        "Endothelial": {"tgfb_signaling": 0.4, "fibroinflammatory_core": 0.3,
                        "myofibroblast_ecm": 0.2, "profibrotic_macrophage": 0.2,
                        "aberrant_epithelial": 0.1, "epithelial_stress": 0.2},
        "Lymphoid":    {"fibroinflammatory_core": 0.4, "profibrotic_macrophage": 0.3,
                        "tgfb_signaling": 0.3, "myofibroblast_ecm": 0.1,
                        "aberrant_epithelial": 0.1, "epithelial_stress": 0.2},
    }
    gene2mod = {}
    for mid, m in mods.items():
        for g in m["genes"]:
            gene2mod.setdefault(g, []).append(mid)

    from .utils import differential_expression
    pb_de = {}
    deg_counts = []
    n_ipf_donors, n_ctl_donors = 14, 10
    for ct in CELL_TYPES:
        # per-cell-type true disease log2FC (localized by module weights)
        fc = pd.Series(rng.normal(0, 0.15, len(all_genes)), index=all_genes)
        for g in up_genes:
            w = max(ct_weight[ct].get(mid, 0.0) for mid in gene2mod.get(g, ["_"]))
            fc[g] = true_fc[g] * w
        for g in down_genes:
            fc[g] = true_fc[g] * (0.6 + 0.4 * rng.random())
        # simulate donor-level pseudobulk samples (log-CPM-like), run true DE
        base = pd.Series(rng.uniform(3, 9, len(all_genes)), index=all_genes)
        noise = 0.55  # donor-to-donor biological variability
        cols, meta = {}, {}
        for i in range(n_ipf_donors):
            s = f"{ct}_IPF_{i:02d}"
            cols[s] = base + fc + rng.normal(0, noise, len(all_genes))
            meta[s] = "IPF"
        for i in range(n_ctl_donors):
            s = f"{ct}_CTL_{i:02d}"
            cols[s] = base + rng.normal(0, noise, len(all_genes))
            meta[s] = "control"
        expr = pd.DataFrame(cols, index=all_genes)
        res = differential_expression(expr, pd.Series(meta), "IPF", "control")
        pb_de[ct] = res
        sig = res[res["fdr"] < 0.05]
        deg_counts.append({"cell_type": ct,
                           "up": int((sig["log2FC"] > 0).sum()),
                           "down": int((sig["log2FC"] < 0).sum())})
    pseudobulk = {"de": pb_de, "deg_counts": pd.DataFrame(deg_counts),
                  "source": "SYNTHETIC-DEMO"}

    # --- snRNA external: module expression shift (IPF - control) by cell type ---
    snrna_rows = []
    for ct in CELL_TYPES:
        for mid, m in mods.items():
            present = [g for g in m["genes"] if g in pb_de[ct].index]
            shift = pb_de[ct].loc[present, "log2FC"].mean()
            # external cohort: add etiology heterogeneity/noise, damp myeloid
            damp = 0.5 if (ct == "Myeloid") else 1.0
            snrna_rows.append({"cell_type": ct, "module": mid,
                               "shift": shift * damp + rng.normal(0, 0.25)})
    snrna_external = pd.DataFrame(snrna_rows)

    # --- LINCS: compound x gene signature matrix ---
    compounds = {
        "nintedanib": ("pdgf", None), "pirfenidone": ("tgfb", None),
        "saracatinib": ("src", None), "imatinib": ("pdgf", None),
        "sirolimus": ("mtor", None), "vorinostat": ("hdac", None),
        "dasatinib": ("src", "cytotoxic_chemo"), "fasudil": ("rho_rock", None),
        "ruxolitinib": ("jak", "immunosuppressant_high"),
        "metformin": (None, None), "simvastatin": (None, None),
        "doxorubicin": (None, "cardiotoxic"), "bleomycin": (None, "pulmonary_toxic"),
        "haloperidol": (None, "neuropsychiatric"),
        "dexamethasone": (None, "steroid_like"), "gefitinib": ("pdgf", None),
        "losartan": ("tgfb", None), "quercetin": (None, None),
        "thalidomide": (None, "neuropsychiatric"), "colchicine": (None, "cytotoxic_chemo"),
    }
    lincs = pd.DataFrame(0.0, index=list(compounds), columns=all_genes)
    disease_vec = true_fc  # what disease does
    for cmpd, (tgt, risk) in compounds.items():
        noise = rng.normal(0, 0.6, len(all_genes))
        # antifibrotic/target compounds reverse the disease vector; others random
        if tgt in ("tgfb", "pdgf", "src", "mtor", "hdac", "rho_rock", "jak"):
            reversal = -disease_vec * rng.uniform(0.6, 1.1)
        elif cmpd in ("metformin", "simvastatin", "quercetin", "losartan"):
            reversal = -disease_vec * rng.uniform(0.3, 0.6)
        else:
            reversal = disease_vec * rng.uniform(-0.2, 0.4)
        lincs.loc[cmpd] = (reversal + noise).values

    drug_meta = pd.DataFrame([
        {"compound": c, "target_class": t or "", "risk_class": r or ""}
        for c, (t, r) in compounds.items()
    ]).set_index("compound")

    return {
        "bulk": bulk,
        "bulk_external": bulk_external,
        "pseudobulk": pseudobulk,
        "snrna_external": snrna_external,
        "lincs": lincs,
        "drug_meta": drug_meta,
        "is_demo": True,
        "up_genes": up_genes,
        "down_genes": down_genes,
    }
