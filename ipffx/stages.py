"""The seven analytical stages of the IPF framework.

Each function takes the data bundle (and prior-stage outputs where needed) and
returns a dict of result tables (pandas objects). Figures are produced separately
in ipffx.figures. No stage claims causality or therapeutic efficacy.
"""
from __future__ import annotations
import numpy as np
import pandas as pd

from . import config
from .data import compute_bulk_de, CELL_TYPES


# --------------------------------------------------------------------------- #
# Stage 1 — cross-cohort bulk differential expression + reproducible signature
# --------------------------------------------------------------------------- #
def stage1_cross_cohort(bundle: dict, fdr_thresh: float = 0.05) -> dict:
    cohorts = bundle["bulk"]
    ids = list(cohorts)
    if len(ids) < 2:
        raise ValueError("Stage 1 needs >=2 discovery cohorts.")
    de = {cid: compute_bulk_de(c) for cid, c in cohorts.items()}
    common = set(de[ids[0]].index)
    for cid in ids[1:]:
        common &= set(de[cid].index)
    common = sorted(common)

    merged = pd.DataFrame(index=common)
    for cid in ids:
        merged[f"log2FC_{cid}"] = de[cid].loc[common, "log2FC"]
        merged[f"fdr_{cid}"] = de[cid].loc[common, "fdr"]
    fc_cols = [f"log2FC_{cid}" for cid in ids]
    fdr_cols = [f"fdr_{cid}" for cid in ids]
    merged["mean_log2FC"] = merged[fc_cols].mean(axis=1)

    same_dir = (merged[fc_cols].gt(0).all(axis=1) | merged[fc_cols].lt(0).all(axis=1))
    sig_all = merged[fdr_cols].lt(fdr_thresh).all(axis=1)
    up = sig_all & merged[fc_cols].gt(0).all(axis=1)
    down = sig_all & merged[fc_cols].lt(0).all(axis=1)
    merged["replicated_up"] = up
    merged["replicated_down"] = down
    merged["same_direction"] = same_dir

    summary = {
        "n_tested_both": len(common),
        "n_same_direction": int(same_dir.sum()),
        "n_replicated_up": int(up.sum()),
        "n_replicated_down": int(down.sum()),
    }
    return {"de": de, "merged": merged, "summary": summary,
            "replicated_up": merged.index[up].tolist(),
            "replicated_down": merged.index[down].tolist(),
            "cohort_ids": ids}


# --------------------------------------------------------------------------- #
# Stage 2 — donor-level cell-type pseudobulk localization
# --------------------------------------------------------------------------- #
def stage2_pseudobulk(bundle: dict, stage1: dict, fdr_thresh: float = 0.05) -> dict:
    pb = bundle["pseudobulk"]
    de = pb["de"]
    deg = pb["deg_counts"].set_index("cell_type")

    # How much of the replicated bulk signature localizes to each cell type?
    rep_up = set(stage1["replicated_up"])
    rep_down = set(stage1["replicated_down"])
    loc_rows = []
    for ct, res in de.items():
        sig = res[res["fdr"] < fdr_thresh]
        ct_up = set(sig.index[sig["log2FC"] > 0])
        ct_down = set(sig.index[sig["log2FC"] < 0])
        loc_rows.append({
            "cell_type": ct,
            "bulk_up_localized": len(rep_up & ct_up),
            "bulk_down_localized": len(rep_down & ct_down),
        })
    localization = pd.DataFrame(loc_rows).set_index("cell_type")
    return {"deg_counts": deg, "localization": localization, "de": de}


# --------------------------------------------------------------------------- #
# Stage 3 — curated module construction and scoring
# --------------------------------------------------------------------------- #
def stage3_modules(bundle: dict) -> dict:
    mods = config.load_modules()
    pb_de = bundle["pseudobulk"]["de"]

    # module shift = mean disease log2FC across present module genes, per cell type
    shift = pd.DataFrame(index=list(mods), columns=list(pb_de), dtype=float)
    n_present = pd.DataFrame(index=list(mods), columns=list(pb_de), dtype=float)
    for mid, m in mods.items():
        for ct, res in pb_de.items():
            present = [g for g in m["genes"] if g in res.index]
            shift.loc[mid, ct] = res.loc[present, "log2FC"].mean() if present else np.nan
            n_present.loc[mid, ct] = len(present)

    # bulk sample-level module scores (mean module-gene expression), IPF - control
    bulk_scores = {}
    for cid, c in bundle["bulk"].items():
        expr, meta = c["expr"], c["meta"]
        rows = {}
        for mid, m in mods.items():
            present = [g for g in m["genes"] if g in expr.index]
            score = expr.loc[present].mean(axis=0)
            rows[mid] = score
        smat = pd.DataFrame(rows)  # samples x modules
        diff = (smat[meta == "IPF"].mean() - smat[meta == "control"].mean())
        bulk_scores[cid] = diff
    bulk_module_diff = pd.DataFrame(bulk_scores)

    # correlation of TGFb coordinate with other modules across cell types
    corr = shift.T.corr(method="spearman")

    labels = {mid: m["label"] for mid, m in mods.items()}
    return {"module_shift": shift, "n_present": n_present,
            "bulk_module_diff": bulk_module_diff, "module_corr": corr,
            "labels": labels}


# --------------------------------------------------------------------------- #
# Stage 4 — focused ligand-receptor prioritization
# --------------------------------------------------------------------------- #
def stage4_ligand_receptor(bundle: dict) -> dict:
    axes = config.load_ligand_receptor()
    pb_de = bundle["pseudobulk"]["de"]

    def gene_fc(gene, cell_types):
        vals = []
        for ct in cell_types:
            res = pb_de.get(_match_ct(ct, pb_de))
            if res is not None and gene in res.index:
                vals.append(res.loc[gene, "log2FC"])
        return np.mean(vals) if vals else np.nan

    rows = []
    for ax in axes:
        lig_fc = gene_fc(ax["ligand"], ax["sender"])
        rec_fc = gene_fc(ax["receptor"], ax["receiver"])
        # evidence score: ligand up + receptor present/up + curated route
        score = 0.0
        if not np.isnan(lig_fc):
            score += max(lig_fc, 0) * 1.0
        if not np.isnan(rec_fc):
            score += max(rec_fc, 0) * 0.5 + 0.25  # receptor present credit
        score += 0.25  # curated-route credit
        rows.append({
            "axis": f"{ax['ligand']}-{ax['receptor']}",
            "route": ax["route"],
            "ligand_log2FC": lig_fc, "receptor_log2FC": rec_fc,
            "sender": ",".join(ax["sender"]), "receiver": ",".join(ax["receiver"]),
            "evidence_score": score,
        })
    tab = pd.DataFrame(rows).sort_values("evidence_score", ascending=False)
    return {"axes": tab.reset_index(drop=True)}


def _match_ct(name: str, de: dict) -> str | None:
    """Map an L-R compartment label to the closest pseudobulk cell type."""
    name = name.lower()
    mapping = {"macrophage": "Myeloid", "myeloid": "Myeloid",
               "fibroblast": "Fibroblast", "epithelial": "Epithelial",
               "endothelial": "Endothelial", "immune": "Lymphoid",
               "lymphoid": "Lymphoid"}
    target = mapping.get(name)
    if target and target in de:
        return target
    for ct in de:
        if name in ct.lower():
            return ct
    return None


# --------------------------------------------------------------------------- #
# Stage 5 — external validation (bulk + snRNA-seq)
# --------------------------------------------------------------------------- #
def stage5_external(bundle: dict) -> dict:
    mods = config.load_modules()
    ext = bundle["bulk_external"]
    expr, meta = ext["expr"], ext["meta"]
    rows = {}
    for mid, m in mods.items():
        present = [g for g in m["genes"] if g in expr.index]
        score = expr.loc[present].mean(axis=0)
        rows[mid] = score.groupby(meta).mean()
    ext_bulk = pd.DataFrame(rows).T  # modules x groups
    ext_bulk["IPF_minus_control"] = ext_bulk.get("IPF", 0) - ext_bulk.get("control", 0)

    snrna = bundle["snrna_external"].pivot(index="module", columns="cell_type",
                                           values="shift")
    return {"external_bulk": ext_bulk, "snrna_external": snrna}


# --------------------------------------------------------------------------- #
# Stage 6 — LINCS / L1000 safety-aware perturbational prioritization
# --------------------------------------------------------------------------- #
def _connectivity(sig_up, sig_down, lincs: pd.DataFrame) -> pd.Series:
    """Reversal score per compound: positive => opposes the disease signature.

    score = mean(compound z over disease-down genes) - mean(z over disease-up genes)
    i.e. a good reverser is LOW where disease is HIGH and HIGH where disease is LOW.
    """
    up = [g for g in sig_up if g in lincs.columns]
    dn = [g for g in sig_down if g in lincs.columns]
    if not up and not dn:
        return pd.Series(0.0, index=lincs.index)
    up_mean = lincs[up].mean(axis=1) if up else 0.0
    dn_mean = lincs[dn].mean(axis=1) if dn else 0.0
    return dn_mean - up_mean


def stage6_perturbation(bundle: dict, stage1: dict, top_n: int = 150) -> dict:
    safety = config.load_drug_safety()
    lincs = bundle["lincs"]
    drug_meta = bundle["drug_meta"]
    mods = config.load_modules()

    # Build the four query signatures (bulk + 3 cell-type pseudobulk)
    merged = stage1["merged"]
    bulk_up = merged.loc[merged["replicated_up"]].sort_values(
        "mean_log2FC", ascending=False).index[:top_n].tolist()
    bulk_down = merged.loc[merged["replicated_down"]].sort_values(
        "mean_log2FC").index[:top_n].tolist()

    pb_de = bundle["pseudobulk"]["de"]
    sigs = {"bulk": (bulk_up, bulk_down)}
    for ct in ["Fibroblast", "Myeloid", "Epithelial"]:
        if ct in pb_de:
            res = pb_de[ct]
            up = res.sort_values("log2FC", ascending=False).index[:top_n].tolist()
            dn = res.sort_values("log2FC").index[:top_n].tolist()
            sigs[ct] = (up, dn)

    conn = pd.DataFrame({name: _connectivity(u, d, lincs)
                         for name, (u, d) in sigs.items()})
    conn["recurrence"] = conn.mean(axis=1)
    conn["multi_signature_support"] = (conn[list(sigs)] > 0).sum(axis=1)

    # module support: does the compound push TGFb + fibroinflammatory genes DOWN?
    mod_support = {}
    for cmpd in lincs.index:
        vals = []
        for mid in ["tgfb_signaling", "fibroinflammatory_core", "myofibroblast_ecm"]:
            genes = [g for g in mods[mid]["genes"] if g in lincs.columns]
            if genes:
                vals.append(-lincs.loc[cmpd, genes].mean())  # negative expr = suppressed
        mod_support[cmpd] = np.mean(vals) if vals else 0.0
    conn["module_support"] = pd.Series(mod_support)

    # known-target interpretability
    def target_credit(cmpd):
        tc = str(drug_meta.loc[cmpd, "target_class"]) if cmpd in drug_meta.index else ""
        return 1.0 if tc and tc in safety["pathway_targets"] else 0.0
    conn["known_target"] = [target_credit(c) for c in conn.index]

    # safety penalty
    def penalty(cmpd):
        rc = str(drug_meta.loc[cmpd, "risk_class"]) if cmpd in drug_meta.index else ""
        return float(safety["risk_classes"].get(rc, 0))
    conn["safety_penalty"] = [penalty(c) for c in conn.index]

    w = safety["weights"]
    # normalize the continuous components to comparable scales
    def z(s):
        s = s.astype(float)
        return (s - s.mean()) / (s.std() + 1e-9)
    conn["total_score"] = (
        w["reversal_recurrence"] * z(conn["recurrence"])
        + w["multi_signature_support"] * z(conn["multi_signature_support"])
        + w["module_support"] * z(conn["module_support"])
        + w["known_target_interpretability"] * conn["known_target"]
        - w["safety_penalty"] * (conn["safety_penalty"] / 4.0)
    )
    conn = conn.sort_values("total_score", ascending=False)

    def categorize(row):
        if row["safety_penalty"] >= 3:
            return "cautionary" if row["recurrence"] > 0 else "deprioritized"
        if row["total_score"] > 0.75 and row["recurrence"] > 0:
            return "prioritized"
        if row["recurrence"] > 0:
            return "supporting"
        return "deprioritized"
    conn["category"] = conn.apply(categorize, axis=1)
    conn.insert(0, "compound", conn.index)

    pos = safety["positive_controls"]
    pos_rank = {p: int(np.where(conn.index == p)[0][0]) + 1
                for p in pos if p in conn.index}
    return {"perturbation": conn.reset_index(drop=True),
            "signatures": {k: {"up": u, "down": d} for k, (u, d) in sigs.items()},
            "positive_control_ranks": pos_rank}


# --------------------------------------------------------------------------- #
# Stage 7 — robustness (leave-one-cohort-out + size-matched random controls)
# --------------------------------------------------------------------------- #
def stage7_robustness(bundle: dict, stage1: dict, n_perm: int = 1000,
                      fdr_thresh: float = 0.05, seed: int = 0) -> dict:
    de = stage1["de"]
    ids = stage1["cohort_ids"]
    rng = np.random.default_rng(seed)

    # Leave-one-cohort-out: define significant genes in cohort A, measure the
    # fraction with same-direction effect in cohort B.
    loco = []
    for a in ids:
        for b in ids:
            if a == b:
                continue
            da, db = de[a], de[b]
            common = da.index.intersection(db.index)
            sig_a = da.loc[common]
            sig_a = sig_a[sig_a["fdr"] < fdr_thresh]
            if len(sig_a) == 0:
                continue
            same = (np.sign(sig_a["log2FC"]) ==
                    np.sign(db.loc[sig_a.index, "log2FC"]))
            loco.append({"discovery": a, "validation": b,
                         "n_sig": len(sig_a),
                         "same_direction_frac": float(same.mean())})
    loco = pd.DataFrame(loco)

    # Size-matched random-signature control: is the replicated signature's
    # cross-cohort concordance greater than random gene sets of the same size?
    merged = stage1["merged"]
    fc_cols = [c for c in merged.columns if c.startswith("log2FC_")]
    obs_sig = merged.index[merged["replicated_up"] | merged["replicated_down"]]
    n = max(len(obs_sig), 1)
    all_genes = merged.index.to_numpy()

    def concordance(idx):
        sub = merged.loc[idx, fc_cols]
        return float((sub.gt(0).all(axis=1) | sub.lt(0).all(axis=1)).mean())
    obs_conc = concordance(obs_sig) if len(obs_sig) else np.nan
    null = np.array([concordance(rng.choice(all_genes, n, replace=False))
                     for _ in range(n_perm)])
    emp_p = (np.sum(null >= obs_conc) + 1) / (n_perm + 1)
    emp_p_str = "<0.001" if emp_p < 1 / (n_perm + 1) + 1e-12 else f"{emp_p:.3g}"

    return {"loco": loco,
            "random_control": {"observed_concordance": obs_conc,
                               "null_mean": float(null.mean()),
                               "null_std": float(null.std()),
                               "n_perm": n_perm,
                               "empirical_p": emp_p, "empirical_p_str": emp_p_str,
                               "null": null}}
