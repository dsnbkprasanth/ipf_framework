"""End-to-end orchestration of the 7-stage IPF framework."""
from __future__ import annotations
from pathlib import Path
import json
import pandas as pd

from . import config, figures
from . import stages as S


def build_table1(outdir: Path) -> pd.DataFrame:
    ds = config.load_datasets()["datasets"]
    rows = [{"Resource": d["id"], "Name": d["name"], "Data type": d["data_type"],
             "Analytical role": d["role"],
             "Description": " ".join(d["description"].split())} for d in ds]
    tab = pd.DataFrame(rows)
    tab.to_csv(outdir / "tables" / "Table1_datasets.csv", index=False)
    return tab


def run_all(bundle: dict, outdir: Path | str | None = None,
            n_perm: int = 1000, progress=None) -> dict:
    """Run every stage, emit figures + tables. Returns a results dict.

    progress : optional callable(step_name, fraction) for GUIs.
    """
    outdir = Path(outdir) if outdir else config.default_output_dir()
    (outdir / "figures").mkdir(parents=True, exist_ok=True)
    (outdir / "tables").mkdir(parents=True, exist_ok=True)
    demo = bundle.get("is_demo", False)

    def tick(name, frac):
        if progress:
            progress(name, frac)

    res = {"outdir": str(outdir), "is_demo": demo}

    tick("Table 1 / workflow", 0.05)
    res["table1"] = build_table1(outdir)
    res["fig1"] = figures.fig1_workflow(outdir / "figures", demo=demo)

    tick("Stage 1: cross-cohort bulk", 0.15)
    s1 = S.stage1_cross_cohort(bundle)
    res["stage1"] = s1
    s1["merged"].to_csv(outdir / "tables" / "S_stage1_replicated_signature.csv")

    tick("Stage 2: pseudobulk localization", 0.30)
    s2 = S.stage2_pseudobulk(bundle, s1)
    res["stage2"] = s2
    s2["deg_counts"].to_csv(outdir / "tables" / "S_stage2_pseudobulk_deg_counts.csv")
    s2["localization"].to_csv(outdir / "tables" / "S_stage2_localization.csv")

    tick("Stage 3: module scoring", 0.45)
    s3 = S.stage3_modules(bundle)
    res["stage3"] = s3
    s3["module_shift"].to_csv(outdir / "tables" / "S_stage3_module_shift.csv")

    tick("Stage 4: ligand-receptor", 0.55)
    s4 = S.stage4_ligand_receptor(bundle)
    res["stage4"] = s4
    s4["axes"].to_csv(outdir / "tables" / "S_stage4_ligand_receptor.csv", index=False)

    tick("Stage 5: external validation", 0.65)
    s5 = S.stage5_external(bundle)
    res["stage5"] = s5
    s5["external_bulk"].to_csv(outdir / "tables" / "S_stage5_external_bulk.csv")
    s5["snrna_external"].to_csv(outdir / "tables" / "S_stage5_snrna_external.csv")

    tick("Stage 6: perturbational prioritization", 0.80)
    s6 = S.stage6_perturbation(bundle, s1)
    res["stage6"] = s6
    s6["perturbation"].to_csv(outdir / "tables" / "S_stage6_perturbation.csv", index=False)

    tick("Stage 7: robustness", 0.90)
    s7 = S.stage7_robustness(bundle, s1, n_perm=n_perm)
    res["stage7"] = s7
    s7["loco"].to_csv(outdir / "tables" / "S_stage7_loco.csv", index=False)

    tick("Rendering figures", 0.95)
    res["fig2"] = figures.fig2_crosscohort(outdir / "figures", s1, s2, s7)
    res["fig3"] = figures.fig3_modules(outdir / "figures", s3)
    res["fig4"] = figures.fig4_ligand_receptor(outdir / "figures", s4)
    res["fig5"] = figures.fig5_tgfb(outdir / "figures", s3)
    res["fig6"] = figures.fig6_external(outdir / "figures", s3, s5, s7)
    res["fig7"] = figures.fig7_perturbation(outdir / "figures", s6)

    # machine-readable summary
    summary = {
        "is_demo": demo,
        "stage1_summary": s1["summary"],
        "pseudobulk_deg_counts": s2["deg_counts"].reset_index().to_dict("records"),
        "positive_control_ranks": s6["positive_control_ranks"],
        "robustness": {k: v for k, v in s7["random_control"].items() if k != "null"},
        "top_prioritized": s6["perturbation"]
            .query("category == 'prioritized'")["compound"].tolist(),
    }
    with open(outdir / "run_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2, default=str)
    res["summary"] = summary
    tick("Done", 1.0)
    return res
