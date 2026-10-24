"""Streamlit GUI for the IPF transcriptomic framework.

Launch:  streamlit run app/streamlit_app.py
Then open the URL it prints (usually http://localhost:8501).
"""
from __future__ import annotations
import sys
from pathlib import Path
import io

import pandas as pd
import streamlit as st

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from ipffx import data, pipeline, config  # noqa: E402

st.set_page_config(page_title="IPF Transcriptomic Framework", layout="wide",
                   page_icon="🫁")

STAGE_INFO = [
    ("1. Cross-cohort bulk", "Reproducible IPF signature across >=2 bulk cohorts."),
    ("2. Pseudobulk localization", "Donor-level cell-type DEGs (which compartment)."),
    ("3. Module scoring", "Curated epithelial / fibrotic / TGF-b module shifts."),
    ("4. Ligand-receptor", "Focused matrix-integrin / cytokine axis prioritization."),
    ("5. External validation", "Bulk + snRNA external support, etiology heterogeneity."),
    ("6. Perturbational (LINCS)", "Safety-aware drug-repositioning hypotheses."),
    ("7. Robustness", "Leave-one-cohort-out + size-matched random controls."),
]

st.title("🫁 IPF Transcriptomic Framework")
st.caption("Reproducible cross-cohort bulk → donor-level pseudobulk localization → "
           "curated module scoring → focused ligand-receptor → external validation → "
           "LINCS safety-aware perturbational prioritization → robustness. "
           "All perturbational/module results are computational hypotheses, not "
           "therapeutic evidence.")

# --------------------------------------------------------------------------- #
# Sidebar: data source
# --------------------------------------------------------------------------- #
with st.sidebar:
    st.header("1  Data source")
    mode = st.radio("Choose input", ["Demo (synthetic, no downloads)", "Real data (upload files)"])
    n_perm = st.slider("Random-control permutations", 100, 5000, 1000, step=100)
    outdir = st.text_input("Output directory", "results_gui")

    bundle = None
    if mode.startswith("Demo"):
        seed = st.number_input("Demo random seed", 0, 9999, 7)
        if st.button("Load demo data", type="primary"):
            st.session_state["bundle"] = data.make_demo_bundle(seed=int(seed))
            st.session_state.pop("results", None)
            st.success("Demo data loaded.")
    else:
        st.markdown("**Bulk discovery cohorts** (need ≥2). Each: expr (genes×samples) "
                    "+ metadata with a group column.")
        st.info("See scripts/DOWNLOAD_DATA.md for where to get GSE47460, GSE150910, "
                "GSE135893, GSE136831, GSE134692 and LINCS.")
        group_col = st.text_input("Metadata group column", "diagnosis")
        case = st.text_input("Case label (IPF)", "IPF")
        control = st.text_input("Control label", "control")
        cohorts = {}
        for i in range(1, 4):
            with st.expander(f"Bulk cohort {i}" + (" (external)" if i == 3 else "")):
                cid = st.text_input(f"Cohort {i} id", ["GSE47460", "GSE150910", "GSE134692"][i-1],
                                    key=f"cid{i}")
                ef = st.file_uploader(f"expr {i} (csv/tsv)", key=f"ef{i}")
                mf = st.file_uploader(f"meta {i} (csv/tsv)", key=f"mf{i}")
                if ef and mf:
                    ep = Path(outdir) / f"_upload_expr{i}.csv"
                    mp = Path(outdir) / f"_upload_meta{i}.csv"
                    Path(outdir).mkdir(parents=True, exist_ok=True)
                    ep.write_bytes(ef.getvalue()); mp.write_bytes(mf.getvalue())
                    cohorts[i] = (cid, str(ep), str(mp))
        h5 = st.file_uploader("Atlas .h5ad (pseudobulk)", type=["h5ad"])
        donor = st.text_input("h5ad donor column", "donor")
        ctcol = st.text_input("h5ad cell-type column", "cell_type")
        condcol = st.text_input("h5ad condition column", "disease")
        lincs_f = st.file_uploader("LINCS signature matrix (csv: compounds×genes)")
        if st.button("Load real data", type="primary"):
            try:
                b = {"bulk": {}}
                for i in (1, 2):
                    cid, ep, mp = cohorts[i]
                    b["bulk"][cid] = data.load_bulk_cohort(cid, ep, mp, group_col, case, control)
                cid, ep, mp = cohorts[3]
                b["bulk_external"] = data.load_bulk_cohort(cid, ep, mp, group_col, case, control)
                hp = Path(outdir) / "_atlas.h5ad"; hp.write_bytes(h5.getvalue())
                b["pseudobulk"] = data.pseudobulk_from_h5ad(str(hp), donor, ctcol, condcol, case, control)
                mods = config.load_modules()
                rows = []
                for ct, res in b["pseudobulk"]["de"].items():
                    for mid, m in mods.items():
                        present = [g for g in m["genes"] if g in res.index]
                        rows.append({"cell_type": ct, "module": mid,
                                     "shift": res.loc[present, "log2FC"].mean() if present else 0.0})
                b["snrna_external"] = pd.DataFrame(rows)
                lp = Path(outdir) / "_lincs.csv"; lp.write_bytes(lincs_f.getvalue())
                b["lincs"] = data.load_lincs_matrix(str(lp))
                b["drug_meta"] = pd.DataFrame(columns=["target_class", "risk_class"],
                                              index=b["lincs"].index).fillna("")
                b["is_demo"] = False
                st.session_state["bundle"] = b
                st.session_state.pop("results", None)
                st.success("Real data loaded.")
            except Exception as e:
                st.error(f"Could not load: {e}")

# --------------------------------------------------------------------------- #
# Run
# --------------------------------------------------------------------------- #
bundle = st.session_state.get("bundle")
if bundle is None:
    st.info("← Load data in the sidebar to begin. Start with **Demo** to see the whole "
            "workflow run instantly.")
    st.subheader("The seven stages")
    for name, desc in STAGE_INFO:
        st.markdown(f"**{name}** — {desc}")
    st.stop()

if bundle.get("is_demo"):
    st.warning("DEMO MODE — synthetic data. Results illustrate the workflow only; "
               "they are not real IPF findings.")

if st.button("▶ Run full pipeline", type="primary"):
    prog = st.progress(0.0, text="Starting…")
    def cb(name, frac):
        prog.progress(min(frac, 1.0), text=name)
    st.session_state["results"] = pipeline.run_all(bundle, outdir=outdir, n_perm=n_perm,
                                                    progress=cb)
    prog.progress(1.0, text="Done")

res = st.session_state.get("results")
if res is None:
    st.stop()

# --------------------------------------------------------------------------- #
# Results tabs
# --------------------------------------------------------------------------- #
s = res["summary"]
c1, c2, c3, c4 = st.columns(4)
c1.metric("Replicated up genes", s["stage1_summary"]["n_replicated_up"])
c2.metric("Replicated down genes", s["stage1_summary"]["n_replicated_down"])
c3.metric("Same-direction genes", s["stage1_summary"]["n_same_direction"])
c4.metric("Random-control P", s["robustness"]["empirical_p_str"])

tabs = st.tabs(["Workflow", "1 Bulk", "2 Pseudobulk", "3 Modules",
                "4 Ligand-receptor", "5 External", "6 Perturbation", "7 Robustness",
                "Downloads"])

with tabs[0]:
    st.image(res["fig1"])
    st.dataframe(res["table1"], use_container_width=True)

with tabs[1]:
    st.image(res["fig2"])
    st.dataframe(res["stage1"]["merged"].head(200), use_container_width=True)

with tabs[2]:
    st.image(res["fig2"])
    st.dataframe(res["stage2"]["deg_counts"], use_container_width=True)
    st.caption("Localization of the replicated bulk signature to each compartment:")
    st.dataframe(res["stage2"]["localization"], use_container_width=True)

with tabs[3]:
    col = st.columns(2)
    col[0].image(res["fig3"]); col[1].image(res["fig5"])
    st.dataframe(res["stage3"]["module_shift"], use_container_width=True)

with tabs[4]:
    st.image(res["fig4"])
    st.dataframe(res["stage4"]["axes"], use_container_width=True)

with tabs[5]:
    st.image(res["fig6"])
    st.dataframe(res["stage5"]["external_bulk"], use_container_width=True)

with tabs[6]:
    st.image(res["fig7"])
    pert = res["stage6"]["perturbation"]
    st.dataframe(pert, use_container_width=True)
    st.caption(f"Positive-control ranks (approved antifibrotics): "
               f"{res['stage6']['positive_control_ranks']}")
    st.info("Prioritized compounds are computational perturbational hypotheses "
            "requiring experimental validation — never evidence of efficacy.")

with tabs[7]:
    rc = s["robustness"]
    st.write(rc)
    st.dataframe(res["stage7"]["loco"], use_container_width=True)

with tabs[8]:
    st.markdown(f"All figures (PNG, 300 dpi) and tables (CSV) were written to "
                f"`{res['outdir']}/`.")
    for key in ["fig1", "fig2", "fig3", "fig4", "fig5", "fig6", "fig7"]:
        p = Path(res[key])
        st.download_button(f"Download {p.name}", p.read_bytes(), file_name=p.name)
    tdir = Path(res["outdir"]) / "tables"
    for csv in sorted(tdir.glob("*.csv")):
        st.download_button(f"Download {csv.name}", csv.read_bytes(), file_name=csv.name)
