#!/usr/bin/env python3
"""Command-line runner for the IPF framework.

Demo (no downloads):
    python scripts/run_pipeline.py --demo --out results_demo

Real data (edit config/datasets.yaml, then provide local file paths):
    python scripts/run_pipeline.py \
        --bulk GSE47460=expr1.csv,meta1.csv,diagnosis,IPF,control \
        --bulk GSE150910=expr2.csv,meta2.csv,diagnosis,IPF,control \
        --bulk-external GSE134692=exprX.csv,metaX.csv,diagnosis,IPF,control \
        --h5ad atlas.h5ad --donor donor --celltype cell_type --condition disease \
        --case IPF --control control \
        --lincs lincs_signatures.csv \
        --out results_real
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from ipffx import data, pipeline  # noqa: E402


def parse_bulk(spec):
    cid, rest = spec.split("=", 1)
    expr, meta, col, case, control = rest.split(",")
    return data.load_bulk_cohort(cid, expr, meta, col, case, control)


def main():
    ap = argparse.ArgumentParser(description="IPF transcriptomic framework")
    ap.add_argument("--demo", action="store_true", help="run on synthetic demo data")
    ap.add_argument("--bulk", action="append", default=[],
                    help="cohort spec id=expr,meta,groupcol,case,control (>=2 for discovery)")
    ap.add_argument("--bulk-external", help="external bulk cohort, same spec format")
    ap.add_argument("--h5ad", help="single-cell atlas .h5ad for pseudobulk")
    ap.add_argument("--donor", default="donor")
    ap.add_argument("--celltype", default="cell_type")
    ap.add_argument("--condition", default="disease")
    ap.add_argument("--case", default="IPF")
    ap.add_argument("--control", default="control")
    ap.add_argument("--lincs", help="LINCS/L1000 signature matrix csv (compounds x genes)")
    ap.add_argument("--drug-meta", help="csv: compound,target_class,risk_class")
    ap.add_argument("--out", default="results", help="output directory")
    ap.add_argument("--n-perm", type=int, default=1000)
    args = ap.parse_args()

    if args.demo:
        bundle = data.make_demo_bundle()
    else:
        if len(args.bulk) < 2:
            ap.error("Provide >=2 --bulk discovery cohorts, or use --demo.")
        bundle = {"bulk": {c["id"]: c for c in (parse_bulk(s) for s in args.bulk)}}
        if args.bulk_external:
            bundle["bulk_external"] = parse_bulk(args.bulk_external)
        else:
            ap.error("--bulk-external is required in real-data mode.")
        if not args.h5ad:
            ap.error("--h5ad atlas is required in real-data mode.")
        bundle["pseudobulk"] = data.pseudobulk_from_h5ad(
            args.h5ad, args.donor, args.celltype, args.condition, args.case, args.control)
        # snRNA external: reuse pseudobulk module shift unless a separate atlas given
        import pandas as pd
        from ipffx import config as cfg
        mods = cfg.load_modules()
        rows = []
        for ct, res in bundle["pseudobulk"]["de"].items():
            for mid, m in mods.items():
                present = [g for g in m["genes"] if g in res.index]
                rows.append({"cell_type": ct, "module": mid,
                             "shift": res.loc[present, "log2FC"].mean() if present else 0.0})
        bundle["snrna_external"] = pd.DataFrame(rows)
        if not args.lincs:
            ap.error("--lincs signature matrix is required in real-data mode.")
        bundle["lincs"] = data.load_lincs_matrix(args.lincs)
        if args.drug_meta:
            bundle["drug_meta"] = pd.read_csv(args.drug_meta, index_col=0)
        else:
            bundle["drug_meta"] = pd.DataFrame(
                columns=["target_class", "risk_class"],
                index=bundle["lincs"].index).fillna("")
        bundle["is_demo"] = False

    def prog(name, frac):
        print(f"[{frac*100:5.1f}%] {name}")

    res = pipeline.run_all(bundle, outdir=args.out, n_perm=args.n_perm, progress=prog)
    print("\nDone. Outputs in:", res["outdir"])
    print("Summary:", res["summary"]["stage1_summary"])
    if res["summary"]["top_prioritized"]:
        print("Prioritized (hypotheses only):", res["summary"]["top_prioritized"])


if __name__ == "__main__":
    main()
