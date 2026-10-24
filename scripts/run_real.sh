#!/usr/bin/env bash
# Full real-data analysis on prepared public inputs.
# Prerequisites: scripts/download_data.sh, prepare_bulk.py, build_atlas_h5ad.py.
set -euo pipefail
cd "$(dirname "$0")/.."

OUT=${1:-results_real}

echo ">>> Running 7-stage pipeline -> $OUT"
python scripts/run_pipeline.py   --bulk GSE47460=data/GSE47460/expr.csv,data/GSE47460/meta.csv,diagnosis,IPF,control   --bulk GSE150910=data/GSE150910/expr.csv,data/GSE150910/meta.csv,diagnosis,IPF,control   --bulk-external GSE134692=data/GSE134692/expr.csv,data/GSE134692/meta.csv,diagnosis,IPF,control   --h5ad data/GSE135893/ipf_atlas.h5ad --donor donor --celltype cell_type --condition disease   --case IPF --control control   --lincs data/lincs/lincs_signatures.csv --drug-meta data/lincs/drug_meta.csv   --out "$OUT" --n-perm 2000

echo
echo "DONE."
echo "  Figures : $OUT/figures/"
echo "  Tables  : $OUT/tables/"
echo "  Summary : $OUT/run_summary.json"
