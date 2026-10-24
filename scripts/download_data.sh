#!/usr/bin/env bash
# One-shot downloader for all real inputs of the IPF framework.
# Run from the project root:   bash scripts/download_data.sh
# Requires: wget, python (with the ipffx env active for the prepare/build steps).
#
# Downloads land in ./data/. Bulk GEO series-matrix + supplementary files are
# fetched from the stable NCBI GEO FTP tree; the atlas and LINCS have their own
# steps below. Big files — expect several GB and a while for the atlas.
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p data
GEO=https://ftp.ncbi.nlm.nih.gov/geo/series

dl_series () {   # $1 = GSE id, $2 = ftp folder bucket (e.g. GSE150nnn)
  local id=$1 bucket=$2
  echo ">>> $id : series matrix + supplementary"
  mkdir -p "data/$id"
  # series matrix (has phenotype for all; has expression for microarray)
  wget -c -q --show-progress -P "data/$id" \
    "$GEO/$bucket/$id/matrix/${id}_series_matrix.txt.gz" || \
    echo "   (no single series_matrix — multi-platform; prepare_bulk.py will handle)"
  # all supplementary files (RNA-seq count tables live here)
  wget -c -q --show-progress -r -np -nd -R "index.html*" \
    -A "*.txt.gz,*.csv.gz,*.tsv.gz,*.tar,*.gz,*.xlsx" \
    -P "data/$id" "$GEO/$bucket/$id/suppl/" || true
}

# ---- Bulk cohorts ----------------------------------------------------------
dl_series GSE47460  GSE47nnn      # LGRC microarray  (discovery)
dl_series GSE150910 GSE150nnn     # RNA-seq          (discovery)
dl_series GSE134692 GSE134nnn     # RNA-seq          (external bulk)

# ---- Single-cell atlas (GSE135893, Habermann 2020) -------------------------
echo ">>> GSE135893 : atlas supplementary (raw matrix + metadata)"
dl_series GSE135893 GSE135nnn
cat <<'MSG'
   Atlas note: the EASIEST route to an .h5ad is the IPF Cell Atlas:
     https://www.ipfcellatlas.com/   (Vanderbilt/Habermann dataset -> download .h5ad)
   Put the file at data/GSE135893/ipf_atlas.h5ad
   Otherwise run:  python scripts/build_atlas_h5ad.py   (assembles from GEO suppl)
MSG

# ---- (optional) external snRNA atlas GSE136831 (Adams 2020) ----------------
# dl_series GSE136831 GSE136nnn

# ---- LINCS L1000 perturbational signatures (public Enrichr GMTs) -----------
echo ">>> LINCS L1000 chem-perturbation signatures (Enrichr, no login)"
python scripts/build_lincs.py --out data/lincs

echo
echo "Download complete. Next:"
echo "  1) python scripts/prepare_bulk.py        # -> data/<GSE>/expr.csv + meta.csv"
echo "  2) (atlas) place ipf_atlas.h5ad or run build_atlas_h5ad.py"
echo "  3) bash scripts/run_real.sh              # full analysis + figures + tables"
