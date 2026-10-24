# Where to get the real data

Nothing patient-derived is bundled with this package. Download the public
resources below and point the CLI / GUI at the local files. A synthetic demo
(`--demo`) lets you run the whole workflow with no downloads.

## Bulk discovery + external cohorts (NCBI GEO)

| Cohort    | Role                | Link |
|-----------|---------------------|------|
| GSE47460  | bulk discovery (LGRC microarray, large ILD/IPF) | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47460 |
| GSE150910 | bulk discovery (IPF/CHP/control RNA-seq)        | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150910 |
| GSE134692 | external bulk validation (IPF vs control)       | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134692 |

For each, download the processed expression matrix (series matrix or supplementary
counts) and build two files:

* `expr.csv`  — genes (rows) × samples (columns). Counts or normalized/log values;
  raw counts are auto log2(x+1)-transformed.
* `meta.csv`  — samples (rows) × metadata; must contain a group column (e.g.
  `diagnosis`) with your case label (e.g. `IPF`) and control label (e.g. `control`).

`GEOparse` (installed in the env) can automate this, e.g.:

```python
import GEOparse
gse = GEOparse.get_GEO("GSE150910", destdir="./data")
# assemble expr (pivot of gse.pivot_samples) and a metadata frame from gse.phenotype_data
```

## Single-cell / single-nucleus atlas (pseudobulk localization)

| Atlas     | Role | Link / DOI |
|-----------|------|-----------|
| GSE135893 | main donor-level pseudobulk (Habermann 2020) | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135893 — DOI 10.1126/sciadv.aba1972 |
| GSE136831 | external snRNA support (Adams 2020)          | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136831 — DOI 10.1126/sciadv.aba1983 |

Provide the atlas as an `.h5ad` (AnnData) with raw counts in `.X` and `.obs`
columns for donor, cell type and condition. Processed objects are also on the
IPF Cell Atlas (www.ipfcellatlas.com). The pipeline aggregates to donor×cell-type
pseudobulk and runs per-cell-type disease DE.

## Perturbational signatures (LINCS / L1000)

* L1000FWD / CLUE: https://maayanlab.cloud/L1000FWD/ and https://clue.io/
* Provide a matrix `lincs.csv` — compounds (rows) × genes (columns), values =
  signed differential expression / z-scores each compound induces.
* Optional `drug_meta.csv` — columns `compound,target_class,risk_class`
  (see `config/drug_safety.yaml` for the recognized classes) to enable
  target-credit and safety-penalty scoring.

## Gene identifiers

Use HGNC gene symbols throughout (the curated modules in `config/modules.yaml`
and the L-R axes use symbols). If your matrices use Ensembl IDs, map them to
symbols before running.
