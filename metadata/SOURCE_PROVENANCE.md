# Source Provenance

## External datasets (not redistributed; retrieved from source)

| ID | Type | Role | n_IPF | n_control | Source | Access |
|---|---|---|---|---|---|---|
| GSE150910 | bulk RNA-seq | discovery | 103 | 103 | NCBI GEO | public |
| GSE47460  | bulk array    | replication | 160 | 108 | NCBI GEO | public |
| GSE134692 | bulk RNA-seq | third cohort / external bulk | 46 | 50 | NCBI GEO | public |
| GSE135893 | scRNA-seq (snRNA) | donor-level pseudobulk | 12 | 10 | NCBI GEO | public |
| GSE132771 | scRNA-seq | multi-atlas replication | - | - | NCBI GEO | public |
| LINCS L1000 | perturbational | reversal screen | - | - | Enrichr `LINCS_L1000_Chem_Pert_Consensus_Sigs` | public |

Retrieval scripts: `scripts/download_data.sh`, `scripts/build_atlas_h5ad.py`,
`scripts/build_lincs.py`.

Processing commit: see `metadata/REPRODUCIBILITY_MANIFEST.json → release_version`
and the git tag on this release.

## No restricted data

No patient-level identifiers, no controlled-access datasets, no proprietary
third-party content are included in this repository. All inputs are retrieved
at run time from public GEO / Enrichr endpoints.
