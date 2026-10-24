# Data dictionary

## results_real/tables/S_stage1_replicated_signature.csv

| Column | Type | Meaning |
|---|---|---|
| gene (Unnamed: 0) | string | HGNC symbol |
| log2FC_GSE47460 | float | DESeq2 log2 fold-change (IPF vs control), GSE47460 |
| fdr_GSE47460 | float | BH-adjusted p-value, GSE47460 |
| log2FC_GSE150910 | float | DESeq2 log2 fold-change (IPF vs control), GSE150910 |
| fdr_GSE150910 | float | BH-adjusted p-value, GSE150910 |
| mean_log2FC | float | mean log2FC across both cohorts |
| replicated_up | bool | FDR<0.05 and log2FC>0 in both cohorts |
| replicated_down | bool | FDR<0.05 and log2FC<0 in both cohorts |
| same_direction | bool | log2FC sign agrees across cohorts |

## results_real/tables/S_stage2_pseudobulk_deg_counts.csv

| Column | Type | Meaning |
|---|---|---|
| cell_type | string | broad compartment (Epithelial/Immune/Endothelial/Mesenchymal) |
| up | int | BH-FDR<0.05 up-regulated gene count |
| down | int | BH-FDR<0.05 down-regulated gene count |

## results_real/tables/S_stage3_module_shift.csv

Wide-form matrix: rows = module name, columns = compartment, cells = donor-level
pseudobulk module-score log2FC (IPF − control).

## results_real/tables/S_stage4_ligand_receptor.csv

| Column | Type | Meaning |
|---|---|---|
| axis | string | ligand-receptor axis label (e.g. CCL2-CCR2) |
| route | string | signalling route (chemokine, TGF-beta, cytokine, PDGF, WNT, matrix-integrin) |
| ligand_log2FC | float | ligand bulk log2FC (discovery cohort) |
| receptor_log2FC | float | receptor bulk log2FC |
| sender | string | sender compartment(s), comma-separated |
| receiver | string | receiver compartment(s), comma-separated |
| evidence_score | float | composite score (ligand log2FC + prior support) |

## results_real/tables/S_stage5_external_bulk.csv

| Column | Type | Meaning |
|---|---|---|
| module (Unnamed: 0) | string | curated module |
| IPF | float | mean module score in IPF samples (GSE134692) |
| control | float | mean module score in control samples |
| IPF_minus_control | float | difference |

## results_real/tables/S_stage6_perturbation.csv

| Column | Type | Meaning |
|---|---|---|
| compound | string | LINCS compound name |
| bulk | float | reversal score vs bulk signature |
| Epithelial | float | reversal score vs epithelial pseudobulk signature |
| recurrence | int | number of signatures the compound reverses |
| multi_signature_support | float | aggregate reversal |
| module_support | float | module-level support |
| known_target | int | 1 if compound hits a signature target gene |
| safety_penalty | float | penalty (toxicity / environmental agents / off-label) |
| total_score | float | composite translational score |
| category | string | prioritized / supporting / cautionary / excluded_safety |

## results_real/tables/S_stage7_loco.csv

| Column | Type | Meaning |
|---|---|---|
| discovery | string | cohort held in for discovery |
| validation | string | cohort used for validation |
| n_sig | int | discovery-significant gene count |
| same_direction_frac | float | fraction of n_sig that preserve direction in validation |

## results_real/tables/S7e_jackknife.csv

| Column | Type | Meaning |
|---|---|---|
| compound | string | LINCS compound name |
| median_rank_jackknife | float | median rank across jackknife iterations |
| top10_frequency_pct | float | % of jackknife iterations where the compound was in the top 10 |
