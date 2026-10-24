# Reproducibility test

## Environment tested

- OS: Linux (Ubuntu 24.04 LTS)
- Hardware: NVIDIA RTX 5070 (CUDA 12.4), 32 GB RAM
- Python: 3.10 (from `environment.yml`)
- Key pinned packages: `pydeseq2==0.4.12`, `anndata<0.11`, `numpy<2`

## Commands

```bash
conda env create -f environment.yml
conda activate ipffx
make demo
make download && make prepare && make atlas && make lincs && make real
```

## Expected outputs

All files listed in `metadata/DATA_MANIFEST.csv` with their SHA-256 digests in
`metadata/CHECKSUMS.sha256`.

## Tolerances

- DESeq2 outputs are reproducible to the precision allowed by `pydeseq2 0.4.12`
  on a fixed BLAS backend. Expected agreement: relative tolerance 1e-4 on
  log2FC estimates across BLAS backends.
- CIBERSORT nu-SVR deconvolution may differ by ±0.01 in proportion estimates
  across scikit-learn minor versions.
- LIANA rank-aggregate is deterministic given the random seed.

## Observed status

Reproduced successfully in the tested environment on 2026-10-24. All
documented tables regenerated with checksums matching the manifest. The
run_summary.json for the primary pipeline is at `results_real/run_summary.json`.
