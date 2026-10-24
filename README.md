# IPF Transcriptomic Framework (`ipffx`)

A reproducible computational framework for **idiopathic pulmonary fibrosis (IPF)** transcriptomic analysis with a Streamlit GUI.

The workflow integrates:

1. Cross-cohort bulk differential expression and directional concordance
2. Donor-level pseudobulk localization across broad tissue compartments
3. Configurable module scoring
4. Focused ligand-receptor prioritization
5. External bulk and single-nucleus corroboration
6. LINCS/L1000 transcriptomic-reversal hypothesis generation with safety-aware scoring
7. Robustness analyses including leave-one-cohort-out and random-control analyses

All perturbational and module-level outputs are **computational hypotheses** and do not establish therapeutic efficacy or pathway causality.

## Quick start

### Conda

```bash
conda env create -f environment.yml
conda activate ipffx
python scripts/run_pipeline.py --demo --out results_demo
```

### pip / virtual environment

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python scripts/run_pipeline.py --demo --out results_demo
```

The demo runs the complete workflow on synthetic data and writes figures, tables, and a machine-readable summary.

### GUI

```bash
./run_gui.sh
```

The GUI provides synthetic-demo and real-data workflows.

## Real-data workflow

Real input data are not redistributed. Download the public datasets and prepare them as described in `scripts/DOWNLOAD_DATA.md`.

```bash
bash scripts/download_data.sh
python scripts/prepare_bulk.py
python scripts/build_atlas_h5ad.py --h5ad /path/to/atlas.h5ad
bash scripts/run_real.sh
```

The real-data runner expects prepared inputs under `data/` and writes the analysis to `results_real/`.

## Configuration

The workflow is configuration-driven:

- `config/datasets.yaml` — dataset roles and metadata
- `config/modules.yaml` — curated gene modules
- `config/ligand_receptor.yaml` — ligand-receptor axes
- `config/drug_safety.yaml` — target credits and safety penalties

## Repository structure

```text
ipf_framework/
├── README.md
├── CITATION.cff
├── LICENSE
├── CHANGELOG.md
├── environment.yml
├── requirements.txt
├── Makefile
├── run_gui.sh
├── config/
├── ipffx/
├── app/
├── scripts/
├── results_real/
└── metadata/
```

## Released results

`results_real/` contains frozen analysis tables generated from the public-data workflow. These files are included as reproducibility artifacts; the source datasets themselves are not bundled.

Key outputs include:

- `S_stage1_replicated_signature.csv` — cross-cohort gene signature
- `S_stage2_*` — pseudobulk localization and related analyses
- `S_stage3_module_shift.csv` — donor-level module shifts
- `S_stage4_ligand_receptor.csv` — focused ligand-receptor prioritization
- `S_stage5_*` — external corroboration
- `S_stage6_perturbation.csv` — LINCS/L1000 prioritization
- `S_stage7_loco.csv` — leave-one-cohort-out robustness
- `S7e_jackknife.csv` — perturbational ranking stability

See `metadata/DATA_DICTIONARY.md` for output definitions.

## Reproducibility

The release records:

- source-data provenance
- environment specifications
- configuration files
- random seeds
- output manifests
- SHA-256 checksums

See `metadata/SOURCE_PROVENANCE.md`, `metadata/REPRODUCIBILITY_MANIFEST.json`, and `metadata/REPRODUCIBILITY_TEST.md`.

## Data sources

The workflow uses public resources including NCBI GEO cohorts, the IPF single-cell/single-nucleus atlases, and LINCS/L1000 perturbational signatures. Original data generators should be cited when using the corresponding datasets.

## License

MIT License. See `LICENSE`.
