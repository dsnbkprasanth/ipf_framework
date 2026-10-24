"""
ipffx — a reproducible computational transcriptomic framework for cell-type-
resolved aberrant-epithelial / myofibroblast (TGF-beta) remodeling in idiopathic
pulmonary fibrosis (IPF).

The package re-implements, for IPF, the seven-stage flow of the reference
heart-failure framework:

  1. Cross-cohort bulk differential expression + reproducible signature
  2. Donor-level cell-type pseudobulk localization
  3. Curated module construction and scoring
  4. Focused ligand-receptor prioritization
  5. External validation (bulk + snRNA-seq)
  6. LINCS / L1000 safety-aware perturbational prioritization
  7. Robustness (leave-one-cohort-out + size-matched random controls)

Nothing here claims drug efficacy or pathway causality. Perturbational and
module results are computational hypotheses requiring experimental validation.
"""

__version__ = "0.1.0"
__all__ = ["config", "data", "pipeline"]
