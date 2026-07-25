# Mechanistic Taxonomy of Oncogene-Negative LUAD — v1 (TCGA)

First evidence-based, per-tumor mechanistic assignment across the mission's 8 classes.
Engine: soft, evidence-weighted posteriors (NOT p-values) — multiple weak sources combined per tumor.
Cohort: TCGA-LUAD RPA-negative, n=117 (data-complete).

## Class assignment (best-explaining class per tumor)
| Class | n (33%…) | mean mass |
|---|---|---|
| C2 Tumor-suppressor convergence (NF1/RASA1/PTEN/STK11/KEAP1 + MAPK/PI3K output) | 39 | 0.30 |
| C6 Microenvironment-dependent (leukocyte/stromal) | 24 | 0.18 |
| C1 Hidden canonical (RTK expr-outlier / EGFR-ERBB2-MET amp / MET-ex14) | 24 | 0.18 |
| C7 High state-entropy (residual: high lineage entropy, no lesion class) | 13 | 0.13 |
| C8 Unknown | 12 | 0.06 |
| C5 Ligand-driven (NRG1/HGF/EGF-family/FGF/IGF expr-outlier) | 5 | 0.07 |
| C4 Lineage-conditioned (mucinous/basal + pathway output) | 0 best | 0.08 |
| C3 Protein-state | — | not in TCGA (needs CPTAC) |

## Positive control (Exp 30 / calibration) — PASSES
Same engine on driver-POSITIVE tumors (n=379): C1 hidden-canonical is the leading class (116),
and C1 mean mass is higher in driver-positive (0.25) than driver-negative (0.18). The engine
recovers canonical activation where it genuinely exists.

## Concordance with prior art
- C2 dominance matches the established tumor-suppressor-convergence view of oncogene-negative LUAD.
- C1 ≈ 20% ≈ Chen 2021's ~33% WGS-reclassification rate (we recover a portion from RNA/CN only).

## Caveats (v1)
- **Class 3 (protein-state) absent** — requires CPTAC phospho/protein (separate ~18-tumor cohort).
- Evidence weights are heuristic; per-class thresholds not yet literature-calibrated.
- Bulk data conflates malignant/stromal (C6 partly admixture; C4/C7 bulk-limited).
- Entropy is modeled as a residual class (only claims tumors no lesion-class explains) to avoid the
  bulk-admixture inflation seen in the uncalibrated first pass.

## Next (to complete Deliverables #1–6)
1. Add CPTAC (Class 3 protein/phospho outliers) on its ~18 oncogene-negative tumors.
2. Per-tumor evidence GRAPHS (Deliverable #2) with explicit evidence items + competing posteriors.
3. Per-mechanism evidence graph (#3): prevalence, reproducibility, support.
4. Rank remaining Unknown tumors (#4) — highest-value for future sequencing.
5. Failure report (#5) — every rejected hypothesis (incl. the whole prior CIN→LAG3 arc, and the
   blocked structural classes).
6. Knowledge gaps (#6): the experiment/dataset that would most reduce uncertainty per class.
