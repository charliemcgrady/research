# Strict-cohort rerun — TCGA taxonomy excluding the 28 WGS-reclassified tumors

The TCGA "oncogene-negative" set carries a soft boundary: 28 of the 118 RPA-negative tumors were assigned
a driver alteration on WGS re-analysis (`Current_RPA_Alteration`). This rerun (`build_tcga_taxonomy_strict.py`)
**excludes those 28** to give a WGS-strict oncogene-negative cohort (**n = 89**), and moves the 28 into the
oncogene-positive control group. Per-tumor posteriors are unchanged (they are computed against cohort-wide
references); what changes is the population and the control.

## Result 1 — the taxonomy is robust to the strict definition

| Class | Full (n=117) | Strict (n=89) |
|---|---|---|
| Tumor-suppressor convergence (C2) | 35 (29.9%) | **28 (31.5%)** |
| High-plasticity / entropy (C7) | 32 (27.4%) | 22 (24.7%) |
| Microenvironment (C6) | 23 (19.7%) | 18 (20.2%) |
| Hidden canonical (C1) | 20 (17.1%) | 15 (16.9%) |
| Ligand / autocrine (C5) | 5 (4.3%) | 4 (4.5%) |
| Lineage-conditioned (C4) | 1 | 1 |
| Unknown (C8) | 1 | 1 |

Every proportion moves ≤3 points; the ranking is identical. Mean posterior masses are essentially
unchanged (C2 0.259→0.262, C7 0.223→0.217, C1 0.171→0.173). **The mechanistic taxonomy does not depend on
the reclassified tumors.** Data: `strict_vs_full_comparison.csv`, `tumor_evidence_TCGA_strict.csv`.

## Result 2 — the 28 excluded tumors validate the engine's limits *and* its hits

Where the engine had placed the 28 (in the full run) vs the WGS driver that reclassified them:

| WGS driver (n) | Engine had called them | Interpretation |
|---|---|---|
| **KRAS point mutation (20)** | C7 entropy (10), C2 (3), C6 (4), C1 (2), C5 (1) | **Invisible by design** — KRAS activation is a point mutation, not an RTK expression/CNA event, so the bulk engine cannot see it. These scattered mostly into the *unexplained / entropy / microenvironment* bins. |
| **Amplification: MAPK1×3, EGFR×1, ARAF×1** | C1 hidden-canonical (EGFR_amp, 2× MAPK1_amp), C2 (MAPK1_amp, ARAF_amp) | **Partially caught** — CNA-visible reclassifications were flagged by the cell-intrinsic classes C1/C2. |
| **NRG1 fusion (1)** | **C5 ligand** | **Direct hit** — the engine flagged elevated NRG1 (ligand class) on the exact tumor that carries an NRG1 fusion. |
| **Suppressor deletion: RASA1_del, NF1_del** | C2 suppressor (RASA1), C6 (NF1) | RASA1 correctly landed in suppressor-convergence. |

Two honest lessons:

1. **KRAS point mutants are the taxonomy's blind spot.** 20 of 28 reclassified tumors were KRAS-mutant, and
   the bulk RNA/CNA engine had no way to see them — they hid disproportionately in the **C7 entropy** and
   **C6 microenvironment** bins. This is the single strongest argument that the "unexplained / high-
   plasticity" tail (Deliverable #4) genuinely conceals detectable, even *targetable* (KRAS G12C) drivers,
   and that WGS is the missing layer (knowledge gap G1). It also means C7/C6 prevalence in any bulk
   oncogene-negative analysis is inflated by cryptic point-mutation drivers.
2. **The cell-intrinsic classes (C1, C2, C5) point at the right tumors when the driver is CNA/expression/
   fusion-visible.** EGFR/MAPK1 amplifications → C1; NRG1 fusion → C5; RASA1 deletion → C2. This is a real
   positive control: when a hidden driver leaves a transcriptomic/copy-number footprint, the taxonomy finds
   it.

## Positive control (strict)

With the 28 added, the oncogene-positive control grows to n=407; mean C1 mass stays higher in positives
(0.228) than in the strict negatives (0.173), preserving the control's direction.

## Bottom line

Excluding the WGS-reclassified tumors leaves the taxonomy structurally unchanged (C2-led, ~31%), so the
mission conclusions hold on the strict cohort. The exercise's real value is diagnostic: it shows the
taxonomy **cannot see KRAS point mutations** (they masquerade as entropy/microenvironment) but **does
correctly flag copy-number-, expression-, and fusion-visible hidden drivers** into C1/C2/C5 — sharpening
both the "sequence the unexplained tail with WGS" recommendation and the confidence in the cell-intrinsic
classes.
