# LuCA general all-LUAD single-cell analysis — the per-cell LAG3 signal does NOT replicate

Definitive single-cell test of the CIN → CD8-LAG3 mechanism in LUAD broadly (not driver-negative;
n=1 there). Removes the bulk-deconvolution ambiguity that limited E1b by measuring CD8 LAG3
directly in real CD8 cells and inferring CIN per tumor from malignant cells.

## Design (as prespecified; tumor = independent unit)
- **Cohort:** LuCA core, **43 eligible primary LUAD tumors** across **9 datasets** (≥50 malignant &
  ≥20 CD8 cells/tumor). All-eligible incl. metastases = 66 (sensitivity).
- **CIN:** infercnvpy on 23,843 malignant ("Tumor cells") vs 15,000 normal-epithelial reference cells
  (AT1/AT2/Club/Ciliated from normal/normal-adjacent); per-tumor CIN = mean CNV magnitude over
  malignant cells.
- **CD8 readout:** patient/tumor-level **pseudobulk** mean log-expr over 25,232 CD8 T cells.
- **Model:** tumor-level OLS `pb_GENE ~ CIN_z + C(dataset)` (cohort-adjusted); TIGIT prespecified comparator.

## Result — null across the board

| Gene | coef (CIN_z) | p (cohort-adj) | p (no cohort adj) |
|---|---|---|---|
| **LAG3** | **+0.026** | **0.56** | 0.63 |
| TIGIT | −0.016 | 0.56 | 0.36 |
| IFNG | −0.033 | 0.39 | 0.27 |

The LAG3 point estimate is **≈0**, not merely non-significant — this is a null, not just underpowered
absence of significance. TIGIT (comparator) and IFN-γ likewise null.

## Interpretation
The bulk "per-cell LAG3 rises with CIN" signal (run6) **does not survive** when per-cell LAG3 is
measured directly in single cell. Combined with the earlier findings, the LAG3 lead fails a graded
validation:
- run6 (bulk, expression-based infiltration score): +2.81, p=0.001 — the original signal.
- E1 (bulk, independent CIBERSORT-abs infiltration): null.
- E1b (bulk): purity-robust & LAG3-specific but **covariate-dependent** (significant only with the
  expression-based score — which is collinear with LAG3).
- E7 interaction: **general, not driver-specific**.
- **LuCA single cell (direct CD8 LAG3, per-tumor CIN): null (p=0.56, coef≈0).**

The most parsimonious reading: the bulk per-cell signal was largely an artifact of adjusting LAG3 for
an expression-derived T-cell score built from the same RNA (residual co-expression), not a true
per-cell CIN→LAG3 relationship. Single cell — the method that removes that circularity — shows no effect.

## Caveats
- n=43 tumors limits power for a *small* true effect, but the near-zero point estimate argues against a
  meaningful missed effect.
- inferCNV CIN is computed on expression-defined malignant cells (mild circularity) and pooled across 9
  platforms (cohort-adjusted). A single normal-epithelial reference set was used.
- All-LUAD, not driver-negative (untestable in single cell; n=1).

## Bottom line for the whole LAG3 arc
No robust, cross-modal evidence supports a CIN→CD8-LAG3 mechanism in LUAD. The driver-negative bulk
observation stands only as a subgroup finding sensitive to CIN-metric and infiltration-covariate choice,
and it does not generalize or replicate at single-cell resolution. Validation did its job:
a promising AutoDS lead was stress-tested and ultimately not supported.
