# TCGA CIN × driver-status interaction — the gating test

**Question:** does driver status genuinely *modify* the CIN→LAG3 association (i.e. is the effect
driver-negative-specific), or is it a general LUAD phenomenon? Formal interaction model on the full
TCGA-LUAD cohort, across both CIN metrics × both infiltration covariates, + subgroup slopes, +
TIGIT comparator.

## Result: WEAK / ABSENT interaction — the effect is NOT driver-specific

**LAG3 interaction term (CIN_z × driver_positive_RPA), all 4 models:**

| CIN metric | infiltration cov | interaction coef | p |
|---|---|---|---|
| cna_mean_abs_log2 | tcell_infiltration | −0.077 | 0.44 |
| cna_mean_abs_log2 | cibersort_CD8_abs | −0.117 | 0.30 |
| aneuploidy_score | tcell_infiltration | +0.061 | 0.56 |
| aneuploidy_score | cibersort_CD8_abs | −0.015 | 0.90 |

All non-significant (p 0.30–0.90). **No evidence that driver status modifies the slope.**

**Subgroup slopes (CIN→LAG3), driver-negative vs driver-positive — they do not differ materially:**

| model | driver-negative β (p) | driver-positive β (p) |
|---|---|---|
| cna × tcell_infil | 0.332 (0.001) | 0.240 (0.0005) |
| cna × cibersort | 0.283 (0.011) | 0.053 (0.48) |
| aneuploidy × tcell_infil | 0.202 (0.029) | 0.213 (0.0002) |
| aneuploidy × cibersort | 0.125 (0.22) | 0.069 (0.27) |

The CIN→LAG3 slope is **positive in BOTH subgroups** (driver-negative numerically a touch larger but
not significantly so). The covariate-dependence noted in E1b persists (stronger with the expression
T-cell score, weaker with CIBERSORT-abs) — but that is orthogonal to the interaction question.

**TIGIT comparator:** interaction ns in all models (p 0.32–0.97); subgroup slopes null or only
significant (negative) in driver-positive for one covariate → confirms the positive per-cell signal
is **LAG3-specific**, not a generic checkpoint effect.

## Verdict & branch taken
Driver status does **not** modify the CIN–LAG3 association → the effect is a **general LUAD**
phenomenon, present in oncogene-negative and oncogene-positive tumors alike. This mirrors the
original AutoDS run2 finding (interaction ns) but now with gold-standard metrics and both covariates.

**Consequences (per plan):**
1. **Proceed to the general all-LUAD LuCA single-cell analysis** (43 eligible primary tumors) —
   framed as validation of a *general* CIN–CD8-LAG3 mechanism, TIGIT as prespecified comparator.
2. **Keep the driver-negative finding as a bulk-only subgroup observation** — real within that
   subgroup, but not specific to it.
3. **Do NOT pursue the driver-stratified multiplex-IHC route** — that was gated on *strong* evidence
   of driver-negative specificity, which the interaction test does not provide.
