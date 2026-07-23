# E1b — disentangling CIN metric, purity, and infiltration covariate

Decisive follow-up to E1. Same driver-negative TCGA-LUAD data, but testing all three CIN metrics
(`cna_mean_abs_log2` = run6 amplitude proxy; `cin_frac_altered_0.3` = FGA proxy; `aneuploidy_score`
= Taylor arm-level) × ±purity, under two infiltration covariates. **The finding splits cleanly into
two arms with opposite fates.**

## Arm 1 — "CIN reduces T-cell infiltration" (run6 β=−4.10, p<1e-4) → **purity artifact**

PART B: CIN → tcell_infiltration_score, ± purity, and each metric's correlation with purity:

| CIN metric | coef (no purity) | coef (+ purity) | corr with purity |
|---|---|---|---|
| cna_mean_abs_log2 | **−0.501, p<0.001** | −0.150, **p=0.092 (ns)** | **r=0.54** |
| cin_frac_altered_0.3 | −0.490, p<0.001 | −0.167, p=0.046 | r=0.54 |
| aneuploidy_score | −0.032, **p=0.74 (ns)** | −0.137, p=0.058 | r=−0.10 |

**Verdict:** the infiltration-exclusion effect is **substantially purity-driven.** The amplitude/FGA
proxies correlate with tumor purity at r=0.54 (high-purity tumors have less stroma/immune); adjusting
for purity collapses the effect (p→0.09), and the purity-robust Aneuploidy Score (r=−0.10 with purity)
shows **no infiltration effect at all** (p=0.74). run6's headline exclusion result does not survive.

## Arm 2 — "per-cell LAG3 rises with CIN" → **survives purity, LAG3-specific, but covariate-dependent**

PART A (LAG3 ~ CIN + **tcell_infiltration_score** + tmb):

| CIN metric | coef (no purity) | coef (+ purity) |
|---|---|---|
| cna_mean_abs_log2 | +0.328, **p=0.001** | +0.364, **p=0.001** |
| cin_frac_altered_0.3 | +0.227, **p=0.023** | +0.245, **p=0.021** |
| aneuploidy_score | +0.187, **p=0.036** | +0.201, **p=0.029** |

- **Positive and significant for all three CIN metrics, and purity does NOT abolish it** (unlike Arm 1).
- **TIGIT (Part A2): null for every metric (p>0.3)** → LAG3-specificity holds.

PART C2 sensitivity (LAG3 ~ CIN + **cibersort_CD8_absolute** + tmb) — swapping the infiltration covariate:

| CIN metric | coef (no purity) | coef (+ purity) |
|---|---|---|
| cna_mean_abs_log2 | +0.136, p=0.228 | +0.311, **p=0.011** |
| cin_frac_altered_0.3 | +0.062, p=0.568 | +0.201, p=0.083 |
| aneuploidy_score | +0.164, p=0.103 | +0.124, p=0.215 |

**Verdict:** the per-cell LAG3 signal is real but **fragile and covariate-dependent** — robust and
significant when adjusting for the expression-based T-cell score (Part A), but only the amplitude
metric (+purity) reaches significance when adjusting for the independent CIBERSORT-absolute estimate
(Part C2). This resolves the E1↔run6 discrepancy: E1 used cibersort_CD8_absolute (mostly null); run6/E1b
Part A used the expression T-cell score (significant).

> ⚠️ **Circularity caveat:** `tcell_infiltration_score` is built from CD8A/CD3D/… *expression* in the
> same RNA sample as LAG3, so adjusting LAG3 for it is expression-on-expression and can leave residual
> co-expression. The more independent CIBERSORT-absolute estimate gives a weaker signal. This is exactly
> the bulk-deconvolution ambiguity that **single-cell data removes.**

## Metric relationships (Part C1)
cna_mean_abs_log2 ↔ cin_frac_altered_0.3: r=0.90; cna ↔ aneuploidy_score: r=0.61; cin_frac ↔ AS: r=0.52.
The two proxies are near-identical; both diverge moderately from arm-level aneuploidy.

## Overall verdict & implication for the plan
- **Infiltration-exclusion arm: retract** — a tumor-purity artifact.
- **Per-cell LAG3 arm: alive but unresolved** — purity-robust and LAG3-specific, but its significance
  hinges on how infiltration is modeled (bulk deconvolution can't cleanly separate per-cell level from
  composition). Not the confident result run6 implied, not the dead result E1 alone implied.
- **This now *justifies* the single-cell tier (E4–E6)** rather than gating it out: single-cell measures
  per-cell LAG3 in CD8 T cells directly, eliminating the deconvolution/covariate ambiguity that is the
  sole remaining source of doubt. That is the definitive arbiter.
