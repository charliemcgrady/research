# Tier-1 DataVoyager results — the per-cell LAG3 finding does NOT survive gold-standard re-testing

Two DataVoyager experiments re-tested the run6 finding (in driver-negative LUAD, CIN reduces
T-cell infiltration but raises LAG3 per T cell) after replacing our crude proxies with
publication-grade metrics. **The effect largely collapses.**

## E1 — Gold-standard re-test (Taylor Aneuploidy Score + absolute CD8 + purity/WGD)
Driver-negative LUAD (n≈115). Standardized predictors; OLS adjusting for infiltration
(CIBERSORT absolute CD8), purity, WGD, TMB.

| Model | aneuploidy_score coef | p |
|---|---|---|
| **LAG3** (with purity) | **−0.027** | **0.85** |
| LAG3 (no purity) | −0.011 | 0.94 |
| **TIGIT** | −0.054 | 0.65 |
| **IFNG** | +0.209 | 0.39 |
| **aneuploidy ↔ infiltration** (Pearson) | **r = −0.026** | **0.79** |

- The per-cell **LAG3 effect vanishes**: run6 had β=+2.81, p=0.001; here β≈−0.03, p=0.85.
- **IFN-γ per cell** also null (was +3.74, p=0.015).
- **Purity is not the culprit** — adding/removing it barely moves LAG3.
- **The most telling result:** with the standard, purity-corrected Aneuploidy Score, CIN has
  **no association with T-cell infiltration** (r=−0.03, p=0.79) — whereas our proxy gave
  β=−4.10, p<1e-4. The entire run6 structure hinged on the crude CIN proxy.

## E2 — Robustness across 7 infiltration-deconvolution methods
Driver-negative LUAD; LAG3~aneuploidy adjusting for each method's CD8 + purity + TMB.

- **LAG3: 0/7 methods significant.** All 7 coefficients *positive* (0.039–0.173) but all p>0.10.
- **TIGIT: 0/7 significant.** All 7 *negative* (−0.009 to −0.137).

So a **faint, non-significant directional echo persists** (LAG3 positive, TIGIT negative across
every method) — the LAG3>TIGIT pattern isn't random — but it does **not** reach significance.

## Verdict

**The strong per-cell LAG3 claim is not supported by rigorous metrics.** With the citable,
purity-aware Aneuploidy Score and independent infiltration estimates, neither the
infiltration-reduction nor the per-cell LAG3/IFN-γ up-regulation replicates. The run6 result was
substantially an **artifact of the crude CIN proxy** (`cna_mean_abs_log2` / fraction-genome-altered)
and/or the self-computed CD8 score — not a robust biological signal.

What survives is weak: a **consistent-but-non-significant** tendency for LAG3 (not TIGIT) to trend
positive with aneuploidy per T cell (E2, 7/7 same direction). Not disproven, but far below the
confidence run6 implied.

## Why the discordance — and the decisive next test
run6's CIN metric (amplitude / fraction-genome-altered) behaves very differently from the arm-level
**Aneuploidy Score**. Two explanations, separable with one cheap experiment:
- **(a) Purity/proxy artifact:** the FGA proxy is compressed by low tumor purity, which co-varies
  with stromal/immune content — manufacturing a spurious CIN↔immune link.
- **(b) Genuine metric-type biology:** amplitude/focal instability ≠ arm-level aneuploidy; the
  effect could be real for one axis and not the other.

**Proposed E1b (run next):** re-fit the run6 models using the *original* proxy `cna_mean_abs_log2`
as the CIN metric, driver-negative, **with vs without purity adjustment**, alongside Aneuploidy
Score in the same table. If purity adjustment kills the proxy-based effect → (a) artifact. If it
survives with purity but AS still shows nothing → (b) metric-type distinction worth pursuing.

## Bearing on the plan
- Tiers 2–4 (single-cell, replication, clinical) are **lower priority** until E1b clarifies whether
  there is any real effect to chase. The single-cell test (E4–E6) remains the true arbiter of a
  *per-cell* claim, but should be gated on E1b rather than run speculatively.
- This is the intended function of validation: a proxy-driven finding was caught before it drove
  expensive downstream work.
