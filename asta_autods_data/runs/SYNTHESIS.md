# AutoDiscovery Synthesis — Does CIN drive immune escape via TIGIT/LAG3 in driver-negative LUAD?

**5 runs · 100 experiments · all SUCCEEDED · TCGA-LUAD (n≈115 driver-negative with complete omics).**

## Headline

**The theorizer's hypothesis is refuted — and in the opposite direction.** Across all
five runs, chromosomal instability (CIN) in driver-negative LUAD is robustly associated
with an immune-**cold** microenvironment: **lower** TIGIT/LAG3 (and other checkpoints) and
strongly **suppressed** interferon-γ response — not the predicted CIN→cGAS-STING→checkpoint-high
exhaustion. This is consistent with the aneuploidy→immune-exclusion arm of the literature
(the "double-edged sword"), not the checkpoint-mediated-escape arm.

> ⚠️ **Critical caveat (bulk RNA):** TIGIT/LAG3 are expressed mainly on infiltrating T/NK
> cells. Lower *bulk* checkpoint expression in CIN-high tumors most likely reflects **fewer
> infiltrating lymphocytes**, not lower exhaustion *per T cell*. The single most important
> follow-up is to normalize checkpoints to a CD8/T-cell-infiltration signature (or deconvolve)
> and re-test — that could still rescue a per-cell version of the original hypothesis.

## Per-gap verdicts

| Gap | Run | Verdict | Key statistic |
|-----|-----|---------|---------------|
| **G1 Association** | run1 | ❌ **Refuted (opposite sign)** — CIN is *negatively* associated with TIGIT/LAG3 | TIGIT vs CN gain β=−3.23 (p=0.009), loss β=−2.81 (p=0.005); no inverted-U (quadratic p=0.09/0.11) |
| **G2 Specificity** | run2 | ❌ **Not specific to driver-negative** — CIN suppresses checkpoints regardless of driver status | CIN main effect on TIGIT β=−3.33 (p<0.001); CIN×driver interaction p=0.74–0.83 (ns) |
| **G3 Selectivity** | run3 | ❌ **No selectivity** — CIN hits TIGIT/LAG3 no differently than PD-1/PD-L1 | Wald test TIGIT p=0.93, LAG3 p=0.85, PD-1 p=0.96 |
| **G4 Mediation** | run4 | ❌ **Not mediated by cGAS-STING/IFN** — direct association, mechanism unresolved | CIN→STING1 β=−2.5 (p=3e-10) & STING1→TIGIT (p=0.03), but indirect effect (ACME) ns; CXCL9 mediation ns |
| **G5 CN events** | run5 | ✅ **Global CIN + STK11 loss drive the cold state; not 9p21** | CIN→IFN-γ β=−2008 (p<0.001); only ~17% via focal dels; STK11 loss significant; CDKN2A/B ns |

## The cross-run picture

1. **Direction:** CIN → immune-cold, robustly and repeatedly. The core association (G1),
   the driver-status contrast (G2), and the mechanism runs (G4) all recover a **strong,
   significant negative** CIN→checkpoint / CIN→IFN-γ effect.
2. **Not selective, not driver-specific:** the effect is a general immune-coldness phenomenon,
   not a TIGIT/LAG3-specific or oncogene-negative-specific one. (This actually *answers* Gaps
   2–3 — just not in the hypothesized direction.)
3. **Mechanism:** cGAS-STING is *suppressed* by CIN (consistent with chronic-CIN STING
   silencing) and correlates with checkpoints, but does **not** formally mediate the
   CIN→checkpoint link in this bulk cohort. Left unresolved.
4. **Genomic drivers (G5):** the cold state is driven by *global* instability, with **STK11
   loss** the standout focal contributor — echoing the well-known STK11/KEAP1 "cold tumor"
   biology. 9p21/CDKN2A loss was **not** a significant driver (contra a common assumption).

## What this means for the driver-negative literature gap

The literature gap we identified is real, but the data fills it with a *negative/reversal*
result: in driver-negative LUAD, CIN does **not** appear to drive TIGIT/LAG3-mediated immune
escape; it tracks with immune exclusion. That is still a novel, publishable direction — and it
sharpens the next question.

## Recommended next AutoDS runs
1. **Infiltration-normalized re-test (highest priority):** regress TIGIT/LAG3 on CIN *after*
   adjusting for a CD8/T-cell or immune-infiltration score (e.g. a cytolytic/`ss_ALLOGRAFT_REJECTION`
   proxy, or add a deconvolution estimate). Tests exhaustion *per T cell*.
2. **STK11-loss cold program:** characterize the STK11-driven immunosuppressive program in
   driver-negative tumors specifically.
3. **CIN-high + immune-hot outliers:** identify the driver-negative tumors that are CIN-high
   yet checkpoint-high (residual analysis) — the exceptions may reveal the escape route.

*Run IDs and configs in `RUNS.md`; per-experiment JSON in `runs/exp_<runid>.json`.*
