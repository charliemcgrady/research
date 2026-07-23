# Driver-Negative LUAD — Literature Gap Analysis (Final)

**Hypothesis under test:** In oncogene/driver-negative lung adenocarcinoma (LUAD),
chromosomal instability (CIN) drives the tumors and mediates immune escape, with
**TIGIT and LAG3** as checkpoints of interest (panel-aware vs PD-1/PD-L1/TIM-3/CTLA-4).

**Method:** 4 Asta Paper Finder searches (`asta literature find`), **all diligent mode**.
**615 unique papers**, ~74% published 2021+. (q1/q2 were first run in fast mode when
diligent hit transient gateway upstream timeouts; both were later completed in diligent
mode on retry — numbers below reflect the deeper diligent corpus.)

| ID | Focus | Mode | Papers |
|----|-------|------|-------:|
| q1 | Driver-negative LUAD landscape & immune evasion | diligent | 104 |
| q2 | CIN → immune evasion mechanism (cGAS-STING/IFN/antigen) | diligent | 281 |
| q3 | TIGIT/LAG3 in NSCLC vs full checkpoint panel | diligent | 130 |
| q4 | CIN/aneuploidy ↔ checkpoint expression in LUAD | diligent | 129 |

## What is ESTABLISHED (not gaps)

1. **CIN → immune escape is a mature, mechanistic field — but almost entirely
   pan-cancer / general, not driver-stratified LUAD.** Anchor papers:
   - *Non-cell-autonomous cancer progression from chromosomal instability* (Nature 2023):
     CIN → chronic **cGAS–STING** activation → metastasis and immune modulation.
   - *Genomic instability as a driver and suppressor of anti-tumor immunity*
     (Front Immunol 2024): the **double-edged sword** — acute CIN is immunogenic,
     chronic CIN drives evasion.
   - *Suppression of tumor antigen presentation during aneuploid tumor evolution*
     (Oncoimmunology 2019): CIN-high tumors are less immunogenic / less infiltrated.
2. **TIGIT and LAG3 as NSCLC checkpoints/targets is well-covered** (TIGIT n=62, LAG3 n=70),
   including reviews, prognostic meta-analyses, and trial-oriented pieces — but as
   PD-1-adjacent targets, **not tied to genomic instability**.
3. **Driver/oncogene-negative LUAD is a recognized but thin niche** (15 of 615 papers).
   The biology anchors:
   - *Combinatorial inactivation of tumor suppressors efficiently initiates lung
     adenocarcinoma...* (Cancer Research 2021/2022): "a large fraction of lung
     adenocarcinomas lacks mutations in known oncogenes... these oncogene-negative tumors
     remain enigmatic" — and they arise from **tumor-suppressor loss** (TP53/RB1/etc.),
     itself a canonical CIN driver.
   - *Lung adenocarcinomas without driver genes converge to common adaptive strategies* (2024).
   - *Fully accessible fitness landscape of oncogene-negative lung adenocarcinoma* (2023).
   These characterize genesis/evolution but say nothing about the immune axis or checkpoints.

## The GAP (intersection tests across all 615 papers)

| Junction | # papers | Note |
|----------|---------:|------|
| cGAS-STING / micronuclei | 267 | very large, general |
| PD-1 / PD-L1 | 170 | very large |
| Chromosomal instability / aneuploidy | 108 | large, general |
| LAG3 | 70 | large, NSCLC |
| TIGIT | 63 | large, NSCLC |
| Oncogene/driver-negative lung | 15 | thin |
| **CIN AND (TIGIT or LAG3)** | **2** | both generic pan-cancer bioinformatics |
| **CIN AND driver-negative** | **1** | closest prior work (see below) — no immune/checkpoint axis |
| **Driver-negative AND (TIGIT/LAG3)** | **0** | — |

**Robust to search depth.** Doubling the corpus (317 → 615) exploded the *general*
mechanism literature — cGAS-STING went 64 → 267 — while the hypothesis junctions barely
moved (CIN∩TIGIT/LAG3 stayed at 2; driver-neg∩TIGIT/LAG3 stayed at 0). The white space is
not a sampling artifact.

**Closest prior work** (the single CIN∩driver-negative hit): *Proteogenomic characterization
identifies clinical subgroups in EGFR and ALK wild-type never-smoker lung adenocarcinoma*
(2024) — profiles wild-type LUAD including copy-number/CIN features, but does **not** address
immune escape, TIGIT, or LAG3. It is the nearest neighbor to your hypothesis and the natural
citation to differentiate against.

**Interpretation.** Three mature literatures — CIN→immune-escape, TIGIT/LAG3, and
(nascent) driver-negative LUAD — run in near-total isolation. The specific, mechanistically
motivated chain you propose sits in a **genuine white space**:

> Tumor-suppressor-loss-driven (driver-negative) LUAD is *expected* to be CIN-high;
> CIN drives immune escape via chronic cGAS–STING; yet **no study tests whether this
> selectively up-regulates TIGIT/LAG3 in the driver-negative subset.** Zero papers close
> that loop.

## Recommended AutoDiscovery scope (ranked)

**Primary (highest novelty, data purpose-built):**
> *In driver-negative LUAD, does chromosomal-instability burden correlate with — and
> potentially drive — TIGIT/LAG3 up-regulation, more strongly than in driver-positive
> LUAD and more selectively than for PD-1/PD-L1?*

Testable end-to-end in the TCGA data:
- **CIN metric** ← `cleaned_cna (1).csv` (fraction genome altered / arm-level aneuploidy).
- **Driver status** ← `ON classification TCGA.xlsx` (`Previous_RPA_Positive`: 0 = driver-neg, n=118; 1 = driver-pos, n=383).
- **Checkpoints** ← `cleaned_gene_expression.csv.gz` (TIGIT ENSG00000181847, LAG3 ENSG00000089692, + PD-1/PD-L1/TIM-3/CTLA-4).
- **Mechanism readout** ← ssGSEA `HALLMARK_INTERFERON_GAMMA/ALPHA_RESPONSE`,
  `INFLAMMATORY_RESPONSE` (proxy for the cGAS-STING/IFN axis).
- **Confounder** ← `cleaned_tmb.csv` (TMB co-varies with both CIN and immunity — must adjust).

**Secondary gaps worth an AutoDS arm each:**
1. Is the CIN→checkpoint effect *specific* to driver-negative tumors? (interaction test)
2. Does CIN elevate TIGIT/LAG3 *selectively* vs PD-1/PD-L1? (differential axis)
3. Is the CIN→checkpoint link *mediated* by the IFN/cGAS-STING hallmark signature?
   (mediation, using the ssGSEA scores)

## Caveats
- All four searches are now **diligent mode** (615 papers). q1/q2 initially fell back to
  fast mode when diligent hit transient gateway upstream timeouts, but both completed on
  retry; the deeper corpus only strengthened the gap (intersections stayed at 0–2 while the
  general literature roughly doubled).
- "Driver-negative" here = RTK/RAS/RAF-pathway wild-type per the TCGA RPA classification;
  other definitions (e.g. including STK11/KEAP1 status) could shift the subset.
- Intersection counts are keyword-based on title+abstract, so they undercount papers that
  discuss a junction only in full text. The *direction* (near-empty intersections against
  large single-concept pools) is robust, but exact counts are lower bounds.
