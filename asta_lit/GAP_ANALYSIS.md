# Driver-Negative LUAD — Literature Gap Analysis (Final)

**Hypothesis under test:** In oncogene/driver-negative lung adenocarcinoma (LUAD),
chromosomal instability (CIN) drives the tumors and mediates immune escape, with
**TIGIT and LAG3** as checkpoints of interest (panel-aware vs PD-1/PD-L1/TIM-3/CTLA-4).

**Method:** 4 Asta Paper Finder searches (`asta literature find`). q3/q4 diligent mode;
q1/q2 fast mode (diligent repeatedly hit the gateway's upstream timeout on these
broader queries). **317 unique papers**, ~80% published 2021+.

| ID | Focus | Mode | Papers |
|----|-------|------|-------:|
| q1 | Driver-negative LUAD landscape & immune evasion | fast | 43 |
| q2 | CIN → immune evasion mechanism (cGAS-STING/IFN/antigen) | fast | 43 |
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
3. **Driver/oncogene-negative LUAD is a recognized but thin niche.** The key biology
   anchor — *Combinatorial inactivation of tumor suppressors efficiently initiates lung
   adenocarcinoma...* (Cancer Research 2022) — states plainly: "a large fraction of lung
   adenocarcinomas lacks mutations in known oncogenes, and the genesis and treatment of
   these oncogene-negative tumors remain enigmatic," and shows these tumors arise from
   **tumor-suppressor loss** (TP53/RB1/etc.) — which is itself a canonical CIN driver.

## The GAP (intersection tests across all 317 papers)

| Junction | # papers | Note |
|----------|---------:|------|
| Chromosomal instability / aneuploidy | 61 | large, general |
| cGAS-STING / micronuclei | 64 | large, general |
| TIGIT | 62 | large, NSCLC |
| LAG3 | 70 | large, NSCLC |
| Oncogene/driver-negative lung | 7 | thin |
| **CIN AND (TIGIT or LAG3)** | **2** | both generic pan-cancer bioinformatics |
| **CIN AND driver-negative** | **0** | — |
| **Driver-negative AND (TIGIT/LAG3)** | **0** | — |

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
- q1/q2 are **fast mode** (43 papers each) — lighter than diligent. Diligent kept hitting
  the gateway's upstream timeout on these broad queries; a later diligent re-run would
  deepen landscape coverage, but the intersection **zeros** are robust and unlikely to move.
- "Driver-negative" here = RTK/RAS/RAF-pathway wild-type per the TCGA RPA classification;
  other definitions (e.g. including STK11/KEAP1 status) could shift the subset.
