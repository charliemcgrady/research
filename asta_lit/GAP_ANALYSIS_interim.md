# Driver-Negative LUAD — Literature Gap Analysis (Interim)

**Hypothesis under test:** In oncogene/driver-negative lung adenocarcinoma (LUAD),
chromosomal instability (CIN) drives the tumors and mediates immune escape, with
TIGIT and LAG3 as checkpoints of interest (panel-aware vs PD-1/PD-L1/TIM-3/CTLA-4).

**Status:** 2 of 4 Asta Paper Finder searches complete (q3, q4). The landscape (q1)
and CIN→immune-mechanism (q2) searches were rate-limited by the gateway and are
scheduled for automatic retry. Findings below are therefore *provisional* and biased
toward the checkpoint/CIN-intersection literature.

## Searches run (Asta `literature find`, diligent mode)

| ID | Query focus | Papers | % from 2021+ |
|----|-------------|-------:|-------------:|
| q3 | TIGIT & LAG3 in NSCLC/LUAD vs PD-1/PD-L1/TIM-3/CTLA-4 | 130 | 78% |
| q4 | CIN/aneuploidy ↔ TIGIT/LAG3 checkpoint expression in LUAD | 129 | 88% |
| q1 | Driver-negative LUAD landscape & immune evasion | *pending retry* | — |
| q2 | CIN → immune evasion mechanism (cGAS-STING, IFN, antigen pres.) | *pending retry* | — |

239 unique papers across q3+q4.

## Concept coverage (of 239 unique papers)

| Concept | # papers |
|---------|---------:|
| PD-1 / PD-L1 | 100 |
| LAG3 | 69 |
| TIGIT | 62 |
| TIM-3 / CTLA-4 | 61 |
| TCGA / bioinformatic signature work | 61 |
| cGAS-STING | 44 |
| Interferon / antigen presentation | 35 |
| Chromosomal instability / aneuploidy | 26 |
| T-cell exhaustion | 25 |
| KRAS/EGFR (driver) context | 21 |
| **Oncogene/driver-negative lung cancer** | **2** |

### Intersection tests (the actual hypothesis junctions)
- **CIN AND (TIGIT or LAG3) in the same paper: 2** — both generic/pan-cancer
  bioinformatics (a pan-solid-tumor LAG3 meta-analysis; a 2026 pan-cancer CD70/CD80/TIGIT
  multi-omics characterization). Neither is mechanistic, neither is driver-negative-specific.
- **CIN AND driver-negative together: 0**
- **Driver-negative mentioned at all: 2** — both in the context of immunotherapy
  *resistance* generally, not CIN.

## Provisional gaps (candidates for AutoDiscovery)

1. **CIN → TIGIT/LAG3 axis is essentially unmapped.** TIGIT and LAG3 are heavily
   studied as checkpoints and CIN is studied as an immune-evasion driver, but the two
   literatures run in parallel. Almost no work asks whether CIN burden *predicts* or
   *mechanistically drives* TIGIT/LAG3 up-regulation — and none in LUAD specifically.

2. **Driver-negative stratification is absent from the checkpoint literature.** The
   checkpoint/CIN work is not stratified by driver status; the oncogene-negative subset
   (n=118 in your TCGA data) is not treated as a distinct immuno-genomic entity.

3. **TIGIT/LAG3 vs the rest of the panel under high CIN.** No evidence yet on whether
   CIN selectively elevates TIGIT/LAG3 relative to PD-1/PD-L1 — directly testable in your
   data and a clean differentiator.

4. **cGAS-STING as the mediating mechanism (44 papers) is never closed to a checkpoint
   readout in driver-negative LUAD.** The mechanistic chain CIN → micronuclei →
   cGAS-STING/chronic IFN → checkpoint induction is proposed in the abstract but not
   demonstrated end-to-end in this subset.

> These gaps are provisional pending q1/q2. In particular, the "driver-negative is
> barely studied" signal is partly an artifact of q3/q4 being checkpoint/CIN-focused —
> q1 (the dedicated driver-negative landscape search) is needed before asserting the
> driver-negative biology itself is a gap.

## How this maps onto your TCGA data (`asta_autods_data/`)
- **CIN score**: derive from `cleaned_cna (1).csv` (fraction genome altered / arm-level aneuploidy).
- **Driver status**: `ON classification TCGA.xlsx` → `Previous_RPA_Positive` (0 = driver-negative, n=118).
- **Immune-escape readouts**: TIGIT (ENSG00000181847), LAG3 (ENSG00000089692) + full panel
  in `cleaned_gene_expression.csv.gz`; Hallmark immune pathways in the ssGSEA file.
- **Confounders to control**: TMB (`cleaned_tmb.csv`), since CIN and TMB both affect immunity.
