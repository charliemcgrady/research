# Deliverable #5 — Failure report: every hypothesis we rejected or could not test

A taxonomy is only trustworthy if its dead ends are documented. This report lists every mechanistic
hypothesis and experiment that was **refuted**, **collapsed under scrutiny**, or **could not be tested on
open data**, across the full arc of this research program (the earlier CIN→immune-escape line and the
current taxonomy build). Rejected hypotheses are as much a deliverable as the taxonomy itself — they
tell the next investigator where *not* to spend effort.

## A. The CIN → immune-escape → TIGIT/LAG3 arc (fully rejected)

The program began from an Asta-Theorizer hypothesis: in oncogene-negative LUAD, chromosomal instability
(CIN) is the oncogenic engine and drives immune escape via the second-generation checkpoints TIGIT and
LAG3. Every load-bearing claim failed under stress-testing (AutoDiscovery → DataVoyager → single-cell):

| # | Hypothesis | Test | Verdict |
|---|---|---|---|
| A1 | CIN up-regulates TIGIT/LAG3 | AutoDS association (bulk) | **Refuted** — CIN-high tumors are checkpoint-**low**; direction reversed. |
| A2 | The effect is driver-negative-**specific** | AutoDS G2; DataVoyager E7 interaction | **Refuted** — same in oncogene-positive; no interaction. |
| A3 | Selective for TIGIT/LAG3 vs PD-1/PD-L1 | AutoDS G3 | **Refuted** — no checkpoint selectivity. |
| A4 | Mediated by cGAS-STING/IFN | AutoDS G4 | **Not supported** — CIN suppresses STING; no formal mediation. |
| A5 | Per-*T-cell* LAG3 ↑ (infiltration-adjusted) | AutoDS run6 (β=2.81, p=0.001) → DataVoyager E1/E1b | **Collapsed** — the "rescue" was a **tumor-purity artifact**; the proxy-CIN metric correlated r≈0.54 with purity. Significant only when adjusted by an expression-derived T-cell score (induced co-expression). |
| A6 | LAG3~CIN replicates per-cell in single-cell | custom inferCNV → CD8 pseudobulk on LuCA (43 primary tumors) | **Null** — β = 0.026, p = 0.56. Clean negative. |

**Net:** CIN does not drive TIGIT/LAG3-mediated immune escape in oncogene-negative LUAD, and this subset
is not a distinct CIN–immune compartment. The one durable correlate of immune "coldness" was **STK11
loss** — a known axis, not a novel finding. This entire line is a documented negative and directly
motivated the pivot to the mechanistic-taxonomy mission.

## B. Methodological failures that changed conclusions (transferable)

- **B1 — CIN-metric choice flips the answer.** An amplitude/FGA-type proxy and the arm-level Taylor
  Aneuploidy Score correlate only r≈0.6 and gave *opposite* significance. Bulk CIN–immunology built on
  FGA proxies is untrustworthy without a purity-aware aneuploidy score.
- **B2 — Tumor purity is a first-order confounder** in every bulk CIN↔immune analysis (inflates proxy-CIN,
  dilutes the RNA immune signal). Must be modeled explicitly (ABSOLUTE purity).
- **B3 — "Per-cell" claims from bulk are fragile.** Adjusting a checkpoint by an RNA-derived infiltration
  score induces residual co-expression that mimics a per-cell effect; deconvolution disagreed with
  regression adjustment; single-cell settled it. This lesson is baked into the current engine — pathway
  *output* in CPTAC uses direct phospho, not an RNA-derived proxy, precisely to avoid B3.

## C. Structural / regulatory hidden-driver classes we could NOT test on open data

The current taxonomy's Class C1 (hidden canonical) and the originally-planned structural experiments
(enhancer hijack, non-coding activation, complex SV, fusion) require WGS + structural variants + matched
RNA at cohort scale. **This does not exist as open data for LUAD.**

| # | Intended experiment | Why blocked |
|---|---|---|
| C1 | Enhancer-hijack / non-coding activation of RTKs | No open LUAD WGS+SV+RNA cohort (all dbGaP/EGA controlled). |
| C2 | Complex structural variants as drivers | Same; CPTAC study exposes **0 structural variants** via cBioPortal. |
| C3 | Fusion detection (ALK/ROS1/RET/NTRK/NRG1) | No open SV/fusion layer; we used an RNA-outlier *proxy* only (conservative, incomplete). |
| C4 | Pre-emption check | Chen et al. 2021 (*Cell Reports*) already performed WGS of RTK/RAS/RAF-negative LUAD — the definitive structural analysis is done and controlled-access. |

**Consequence for the taxonomy:** a tumor whose true driver is a fusion or enhancer hijack will be
mis-assigned (most likely to C3 protein-state, C1, or C7). This is the single largest blind spot and the
top knowledge gap (Deliverable #6, §G1).

## D. Single-cell interrogation of the driver-negative subset (infeasible)

- **D1.** LuCA (Salcher 2022, the largest LUAD atlas) annotates driver status for only ~75/161 patients
  and contains **exactly one** eligible oncogene-negative primary tumor. A single-cell test of any
  *driver-negative-specific* mechanism is impossible with current public data. (Kim GSE131907 full matrix
  was also unretrievable — GEO capped the download at 2.9 GB.)
- **Consequence:** all single-cell use in this project is confined to lineage-signature interpretation
  and the all-LUAD LAG3 null (A6); single cell is **not** used as a validation cohort for the taxonomy,
  per mission constraint.

## E. Taxonomy-engine failure modes we hit and corrected

- **E1 — Entropy over-firing (first pass).** Bulk lineage entropy is high for most tumors, so C7
  initially swept the cohort. Fixed by making C7 a **residual** (`×(1−max intrinsic lesion)`) gated at
  top-quartile entropy. Positive control then passed. Residual risk documented in TAXONOMY.md caveats.
- **E2 — Unknown (C8) under-population.** Because additive evidence often sums >1, little mass reaches C8.
  We therefore report "unexplained" via the **cell-intrinsic-explanation ranking** (Deliverable #4), not
  via C8 count — otherwise the open-problem size (47%) would be hidden.
- **E3 — TCGA positive control is only moderately separated** (C1 mass 0.233 pos vs 0.171 neg; C7 nearly
  ties C1 among positives) because pathway output is an ssGSEA proxy. CPTAC's phospho readout separates
  cleanly. We did **not** paper over this — it is reported as a metric-quality limitation, consistent
  with B1–B3.

## F. Things that could still be wrong (owned uncertainties)

- **F1.** C3 protein-state activation rests on n=3 of 22 in a single cohort. Plausible, not established.
- **F2.** Ligand class (C5) flags an outlier ligand without confirming an active cognate receptor.
- **F3.** Fusion under-detection (C-section) can masquerade as C1/C3/C7.
- **F4.** Bulk C6/C7 calls can be purity/admixture artifacts (B2) rather than true biology.

Every one of these is carried forward into Deliverable #6 as a concrete, fundable next experiment.
