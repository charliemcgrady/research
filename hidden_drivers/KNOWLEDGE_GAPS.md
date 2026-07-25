# Deliverable #6 — Knowledge gaps: the single experiment that would most reduce uncertainty, per class

For each mechanistic class (and each major failure mode), the one dataset or experiment that would move
it from "plausible" to "established" — or that would resolve the tumors we cannot currently explain.
Ordered by expected information gain for the taxonomy as a whole.

## Ranked by impact

### G1 — Open (or requestable) LUAD **WGS + structural-variant + RNA** cohort  *(highest impact)*
- **Resolves:** C1 (hidden canonical), the entire blocked structural class (failure report §C), fusion
  under-detection (F3), and most Tier-1 unexplained tumors (Deliverable #4).
- **Why it's the top gap:** every genomic mechanism *except* SNV/CNA/expression is currently untestable
  on open data. Enhancer hijack, complex SV, and fusions are the most likely true drivers of the
  hidden-canonical and unexplained tumors, and they are exactly what we cannot see.
- **Concrete ask:** dbGaP/EGA access to the Chen et al. 2021 RTK/RAS/RAF-negative WGS cohort, or TCGA/
  CPTAC controlled-tier WGS + Manta/GRIDSS SV calls; re-run C1 with true SV/fusion evidence.

### G2 — A **second proteogenomic LUAD cohort** (CPTAC-3 / APOLLO / ProTrait)
- **Resolves:** C3 (protein-state activation) reproducibility — currently n=3/22, single cohort.
- **Why:** protein-state activation is the study's highest-value, DNA/RNA-invisible mechanism; it can
  only be confirmed proteogenomically. A second cohort tests whether DDR2/EGFR/MET protein-state
  activation recurs and at what prevalence.
- **Concrete ask:** run the identical C3 outlier logic on an independent phosphoproteomic LUAD set;
  cross-check each C3 tumor's RTK phospho against its downstream ERK/AKT/S6 phospho to confirm signaling.

### G3 — **Single-cell / spatial** LUAD with driver annotation
- **Resolves:** C7 (high-plasticity) — distinguishing true single-cell transcriptional plasticity from
  bulk admixture (failure §E1, F4) — and C4 (lineage-conditioned), which we suspect is a real modifier
  masked at bulk resolution.
- **Why:** the unexplained tail is disproportionately proliferative/EMT (Deliverable #4); only single-cell
  can tell whether those are genuinely plastic tumor cells or mixtures. Also the only way to revisit any
  driver-negative-*specific* mechanism (failure §D).
- **Concrete ask:** a driver-annotated scRNA/spatial cohort with ≥10 oncogene-negative primary tumors —
  which, per the LuCA audit, **does not currently exist** and would have to be generated.

### G4 — **Functional dependency screens** (CRISPR / drug) in oncogene-negative models
- **Resolves:** whether C2 (suppressor convergence) tumors are actually MAPK/PI3K-*dependent*, and whether
  C3 protein-state RTKs are druggable dependencies vs bystander phosphorylation.
- **Why:** the taxonomy asserts mechanism from association; dependency is the missing causal layer.
- **Concrete ask:** DepMap-style RTK/pathway dependency in RTK/RAS/RAF-WT LUAD lines stratified by class;
  test DDR2/MET/EGFR inhibitor sensitivity in protein-state-high, DNA-silent models.

### G5 — **Receptor–ligand pairing** for the ligand class
- **Resolves:** C5 (ligand/autocrine) — confirm each ligand outlier drives an active cognate receptor
  (failure §F2).
- **Concrete ask:** in CPTAC, pair each C5 ligand outlier with its receptor's phosphorylation state;
  where both are high, the autocrine call is upgraded from suggestive to supported.

### G6 — **Purity-matched re-analysis** for the microenvironment class
- **Resolves:** C6 (microenvironment) — separate genuine microenvironment-dependent biology from
  low-purity dilution of an intrinsic driver (failure §B2, F4).
- **Concrete ask:** stratify C6 tumors by ABSOLUTE purity; for low-purity C6 tumors, re-run the intrinsic
  classes on tumor-cell-deconvolved expression to recover any masked driver.

## Gap-to-class matrix

| Class | Primary gap | Status without it |
|---|---|---|
| C1 hidden canonical | G1 (WGS/SV) | tail of hidden drivers under-detected |
| C2 suppressor convergence | G4 (dependency) | mechanism asserted, dependency unproven |
| C3 protein-state | G2 (2nd proteogenomic cohort) | n=3, single-cohort hypothesis |
| C4 lineage-conditioned | G3 (single-cell) | can't be promoted from modifier |
| C5 ligand/autocrine | G5 (receptor pairing) | outlier ≠ active signaling |
| C6 microenvironment | G6 (purity) | may hide diluted intrinsic drivers |
| C7 high-plasticity | G3 (single-cell) | plasticity vs admixture unresolved |
| C8 / unexplained tail | G1 then G3 | 47% of tumors below the "explained" line |

## One-line program recommendation

If exactly one dataset could be added, it is **G1 (open LUAD WGS + structural variants + RNA)**: it is the
only layer that is both (a) currently missing entirely from open data and (b) the most probable location
of the true drivers for the hidden-canonical and unexplained tumors that make up the largest uncertainty
in the taxonomy. Proteogenomic breadth (G2) is the close second, because it is the *only* way to confirm
the one mechanism — protein-state activation — that no DNA/RNA study can ever see.
