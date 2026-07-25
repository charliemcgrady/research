# Deliverable #4 — Ranking of unexplained tumors (highest value for future sequencing)

The taxonomy's honest output is not "every tumor is solved." It is a **ranked list of the tumors we
cannot yet mechanistically explain**, ordered so that the ones most worth spending a WGS / single-cell /
spatial run on come first. Full table: `unknown_ranking.csv` (139 rows, all oncogene-negative tumors).

**Ranking metric.** `unexplained_score = 1 − max(C1..C5)`, i.e. one minus the strongest *cell-intrinsic*
explanation. A tumor with a confident hidden-driver / suppressor / protein-state / ligand call scores
low (well-explained); a tumor whose only mass sits on microenvironment, plasticity, or unknown scores
high. Ties broken by ascending best_p (lower confidence first).

## Tiering

| Tier | Definition (max cell-intrinsic evidence) | n (TCGA / CPTAC) | Meaning |
|---|---|---|---|
| **T1 — no intrinsic signal** | < 0.10 | **13** (10 / 3) | No tumor-intrinsic driver evidence at all. Top priority. |
| **T2 — trace only** | 0.10–0.20 | 17 (—) | Negligible intrinsic signal; essentially unexplained. |
| **T3 — weak** | 0.20–0.35 | 36 | Partial hint, not enough to commit to a class. |
| **T4 — explained** | ≥ 0.35 | 73 | A defensible cell-intrinsic mechanism exists. |

**66 of 139 (47%)** oncogene-negative tumors fall below the 0.35 "explained" line — the true size of the
open problem in this subtype. The 13 Tier-1 tumors are the sharpest targets.

## Tier-1 tumors (no cell-intrinsic evidence) — sequence these first

From `unknown_ranking.csv` (top of the sorted list; all have `intrinsic_explanation ≈ 0`):

| Tumor | Cohort | Best (context) call | Dominant lineage | Why it's unexplained |
|---|---|---|---|---|
| TCGA-97-A4LX | TCGA | C7 plasticity (0.54) | EMT/mesenchymal | zero intrinsic lesion; pure high-entropy state |
| C3N-01072 | CPTAC | C6 microenv (0.55) | EMT/mesenchymal | silent genome/transcriptome/**proteome**; context only |
| TCGA-55-8301 | TCGA | C7 plasticity (0.57) | Proliferative | zero intrinsic lesion |
| TCGA-86-8585 | TCGA | C7 plasticity (0.62) | Proliferative | zero intrinsic lesion |
| C3N-00737 | CPTAC | C7 plasticity (0.62) | Proliferative | **even with proteome+phospho, no driver** |
| TCGA-86-6851 | TCGA | C7 plasticity (0.64) | EMT/mesenchymal | zero intrinsic lesion |
| C3N-02003 | CPTAC | C6 microenv (0.75) | Proliferative | proteogenomically silent tumor cell |

(Full 13 in the CSV.) The three **CPTAC** Tier-1 tumors are the most striking: they are unexplained
*despite* having proteome and phosphoproteome — the most genuinely mysterious tumors in the study,
because the usual "we just didn't measure protein" escape hatch is closed for them.

## What kind of tumor is unexplained?

Among the 29 tumors with intrinsic evidence < 0.20, dominant lineage is overwhelmingly:

- **Proliferative** (10) and **EMT/mesenchymal** (9) — 66% of the unexplained tail.
- A minority are AT2/TRU (4), AT1 (2), Club (2).

Best-call is **C7 high-plasticity (19)** or **C6 microenvironment (10)**. This is a coherent picture: the
unexplained tumors are disproportionately **dedifferentiated / high-turnover / mesenchymal** states —
exactly the biology where a single genomic driver is least expected and where transcriptional plasticity
or microenvironmental dependency is the leading hypothesis.

## Recommended assay per unexplained profile (feeds Deliverable #6)

| Unexplained profile | Most informative next assay |
|---|---|
| Tier-1 CPTAC, proteogenomically silent (C3N-00737, -01072, -02003) | **WGS + structural variants** (fusion/enhancer-hijack/complex SV) — the only untested genomic layer |
| Proliferative / EMT, high-entropy (C7) | **Single-cell / spatial** to confirm true plasticity vs admixture; lineage-resolved driver search |
| Microenvironment-dominant (C6), low purity | **Purity-matched or single-cell** re-analysis to recover a diluted intrinsic driver |
| AT2/TRU, quiet | **Deep WGS** for a low-frequency canonical driver missed at panel resolution |

These recommendations are prioritized (Tier-1 first) so a limited sequencing budget goes to the tumors
that would most reduce taxonomy uncertainty.
