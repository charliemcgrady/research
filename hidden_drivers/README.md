# A mechanistic taxonomy of comprehensively oncogene-negative lung adenocarcinoma

**Mission.** Assign *every* comprehensively oncogene-negative (RTK/RAS/RAF wild-type) LUAD tumor to the
most plausible mechanistic class, using convergent multi-source evidence — not one significant p-value.
"Unknown" and "high-plasticity" are legitimate answers. Two independent cohorts; single-cell used only
for lineage interpretation, never as a validation cohort.

## The six deliverables

| # | Deliverable | File(s) |
|---|---|---|
| 1 | **Mechanistic taxonomy** (8 classes, per-cohort prevalence) | [`taxonomy/TAXONOMY.md`](taxonomy/TAXONOMY.md), `taxonomy/mechanism_prevalence.csv` |
| 2 | **Per-tumor evidence graphs** (139 tumors, competing classes + evidence items) | [`taxonomy/PER_TUMOR_EVIDENCE.md`](taxonomy/PER_TUMOR_EVIDENCE.md), `taxonomy/per_tumor_evidence.json`, `taxonomy/tumor_evidence_{TCGA,CPTAC}.csv`, `taxonomy/evidence_items_{TCGA,CPTAC}.csv` |
| 3 | **Per-mechanism evidence graph** (class-level prevalence/reproducibility/support) | [`taxonomy/MECHANISM_EVIDENCE.md`](taxonomy/MECHANISM_EVIDENCE.md) |
| 4 | **Ranking of unexplained tumors** (highest value for future sequencing) | [`taxonomy/UNKNOWN_RANKING.md`](taxonomy/UNKNOWN_RANKING.md), `taxonomy/unknown_ranking.csv` |
| 5 | **Failure report** (every rejected/untestable hypothesis, incl. the CIN→LAG3 arc) | [`FAILURE_REPORT.md`](FAILURE_REPORT.md) |
| 6 | **Knowledge gaps** (the experiment that most reduces uncertainty, per class) | [`KNOWLEDGE_GAPS.md`](KNOWLEDGE_GAPS.md) |

Supporting: `DATASET_ACCESS_REPORT.md`, `dataset_manifest.csv`, `driver_audit_preliminary.csv`,
`FEASIBILITY_GATE.md`, `lineage/lineage_state_table_TCGA.csv`.

## Cohorts & why two

| Cohort | Oncogene-negative n | Adds |
|---|---|---|
| TCGA-LUAD (bulk) | 117 | breadth; RNA/CNA/MAF/ssGSEA/immune/purity |
| **CPTAC-LUAD** (proteogenomic, Gillette 2020) | 22 | **proteome + phosphoproteome** → the only place Class 3 (protein-state activation) and direct ERK/AKT/S6 pathway output are measurable |

## Headline findings

1. **No single mechanism.** Largest cell-intrinsic class is **tumor-suppressor convergence** (~30% TCGA);
   a hidden-canonical tail is ~10–17%.
2. **Protein-state activation is real and DNA/RNA-invisible.** 3/22 CPTAC tumors are driven by an RTK
   (DDR2, EGFR, MET) activated at protein/phospho level with a silent genome and transcriptome — including
   **EGFR with below-average RNA**. Undetectable without proteogenomics.
3. **~47% of oncogene-negative tumors lack any cell-intrinsic explanation** and are best described by
   microenvironment or transcriptional plasticity — ranked for follow-up in Deliverable #4.
4. **Every measurable class reproduces across both cohorts**; positive controls pass (cleaner on CPTAC
   phospho than TCGA ssGSEA proxy).

## Reproduce

```bash
# TCGA engine (bulk): writes tumor_evidence_TCGA.csv + evidence_items_TCGA.csv, prints positive control
uvx --with pandas --with numpy --with openpyxl python3 hidden_drivers/build_tcga_taxonomy.py
# CPTAC engine (adds C3 protein-state): writes tumor_evidence_CPTAC.csv + evidence_items_CPTAC.csv
uvx --with pandas --with numpy --with scipy python3 hidden_drivers/build_cptac_taxonomy.py
# Cross-cohort deliverable products: mechanism_prevalence.csv, per_tumor_evidence.json, unknown_ranking.csv
uvx --with pandas --with numpy python3 hidden_drivers/build_deliverables.py
```

CPTAC raw inputs are open (LinkedOmics `data_download/CPTAC-LUAD/` + cBioPortal `luad_cptac_2020`) and are
gitignored (re-downloadable); see `build_cptac_taxonomy.py` header. TCGA raw inputs are the
patient-provided `asta_autods_data/` set (also gitignored).

## Method notes (carried from the whole program)

- Independent unit is the **tumor**, not the cell or gene.
- Pathway *output* uses **direct phospho** in CPTAC (ERK T185/Y187, AKT, RPS6) to avoid the RNA-derived
  co-expression artifact that sank the earlier per-cell LAG3 claim (see `FAILURE_REPORT.md` §A5/§B3).
- Lineage signatures exclude all candidate driver genes (no circularity).
- Posteriors rank plausibility over a fixed 8-class model; they are **not** null-hypothesis tests — by
  design. Every call ships with its evidence graph so it can be falsified item-by-item.
