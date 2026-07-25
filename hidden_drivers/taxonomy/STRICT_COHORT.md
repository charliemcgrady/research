# Strict-cohort sensitivity analysis — TCGA taxonomy excluding 28 WGS-reclassified tumors

The TCGA "oncogene-negative" set has a soft boundary: 28 of 118 RPA-negative tumors were assigned a driver
alteration on WGS re-analysis (`Current_RPA_Alteration`). This analysis (`build_tcga_taxonomy_strict.py`)
excludes those 28 to give a WGS-strict oncogene-negative cohort (**n = 89**) and moves them into the
oncogene-positive control. Per-tumor posteriors are unchanged (computed against cohort-wide references);
what changes is the population.

Throughout, findings are tagged **[OBSERVATION]** (directly in the data), **[INFERENCE]** (mechanistic
reading), or **[HYPOTHESIS]** (requires future validation).

> **Framing note.** The 28 tumors were **not** held out from threshold selection, feature engineering,
> weighting, or taxonomy construction — the engine's thresholds were set on the full cohort that included
> them. They are therefore a **retrospective calibration** set, **not** a held-out validation set. All
> language below reflects that.

---

## 1. Primary result — robustness  [OBSERVATION]

Excluding the 28 subsequently-reclassified tumors does not change the inferred mechanistic landscape:

| Class | Full (n=117) | Strict (n=89) | Δ pts |
|---|---|---|---|
| Tumor-suppressor convergence (C2) | 29.9% | **31.5%** | +1.6 |
| Residual / state-entropy (C7) | 27.4% | 24.7% | −2.7 |
| Microenvironment (C6) | 19.7% | 20.2% | +0.5 |
| Hidden canonical (C1) | 17.1% | 16.9% | −0.2 |
| Ligand / autocrine (C5) | 4.3% | 4.5% | +0.2 |
| Lineage-conditioned (C4) | 0.9% | 1.1% | — |
| Unknown (C8) | 0.9% | 1.1% | — |

- Class ordering is unchanged.
- Every class shifts ≤ 3 percentage points.
- Posterior mass is essentially unchanged (C2 0.259→0.262; C7 0.223→0.217; C1 0.171→0.173).
- **C2 remains the dominant mechanistic class.**

**This supports that the taxonomy is not an artifact of the 28 subsequently reclassified tumors.** Data:
`strict_vs_full_comparison.csv`, `tumor_evidence_TCGA_strict.csv`.

---

## 2. Retrospective calibration — what the 28 excluded tumors reveal

### 2a. Drivers visible to the engine's input modalities are assigned to concordant classes  [OBSERVATION]

The engine's inputs are RNA, CNA, and (proxy) fusion signal. For the 8 excluded drivers whose consequences
should be visible in that feature space, we predefined a "correct mechanistic match" **before** inspecting
predictions:

| True WGS mechanism | Predefined concordant class(es) |
|---|---|
| RTK amplification (EGFR) | C1 |
| MAPK-pathway amplification (MAPK1, ARAF) | C1 or C2 |
| Ligand fusion (NRG1) | C5 |
| RAS-suppressor deletion (NF1, RASA1) | C2 |

**Capture rate = 7 / 8 = 87.5% (95% Clopper-Pearson CI 47.3–99.7%).** [OBSERVATION]

| Tumor | Alteration | Predicted | Concordant? |
|---|---|---|---|
| TCGA-50-5939 | EGFR_amp | C1 | ✅ |
| TCGA-78-7155 | MAPK1_amp | C1 | ✅ |
| TCGA-86-8358 | MAPK1_amp | C1 | ✅ |
| TCGA-86-7955 | MAPK1_amp | C2 | ✅ |
| TCGA-69-8255 | ARAF_amp | C2 | ✅ |
| TCGA-55-8614 | RASA1_del | C2 | ✅ |
| **TCGA-86-7954** | **NRG1_Fusion** | **C5** | ✅ |
| TCGA-49-4507 | NF1_del | C6 | ❌ (expected C2) |

The **NRG1 fusion** case is the strongest single demonstration: the engine assigned this tumor to the
ligand-driven class **from elevated NRG1 signalling alone**, before any fusion annotation was used.
[OBSERVATION] This indicates the taxonomy captures a biologically meaningful *consequence* (ligand-pathway
activation) rather than memorising a genomic label. [INFERENCE]

The wide confidence interval (n = 8) means the point estimate of 87.5% is **not precise**; it establishes
concordance in principle, not a validated capture rate. [OBSERVATION]

### 2b. Drivers invisible to the feature space fall into nonspecific classes  [OBSERVATION]

20 of the 28 excluded tumors carry **KRAS point mutations**, whose activating consequence is not directly
observable from RNA/CNA/fusion evidence. Within the 28 excluded tumors:

|  | nonspecific (C6/C7) | specific (C1/C2/C5) |
|---|---|---|
| KRAS point-mutant | 14 | 6 |
| non-KRAS (visible) | 1 | 7 |

**Fisher exact (one-sided) p = 0.0087; Haldane odds ratio = 11.2 (95% CI 1.5–81.0).** [OBSERVATION]

**Correct statement of what this shows** [INFERENCE]:
> Tumors whose initiating lesion is invisible to the current feature space (e.g. KRAS point mutations)
> preferentially fall into nonspecific mechanistic classes such as C7 and C6.

This is a **calibration result that defines the limits of the inference engine.** It does **not** establish
that the remaining C7/C6 tumors in the strict cohort contain hidden point mutations — that is an open
hypothesis, not a demonstrated fact. [HYPOTHESIS] The wide OR interval (small n) further cautions against
a strong quantitative claim.

### 2c. Posterior calibration  [OBSERVATION]

Among the 8 visible drivers, concordant assignments had mean best_p 0.59 (range 0.46–0.80) vs 0.52 for the
single discordant assignment. With only **one** discordant case, this is **underpowered** and cannot
establish that correct assignments receive higher posteriors. **Consequently, posterior values in this
taxonomy should be treated as ranking scores, not calibrated probabilities.**

### 2d. Mechanistic confusion matrix  [OBSERVATION]

Rows = true WGS mechanism, columns = predicted class (`strict_confusion_matrix.csv`):

| True mechanism \ Predicted | C1 | C2 | C5 | C6 | C7 |
|---|---|---|---|---|---|
| RTK amplification | **1** | 0 | 0 | 0 | 0 |
| MAPK-pathway amplification | **2** | **2** | 0 | 0 | 0 |
| Ligand fusion | 0 | 0 | **1** | 0 | 0 |
| RAS-suppressor deletion | 0 | **1** | 0 | 1 | 0 |
| **RAS point mutation** | 2 | 4 | 0 | **4** | **10** |

The blind spot (KRAS point mutation → C6/C7) and the concordant diagonal for amplifications/fusion/deletion
are both immediately visible.

---

## 3. Interpretation of C7 (revised)

C7 is best described **not** as "high-plasticity / unexplained biology" but as a **residual class
representing tumors whose mechanism cannot be confidently inferred from currently available public RNA/CNA/
protein evidence.** Such tumors may harbour any of:

- hidden point mutations (as the KRAS calibration shows is possible),
- regulatory / non-coding mechanisms,
- multiple weak cooperating alterations,
- measurement limitations (purity, admixture), or
- genuinely novel biology.

**The present analysis cannot distinguish among these possibilities.** C7 (and C6) are therefore flagged as
**high-priority targets for future genomic characterisation**, not as a positive biological finding.

---

## 4. Manuscript summary

> The strict-cohort sensitivity analysis strengthens confidence in the proposed taxonomy. Excluding tumors
> later reclassified by WGS leaves the inferred mechanistic landscape essentially unchanged (class ordering
> preserved, every class shifting ≤3 percentage points, posterior mass essentially constant, C2 remaining
> dominant), indicating that the taxonomy is robust to stricter definitions of oncogene-negative LUAD. The
> excluded tumors also provide a valuable retrospective calibration dataset. Drivers with consequences
> visible to the engine's input modalities — amplifications, deletions, and fusion-associated signalling —
> are assigned to biologically concordant mechanistic classes (7/8; NRG1 fusion→C5 recognised from ligand
> signalling alone). In contrast, tumors driven by KRAS point mutations preferentially occupy residual
> classes C6/C7 (Fisher p=0.009; OR 11.2, 95% CI 1.5–81), consistent with the engine's inability to observe
> point mutations from RNA/CNA-derived evidence. These results define both the strengths and the
> modality-dependent limitations of the taxonomy while identifying C6 and C7 as high-priority targets for
> future genomic characterisation. Given the small calibration set and single discordant case, posteriors
> should be interpreted as ranking scores rather than calibrated probabilities.

Artifacts: `strict_vs_full_comparison.csv`, `strict_confusion_matrix.csv`, `excluded28_audit.csv`,
`tumor_evidence_TCGA_strict.csv`, `evidence_items_TCGA_strict.csv`.
