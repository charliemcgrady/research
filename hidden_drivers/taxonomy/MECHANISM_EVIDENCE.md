# Deliverable #3 — Per-mechanism evidence graph

For each of the eight classes: how prevalent it is, how reproducible across the two platforms, what
evidence typically supports it, how strong that evidence is, and where it is weakest. This is the
"class-level" complement to the per-tumor graphs — read it as a portrait of each mechanism as a
population phenomenon. Prevalence numbers are from `mechanism_prevalence.csv`.

Legend: **best** = tumors whose top call is this class; **support** = tumors with material mass (>0.1) for
this class even if not the winner; **mass** = mean posterior mass across the cohort.

---

## C1 — Hidden canonical driver
- **Prevalence.** TCGA 20 best (17.1%), CPTAC 2 (9.1%). Mean mass 0.171 / 0.154.
- **Typical evidence.** RTK RNA over-expression outlier; sub-threshold focal amplification of EGFR/ERBB2/
  MET/FGFR1 (below the oncogene-positive cut). Occasionally a lone dominant amp item (C3N-02000: w=0.97).
- **Reproducibility.** ✅ both cohorts; positive-control validated (mass higher in oncogene-*positive*
  tumors, 0.23–0.25, than negative).
- **Strength / limitation.** Strong when an amplification is present; weaker when resting on an RNA
  outlier alone (could be lineage-driven expression). This is the class most likely to be *confirmed* by
  deeper genomics — these tumors probably do have a canonical driver we can't see at panel resolution.

## C2 — Tumor-suppressor convergence  ← largest cell-intrinsic class
- **Prevalence.** TCGA 35 best (**29.9%**, mass 0.259) — the single largest intrinsic class. CPTAC 2 (9.1%).
- **Typical evidence.** ≥1 of NF1/RASA1/PTEN/STK11/KEAP1 mutated or deep-lost, combined with elevated
  MAPK/PI3K output (ssGSEA in TCGA; **direct ERK/AKT/S6 phospho** in CPTAC). Textbook case C3N-01823 hits
  NF1 + STK11 + KEAP1 simultaneously.
- **Reproducibility.** ✅ both. The CPTAC phospho readout makes the "loss → pathway output" logic
  mechanistically tighter than the TCGA proxy.
- **Strength / limitation.** The most defensible "these are RAS-pathway tumors after all" story;
  converges with Chen 2021 WGS reclassification. Limitation: suppressor loss is common in LUAD generally,
  so convergence is suggestive of mechanism, not proof of dependency — functional validation needed.

## C3 — Protein-state activation  ← the proteogenomic-only class
- **Prevalence.** TCGA **not measurable** (no proteome). CPTAC 3 best (**13.6%**, mass 0.070); 4 tumors
  with support >0.1.
- **Typical evidence.** An RTK (DDR2, EGFR, MET) or downstream node activated at the protein/phospho
  level (z>2) while DNA (no amp/mut) and RNA (z<1.5) are silent. The three winners: DDR2 z=3.6 (RNA 0.1),
  EGFR z=3.2 (RNA −0.3), MET z=3.0.
- **Reproducibility.** By construction CPTAC-exclusive; correctly enriched in oncogene-negative vs
  -positive tumors. Cross-cohort reproduction is impossible without a second proteogenomic LUAD cohort.
- **Strength / limitation.** Scientifically the highest-value class: these tumors are **invisible** to
  every DNA/RNA-only study, so this mechanism can only be discovered proteogenomically. Limitation:
  n=3 on n=22; a hypothesis to be reproduced (CPTAC-3, APOLLO), not an established prevalence.

## C4 — Lineage-conditioned signaling
- **Prevalence.** TCGA 1 best (0.9%), CPTAC 0. But contributes partial mass across many mucinous/basal/
  EMT/proliferative tumors.
- **Typical evidence.** Strong non-AT2 dominant lineage × pathway output.
- **Reproducibility.** ➖ Rarely wins outright in either cohort.
- **Interpretation.** We read the thinness as a real result: in this data, lineage state behaves as a
  **modifier** of other mechanisms, not a standalone driver. It rides along with C2/C6/C7 rather than
  dominating. (A single-cell / lineage-resolved cohort could promote it — see knowledge gaps.)

## C5 — Ligand / autocrine signaling
- **Prevalence.** TCGA 5 (4.3%, mass 0.062), CPTAC 4 (**18.2%**, mass 0.162 — notably higher on protein).
- **Typical evidence.** NRG1/HGF/EGF-family/FGF/IGF/PDGF outlier at RNA or protein level.
- **Reproducibility.** ✅ both, and *stronger* in CPTAC where ligand protein can be measured directly —
  a hint that autocrine ligand mechanisms are under-counted when only RNA is available.
- **Strength / limitation.** Suggestive; a ligand outlier does not prove autocrine signaling without the
  cognate receptor being active. In CPTAC these could be cross-checked against receptor phospho (future).

## C6 — Microenvironment-dependent
- **Prevalence.** TCGA 23 (19.7%, mass 0.195), CPTAC 6 (**27.3%**, mass 0.196). Highly consistent.
- **Typical evidence.** High leukocyte fraction (TCGA GDC) / immune+stromal signature (CPTAC RNA).
- **Reproducibility.** ✅ one of the most stable classes across platforms.
- **Interpretation.** A large fraction of oncogene-negative tumors are best distinguished not by a tumor-
  intrinsic driver but by immune/stromal context. This is a real biological statement *and* a caution:
  a bulk "microenvironment" call can mask a low-purity tumor whose intrinsic driver is simply diluted.

## C7 — High-plasticity / state-entropy
- **Prevalence.** TCGA 32 (27.4%, mass 0.223), CPTAC 5 (22.7%, mass 0.233). Consistent.
- **Typical evidence.** Top-quartile lineage-signature entropy, scored only as a **residual** after
  cell-intrinsic lesions are accounted for; concentrated in EMT/mesenchymal and proliferative tumors.
- **Reproducibility.** ✅ both.
- **Interpretation (revised).** C7 is a **residual / holding class — tumors whose mechanism cannot be
  confidently inferred from currently available public RNA/CNA/protein evidence** — *not* a positive claim
  of transcriptional plasticity. It may contain hidden point mutations, regulatory/non-coding mechanisms,
  weak cooperating alterations, measurement limitations (purity/admixture), or genuinely novel biology; the
  present analysis cannot distinguish these. The strict-cohort calibration (`STRICT_COHORT.md`) shows KRAS
  point mutations — invisible to the feature space — preferentially land here (Fisher p=0.009, OR 11.2),
  which *demonstrates the possibility* of hidden drivers in C7 without proving it for any specific tumor.
  C7 (with C6) is therefore a **high-priority target for future genomic characterisation**, and any
  individual C7 call needs single-cell/WGS to resolve.

## C8 — Unknown
- **Prevalence.** TCGA 1 best, CPTAC 0; small residual mass elsewhere.
- **Interpretation.** Deliberately hard to reach (residual only). Its near-emptiness is *not* a claim
  that every tumor is explained — it reflects that C6/C7 absorb context/ambiguity. The honest "we can't
  explain this tumor" population is captured by the **cell-intrinsic-explanation ranking** in
  Deliverable #4, not by C8 alone.

---

## Cross-mechanism reading

- **Cell-intrinsic vs context.** Summing best-calls: cell-intrinsic driver classes (C1–C5) win in
  **~54% of TCGA and ~45% of CPTAC** oncogene-negative tumors; the rest are context (C6) or plasticity
  (C7). Oncogene-negative LUAD is roughly half "hidden driver of some kind" and half "not a
  tumor-intrinsic driver problem at all."
- **The platform changes the answer.** C3 exists only on protein; C5 and C2 are mechanistically tighter
  on phospho; C1's positive control is cleaner on phospho. The taxonomy is only as complete as the assay
  — the single strongest argument for proteogenomics in this subtype.
- **Reproducible core:** C1, C2, C5, C6, C7 recur across both independent cohorts and platforms. That
  five-class core is the robust skeleton of the taxonomy; C3 extends it wherever proteome exists; C4/C8
  are the thin, caveated edges.
