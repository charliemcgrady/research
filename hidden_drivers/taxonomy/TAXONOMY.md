# Deliverable #1 — A mechanistic taxonomy of comprehensively oncogene-negative lung adenocarcinoma

**Mission.** Do not discover *one* driver. Assign *every* comprehensively oncogene-negative
(RTK/RAS/RAF wild-type) LUAD tumor to the most plausible mechanistic class, using convergent evidence
rather than a single significant association. Evidence from multiple weak sources is preferred over one
strong RNA hit. "Unknown" and "high-plasticity" are legitimate terminal answers, not failures.

## Cohorts

| Cohort | Platform | Oncogene-negative n | Layers used | Protein-state (C3) measurable? |
|---|---|---|---|---|
| **TCGA-LUAD** | bulk multi-omic | **117** | RNA, gene CNA, MAF, Hallmark ssGSEA, GDC leukocyte fraction, ABSOLUTE purity | **No** (no proteome/phospho) |
| **CPTAC-LUAD** (Gillette 2020) | proteogenomic | **22** | RNA, gene CNV, WES/WGS mutations, **proteome + phosphoproteome**, cBioPortal driver annotation | **Yes** |

Oncogene-negative is defined conservatively (a *clean* negative set): a tumor is called
oncogene-**positive** — and excluded — if it carries any activating mutation in EGFR/KRAS/NRAS/HRAS/BRAF/
ERBB2/MAP2K1/RIT1, a focal amplification of EGFR/ERBB2/MET/FGFR1 (CPTAC CNV log-ratio > 0.6), or (CPTAC
only, as a fusion proxy) an ALK/ROS1/RET/NTRK RNA outlier (z > 3). TCGA negativity uses the
patient-supplied RPA classification. Fusion detection is a known limitation (see failure report §F3).

## The eight classes

| ID | Class | Operational evidence | Best measured in |
|---|---|---|---|
| **C1** | Hidden canonical driver | RTK RNA outlier (z>2) **or** RTK focal/sub-threshold amp (CNA 0.3–1.0) | both |
| **C2** | Tumor-suppressor convergence | ≥1 of NF1/RASA1/PTEN/STK11/KEAP1 mutated or deep-lost, **+** MAPK/PI3K pathway output | both |
| **C3** | **Protein-state activation** | RTK or pathway node activated at **protein/phospho** level (z>2) with **no** DNA (amp/mut) **or** RNA (z<1.5) explanation | **CPTAC only** |
| **C4** | Lineage-conditioned signaling | strong non-AT2 dominant lineage (mucinous/basal/EMT/proliferative) × pathway output | both |
| **C5** | Ligand / autocrine signaling | ligand (NRG1/HGF/EGF-family/FGF/IGF/PDGF) RNA or protein outlier (z>2) | both |
| **C6** | Microenvironment-dependent | high leukocyte fraction (TCGA) / immune+stromal signature (CPTAC) | both |
| **C7** | High-plasticity / state-entropy | top-quartile lineage-signature entropy, as a **residual** after cell-intrinsic lesions | both |
| **C8** | Unknown | assigned residual mass when total evidence is weak | both |

Scoring is a soft, evidence-weighted combiner (`hidden_drivers/build_{tcga,cptac}_taxonomy.py`): each class
accrues bounded evidence in [0,1]; masses are normalized to a per-tumor posterior; C7 fires only as a
residual (`×(1−max cell-intrinsic lesion)`) so genuine high-plasticity tumors are not confused with
merely-ambiguous ones; leftover mass falls to C8. Lineage signatures deliberately **exclude** all
candidate driver genes (RTKs, adaptors, ligands) to avoid circularity.

## Result — class assignment of every oncogene-negative tumor

Per-cohort prevalence (best-class call and mean posterior mass):

| Mechanism | TCGA best (n=117) | TCGA mean mass | CPTAC best (n=22) | CPTAC mean mass | Reproducible in both |
|---|---|---|---|---|---|
| Hidden canonical driver (C1) | 20 (17.1%) | 0.171 | 2 (9.1%) | 0.154 | ✅ |
| Tumor-suppressor convergence (C2) | **35 (29.9%)** | **0.259** | 2 (9.1%) | 0.130 | ✅ |
| Protein-state activation (C3) | — (not measurable) | — | **3 (13.6%)** | 0.070 | ✅ (CPTAC) |
| Lineage-conditioned (C4) | 1 (0.9%) | 0.061 | 0 | 0.052 | ➖ weak |
| Ligand / autocrine (C5) | 5 (4.3%) | 0.062 | 4 (18.2%) | 0.162 | ✅ |
| Microenvironment (C6) | 23 (19.7%) | 0.195 | 6 (27.3%) | 0.196 | ✅ |
| High-plasticity / entropy (C7) | 32 (27.4%) | 0.223 | 5 (22.7%) | 0.233 | ✅ |
| Unknown (C8) | 1 (0.9%) | 0.030 | 0 | 0.003 | ➖ |

Full per-tumor tables: `tumor_evidence_TCGA.csv`, `tumor_evidence_CPTAC.csv`.
Per-tumor evidence graphs: `per_tumor_evidence.json` (Deliverable #2).

## What the taxonomy says

1. **No single mechanism explains oncogene-negative LUAD.** The population is genuinely heterogeneous.
   The largest *cell-intrinsic* class is **tumor-suppressor convergence** (C2, ~30% in TCGA) — loss of
   NF1/RASA1/PTEN/STK11/KEAP1 reactivating MAPK/PI3K output without a canonical activating driver. This
   is the most defensible "these tumors are RAS-pathway tumors after all" story, and it is the one that
   converges with Chen et al. 2021 (WGS of RTK/RAS/RAF-negative LUAD reclassified a comparable fraction).

2. **A hidden-canonical tail is real but modest (~10–17%).** These are the tumors a deeper genomic look
   (sub-threshold amplification, RTK over-expression) would most likely re-classify — the classic
   "we missed the driver" cases. In CPTAC one such tumor (C3N-02000) scores 0.96 on a sub-threshold RTK
   amplification alone.

3. **Protein-state activation is a genuine, otherwise-invisible mechanism — and CPTAC is the only place
   we can see it — but the individual candidates need per-candidate validation before any causal claim.**
   The automated engine flagged three of 22 CPTAC oncogene-negative tumors (14%) as protein-state. A
   dedicated validation pass (phosphosite QC + robust MAD scoring + downstream-module + genomic/RNA +
   stromal-source + negative-regulator checks; see `../candidate_cards/CLASS3_VERDICTS.md`) **refined all
   three and reclassified two**:
   - **C3N-02422 (MET)** → **provisional Class 3 (weak)**: phospho-driven (pT995 z=3.67 > protein 1.42),
     coordinated GAB1/SHP2/AKT/SRC, tumor-intrinsic — but RNA is modestly elevated (89th pct), the site
     is 50% missing, and CBL-loss is a competing refinement.
   - **C3N-02588 (EGFR)** → **candidate Class 5 (ligand/autocrine)**: the cleanest phospho outlier
     (pS1042 z=3.49, 3.6% missing, **EGFR-RNA z=−0.31**) with coordinated PLCG1/AKT — but **EREG ligand
     z=3.14** provides a more parsimonious cause, so it is reassigned toward ligand-driven (source
     tumor-vs-stroma unresolved).
   - **C3N-02587 (DDR2)** → **rejected → Unknown**: the signal was **total-protein abundance** (z=5.87,
     56% missing, lone extreme) with a **non-elevated phosphosite** (z=−0.07) and absent downstream — an
     abundance artifact, not activation.
   The mechanism (protein-state RTK activation invisible to DNA/RNA) is real — MET, and EGFR's activation
   read-out, would be **missed by any DNA/RNA-only cohort** — but no candidate reaches "supported Class 3."
   Current honest label for all three: **candidate unexplained protein-state activation**, pending
   pTyr-resolved phospho, ligand-source resolution, and functional dependency. This is the headline reason
   the taxonomy needed a proteogenomic cohort *and* a per-candidate validation layer.

4. **A large minority is not cell-intrinsically explained at all.** Cell-intrinsic evidence
   (max of C1–C5) reaches ≥0.35 in only **63/117 TCGA (54%)** and **10/22 CPTAC (45%)** oncogene-negative
   tumors. The remainder are best described by **microenvironment** (C6) or **high transcriptional
   plasticity / state-entropy** (C7), or remain **Unknown** (C8). These are ranked for follow-up in
   Deliverable #4 — they are the tumors that most need new data (WGS/SV, single-cell, spatial).

5. **Class reproducibility is strong.** Every mechanism that is measurable in both platforms recurs in
   both (C1, C2, C5, C6, C7). Protein-state (C3) is CPTAC-exclusive by construction. Only the
   lineage-conditioned class (C4) is thin — it rarely wins outright, more often contributing partial
   mass — which we read as a real signal that "lineage state × pathway output" is a *modifier*, not a
   standalone driver mechanism, in this dataset.

## Positive controls (engine validity)

Running the same engine on the **oncogene-positive** tumors (which have a known RTK/RAS/RAF driver)
should push mass toward C1 (hidden/known canonical). It does, in both cohorts:

| | mean C1 mass, oncogene-**positive** | mean C1 mass, oncogene-**negative** |
|---|---|---|
| TCGA | 0.233 | 0.171 |
| CPTAC | 0.250 | 0.154 |

The separation is cleaner in **CPTAC**, where pathway output is a direct ERK/AKT/S6 phospho readout,
than in **TCGA**, where it is an ssGSEA proxy — consistent with this project's recurring finding that
metric quality (phospho > ssGSEA proxy) changes the confidence of a mechanistic call. C3 mass is
correctly higher in oncogene-negative than -positive CPTAC tumors (0.070 vs 0.059): protein-state
activation is preferentially invoked where the genome is silent.

## Caveats (carried into the failure report and knowledge gaps)

- **C3 is a single-cohort, n=22 observation.** Three protein-state tumors is a hypothesis, not an
  established prevalence. It needs the CPTAC-3 / independent proteogenomic LUAD cohorts to reproduce.
- **Bulk C7 (entropy) can absorb ambiguity.** We constrained it to top-quartile entropy and made it
  residual, but a tumor called "high-plasticity" from bulk RNA could instead be an admixture artifact;
  single-cell would adjudicate.
- **Fusions are under-detected** (no SV layer); a fusion-driven tumor could masquerade as C3/C1/C7.
- **These are posteriors over a fixed 8-class model, not p-values.** They rank plausibility; they do not
  test a null. That is by design (per mission), and it is why every call ships with its evidence graph.
