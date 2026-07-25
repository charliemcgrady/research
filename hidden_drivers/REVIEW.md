# Critical review — `hidden_drivers` mechanistic taxonomy of oncogene-negative LUAD

*Reviewer stance: senior computational oncology, preparing for publication. Not optimized for optimism.
Empirical checks were run against the committed outputs; numbers below are from those checks.*

---

## 1. Executive summary

The project builds a per-tumor, evidence-graph-based mechanistic classification of comprehensively
oncogene-negative LUAD across two cohorts (TCGA bulk n=117; CPTAC proteogenomic n=22), with a
soft 8-class posterior, a strict-cohort sensitivity analysis, and a per-candidate validation of the
protein-state (Class 3) hits.

**The engineering is sound and unusually transparent; the biological claims are not yet publishable as
discovery.** The single most important finding of this review is empirical: **the class assignments are
strongly confounded by tumor purity** (Spearman purity vs C6 ρ=−0.77 TCGA / −0.82 CPTAC; vs C2 +0.59/+0.52;
vs C1 +0.45/+0.48; vs C7 −0.35/ns — all p<0.001 where flagged). A large share of the "mechanistic
landscape" is a restatement of how much tumor is in the sample. Compounding this, **~⅓ of TCGA and ~40% of
CPTAC assignments are statistical near-ties** (best-minus-second posterior margin <0.10), and the posteriors
are explicitly uncalibrated. The project's own earlier lesson — *purity is a first-order confounder in bulk
multi-omics* — was not internalized by the taxonomy engine.

Net: this is a strong **framework + honest-limitations** contribution and a **proteogenomic proof-of-concept
(Class 3)**, not a discovery of new oncogene-negative biology. It becomes publishable if (a) purity is
adjusted and the taxonomy is shown to survive it, (b) the structure replicates in an independent cohort, and
(c) claims are scoped to what bulk data can resolve.

---

## 2. Major strengths

1. **Reproducibility is genuinely good.** Engines are deterministic and saved; independent outputs agree
   100% (engine best_class vs `tcga_tumor_table.csv`; full vs strict on shared tumors). Raw data is
   gitignored but re-downloadable with documented provenance. This is better than most published pipelines.
2. **Evidence transparency.** Every tumor ships a per-item evidence graph with independence groups,
   confounders-checked, and QC status. This is the paper's most novel and defensible asset — falsifiability
   at the level of the individual tumor.
3. **Intellectual honesty under pressure.** The Class-3 validation *overturned two of three* headline
   candidates (DDR2 rejected as an abundance artifact; EGFR reassigned to ligand-driven), and the write-ups
   were corrected. The strict-cohort analysis was re-labeled "retrospective calibration" rather than
   oversold as validation. This is the behavior reviewers reward.
4. **Class 3 is a real methodological idea.** Using phosphoproteomics to surface DNA/RNA-invisible RTK
   activation is legitimate and is the one place the work can see something no DNA/RNA cohort can.
5. **The calibration/blind-spot analysis is a contribution in its own right** — quantifying that
   feature-space-invisible lesions (KRAS point mutations) fall into residual classes is an honest map of
   what bulk multi-omics cannot do.

---

## 3. Major weaknesses

1. **Purity confounding (critical).** As above. C6 is close to a purity readout; C1/C2 are easier to assign
   at high purity; C7 rises at low purity. Without purity adjustment, "mechanism" and "tumor content" are
   conflated. This threatens the core claim that these are *mechanistic* classes.
2. **Heuristic, hand-set scoring.** Every threshold (`nz(z,2,4)`, amp 0.3–1.0, C2's 0.7/0.3 mix, the
   top-quartile entropy gate) is hand-chosen, never fit or cross-validated. The "posterior" is a normalized
   sum of non-independent heuristic scores — a **ranking heuristic, not a probability model**. Sensitivity
   to these constants is untested.
3. **Assignment fragility.** Median best_p ≈0.49; 33% (TCGA) / 41% (CPTAC) of tumors have a <0.10 margin to
   the runner-up. For a large minority the "class" is a coin-flip between two mechanisms.
4. **Tiny proteogenomic cohort.** CPTAC oncogene-negative n=22; after validation, Class 3 is effectively
   **n=1 weak (MET)**. No independent proteogenomic replication.
5. **No external replication of the taxonomy structure.** The two cohorts are TCGA and CPTAC — overlapping
   in era, platform philosophy, and partly in patients' population. There is no third, independent bulk
   cohort showing the same class proportions.
6. **Untested confounders.** Smoking, stage, and TMB were not modeled; batch/platform was not modeled.
   Smoking status in particular co-varies with KRAS, TMB, and the TRU/non-TRU axis and could drive C7/C2.
7. **Circularity risks (contained but present).** Lineage signatures correctly exclude driver genes, but
   C4 and C7 both derive from the same RNA lineage scoring, and C4/C2 share the same pathway-output term —
   so some "distinctness" is normalization competition, not independent biology.

---

## 4. Statistical concerns

- **Posteriors are not calibrated** (acknowledged) and, given §3.3, should be reported everywhere as
  ranking scores. The current calibration test is underpowered (1 discordant visible driver).
- **The KRAS-enrichment test is near-tautological.** It compares point mutations (no RNA/CNA footprint,
  by definition) to amplifications/fusions (footprint by definition) and finds the former in residual
  classes. Fisher p=0.009 / OR 11.2 is real but the 95% CI is 1.5–81 (n=28) and the effect is largely
  built into the feature space. Report it as a *calibration illustration*, not a hypothesis test.
- **7/8 capture rate** rests on n=8 (95% CI 47–99.7%) and the concordance rules, while stated before
  looking at predictions in *this* analysis, were defined after the engine was tuned on the full cohort.
  It is retrospective and imprecise; do not headline it.
- **Multiple-comparison / model-degrees-of-freedom** are unaccounted for across the many hand-set knobs.
- **NRG1 "blind recovery" is expected by construction.** An NRG1 fusion over-expresses NRG1; a class that
  keys on ligand-RNA outliers will capture it. It is a valid *sanity check* (n=1), not a surprising
  demonstration of emergent biology — the earlier "one of the strongest demonstrations" wording overstates it.

---

## 5. Biological concerns

- **Are the classes distinct biology?** Partly. C1 (hidden canonical) and C2 (suppressor convergence) map to
  real, known axes and have positive-control support. C6/C7 are heavily purity/composition-driven. C4 and C8
  are near-empty. So the taxonomy is ~2 well-grounded intrinsic classes + a purity/composition axis + a
  residual, dressed as 8 mechanisms.
- **C2, C6, C7 largely reorganize known programs.** C2 = tumor-suppressor-loss RAS reactivation (Skoulidis;
  Chen 2021 WGS). C6 = the long-known LUAD immune/stromal axis. C7 = the TRU-vs-non-TRU / proliferative /
  EMT transcriptional gradient. None is new biology on its own.
- **Tumor-as-unit is respected in bulk, but bulk cannot separate cell-intrinsic plasticity (C7) from
  admixture** — a limitation the text acknowledges but the classification cannot resolve.
- **Class 3 after validation:** MET is a weak single tumor (non-canonical S/T site, 50% missing, RNA 89th
  pct, CBL-loss competing); DDR2 is an abundance artifact; EGFR is ligand-driven (C5). The mechanism concept
  survives; the specific biology does not yet.

---

## 6. Novelty assessment

- **Not novel as a subtype taxonomy of LUAD.** It does not clearly beat existing molecular/transcriptional
  classifications, and its individual classes recapitulate known programs.
- **Novel as a framework.** The per-tumor, competing-hypothesis, evidence-graph mechanistic *triage* — with
  explicit independence groups and a residual "we can't tell" class — is a genuinely useful and, to my
  knowledge, uncommon way to present oncogene-negative heterogeneity. This is the defensible novelty.
- **Novel as a proteogenomic proof-of-concept.** "Screen for DNA/RNA-invisible protein-state RTK activation"
  is a real idea; the validation methodology (phospho + downstream-module + stromal-source + negative-
  regulator, receptor excluded) is a reusable template even though the current hits are weak.
- **Novel as an honest limits-of-bulk analysis.** The demonstration that ~half of oncogene-negative tumors
  are not cell-intrinsically resolvable from bulk, and that point-mutation drivers hide in the residual, is
  a publishable *negative/limitations* result.

---

## 7. Publication potential

- **Nature Cancer / Cancer Discovery (discovery paper): not competitive as-is.** No new validated driver,
  severe purity confound, n=22 proteogenomics, no independent replication.
- **ICML/NeurIPS (method paper): weak fit.** The engine is heuristic with no learning, benchmark, or
  theoretical contribution. Would need to be reframed as a learned/calibrated probabilistic model with
  held-out evaluation.
- **Realistic strong fit:** a **methods/resource paper** (e.g. *Cell Reports Methods*, *Genome Medicine*,
  *npj Precision Oncology*, *Bioinformatics*) titled around: *a reproducible, evidence-transparent
  mechanistic-triage framework for oncogene-negative LUAD, and an honest accounting of what bulk multi-omics
  can and cannot resolve* — with Class 3 as a proteogenomic proof-of-concept and the calibration analysis as
  a headline figure. This is achievable **after** the purity fix + one external replication.

---

## 8. Highest-priority next analyses (ranked by expected information gain, not effort)

| # | Analysis | Info gain | P(changes conclusions) | Feasible on public data? | Needs controlled/wet-lab? |
|---|---|---|---|---|---|
| 1 | **Purity adjustment / deconvolution, then re-run the whole taxonomy.** Residualize RNA-derived evidence on ABSOLUTE/TSNet purity (or use tumor-cell-deconvolved expression); recompute classes. | **Highest** — tests whether the taxonomy survives its own dominant confounder | **High** (could collapse C6, shrink C1/C2, reshuffle ~⅓ near-ties) | ✅ immediate | no |
| 2 | **External replication of the class structure** on an independent bulk LUAD cohort (e.g. OncoSG/Chen 2020, or a GEO LUAD with driver annotation). Predefine the class proportions and test concordance. | High — generalization is currently unproven | Medium | ✅ (some controlled) | partly |
| 3 | **Sensitivity/stability analysis of the heuristic knobs** (perturb every threshold/weight; bootstrap tumors; report assignment stability and a per-tumor confidence). | High — converts "posterior" into a defensible ranking with stated stability | Medium | ✅ immediate | no |
| 4 | **Model smoking/stage/TMB/batch** as covariates; test whether any class is explained by them. | Medium-high | Medium | ✅ immediate (CPTAC .tsi has smoking/stage; TCGA clinical) | no |
| 5 | **Second proteogenomic cohort (CPTAC-3 / APOLLO)** to test Class 3 reproducibility. | High for C3 only | Medium | partly | controlled-access |
| 6 | **WGS+SV on the C6/C7 residual** to test the hidden-driver hypothesis directly. | High but expensive | Medium | ❌ | controlled-access |
| 7 | **Functional dependency (DepMap/CRISPR/drug) for MET & EGFR/EREG.** | High causal value | Medium | DepMap ✅ (was blocked in-env); wet-lab for EREG source | DepMap public; source needs single-cell/wet-lab |

Tiers:
- **Immediate (public, no new data): #1, #3, #4, and the DepMap portion of #7.** These should be done before
  any submission.
- **Additional public/lightly-controlled data: #2, #5.**
- **Controlled-access genomics: #6, part of #5.**
- **Experimental validation: EREG-source and functional confirmation for #7.**

---

## 9. Highest-risk assumptions

1. **That the classes are mechanistic rather than purity/composition strata.** (Directly at risk; test #1.)
2. **That bulk lineage entropy (C7) reflects tumor-cell plasticity** rather than admixture/purity.
3. **That the hand-set thresholds are near-optimal.** Untested; test #3.
4. **That TCGA+CPTAC generalize to LUAD at large.** No external cohort yet; test #2.
5. **That Class 3 protein-state activation is a real recurrent mechanism** — currently n=1 weak.
6. **That "oncogene-negative" is clean.** 28/118 TCGA tumors were WGS-reclassified; CPTAC has no SV layer, so
   fusions are undetected. The negative set is porous.

---

## 10. Concrete recommendations for the next phase

1. **Do the purity fix first (#1).** Re-run and report the taxonomy on purity-residualized evidence.
   Explicitly reframe C6 as, in part, a *low-purity / composition* flag, and re-examine whether low-purity
   C6 tumors hide an intrinsic driver (this also strengthens Deliverable #4). This is the make-or-break
   analysis; everything else is contingent on it.
2. **Add a stability layer (#3).** Ship a per-tumor confidence = assignment stability under threshold
   perturbation + bootstrap. Demote near-tie tumors (~⅓) to "ambiguous" explicitly.
3. **Replicate the structure externally (#2)** before claiming the taxonomy generalizes.
4. **Reframe the paper as framework + limits-of-bulk + proteogenomic proof-of-concept**, not as a discovery
   of new oncogene-negative biology. Make the calibration/blind-spot analysis and the evidence-graph the
   headline; present Class 3 honestly as a template with one weak lead.
5. **Scope every class claim.** Merge C4 into "lineage modifier" (drop as a standalone class); treat C8 and
   C4 as non-classes; present C6/C7 as *residual uncertainty*, not mechanisms.
6. **Run DepMap for MET/EGFR now** (it was only blocked by the environment, not by data access) — cheap,
   and materially changes the Class-3 story either way.

---

## Direct answers to the specific questions

**Taxonomy — publishable?** As a framework/resource with the purity fix and one external replication: yes,
in a methods/resource venue. As a discovery taxonomy in a top oncology journal: no, not yet.
**Sufficiently different from prior work?** The *framework* is; the *classes* mostly are not.
**Strongest classes:** C1 (positive-control-supported) and C2 (known but well-supported, tighter on CPTAC
phospho). **Weakest:** C3 (n=1 after validation), C4 (near-absent), C8 (empty). **Merge/remove:** drop C4 as
standalone (fold into a lineage modifier); treat C8 as bookkeeping; **relabel C6/C7 as residual/composition
strata**, not mechanisms, until purity-adjusted.

**Class 3.**
- *MET* — supports a **weak** mechanistic claim; keep as candidate. Greatest single confidence gain:
  **pTyr-enriched phosphoproteomics** (confirm Y1234/5) — or, if unavailable, **MET functional dependency
  (DepMap/CRISPR)**.
- *EGFR* — does **not** support an intrinsic protein-state claim; it is a **ligand (EREG) candidate (C5)**.
  Greatest gain: **single-cell/spatial to localize EREG source** (tumor vs stroma).
- *DDR2* — does **not** support a claim; **reject/Unknown**. Greatest gain: **pTyr DDR2 (Y740)**; if flat,
  close it.

**C6/C7 interpretation.** Neither "unexplained biology" nor "hidden-driver candidates." Best described as
**residual uncertainty that is partly purity/composition-driven**. Before any stronger claim, they require:
(a) purity adjustment (#1), (b) single-cell to separate plasticity from admixture, and (c) WGS+SV to test
for hidden lesions. The KRAS calibration shows they *can* hide drivers; it does not show that they *do*.

**Overall direction.**
1. *Three more months, public data only:* purity-adjust + stability + external replication (#1–#4) and the
   DepMap runs — i.e., prove the taxonomy survives its confounders and generalizes, and firm up MET. That,
   plus honest reframing, is a submittable methods/resource paper.
2. *If controlled-access opens:* WGS+SV on the residual (test whether C6/C7 hide drivers — the real
   biological question) and a second proteogenomic cohort (make Class 3 real). This is what could upgrade
   the work toward a genuine oncology-journal contribution.
3. *One wet-lab collaboration:* functional dependency + drug response in RTK/RAS/RAF-WT LUAD models
   stratified by class — specifically test MET dependency in the protein-state/CBL-low state and EGFR/EREG
   autocrine dependency. Causal evidence is the missing layer everywhere.
4. *Central contribution for a paper:* **not** a new driver. It is *"a reproducible, evidence-transparent
   framework that assigns each oncogene-negative LUAD tumor a most-plausible mechanism with explicit
   competing hypotheses and confidence, and a rigorous accounting of the fraction (~half) that bulk
   multi-omics cannot resolve — with proteogenomics as the route to the invisible (protein-state) subset."*
   The scientific value is the honest triage + limits map, not a discovery.
