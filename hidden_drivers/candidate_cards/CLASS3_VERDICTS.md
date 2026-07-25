# Class-3 candidate validation — verdicts (MET, DDR2, EGFR)

## Provisional summary (required wording)

> MET, DDR2, and EGFR are candidate Class-3 protein-state activation events. Each shows receptor-level
> protein or phosphoprotein activation without obvious amplification and without corresponding RNA
> overexpression. EGFR is particularly discordant, with RNA near the cohort mean despite phosphoprotein
> activation. These findings justify pathway-level and tumor-source validation but do not yet establish
> post-transcriptional driver mechanisms.

Per-candidate validation (checks 1–8, robust median/MAD scoring across 109 CPTAC-LUAD tumors) refines the
automated first-pass calls. **None upgrades to "supported Class 3"; two are reclassified.** This is the
intended function of the validation layer — the automated `max(phospho, protein)` C3 score over-calls, and
per-candidate phospho + downstream + source checks are required before any causal language.

## Verdicts

| Candidate | Tumor | Automated | **After validation** | Decisive evidence |
|---|---|---|---|---|
| **MET** | C3N-02422 | C3 (0.42) | **Provisional Class 3 (weak)** | phospho-driven (pT995 z=3.67 > protein 1.42); coordinated GAB1/SHP2/AKT/SRC; tumor-intrinsic (CAF low, purity 0.85); no genomic driver. **Caveats:** RNA modestly ↑ (89th pct); 50% site missingness; S/T site (pY not captured); CBL/SPRY2/DUSP6 low → CBL-loss competing refinement; MET fusion untested. |
| **DDR2** | C3N-02587 | C3 (0.58) | **REJECTED → Unknown** | signal is total-protein abundance (z=5.87, lone extreme, 56% missing), **phosphosite not elevated (z=−0.07)**; downstream absent (only weak SHP2); collagen ligand absent (z=−3). Fails stop conditions "activation follows total protein" + "downstream absent." |
| **EGFR** | C3N-02588 | C3 (0.32) | **Candidate Class 5 (ligand/autocrine); C3 read-out** | cleanest phospho outlier (pS1042 z=3.49, 99th pct, 3.6% missing, phospho-driven, EGFR-RNA z=−0.31); coordinated AKT/adaptor/**PLCG1**; tumor-intrinsic — **but EREG ligand z=3.14** (+EGF/AREG/BTC) ⇒ reassign per rules. Ligand source (tumor vs stroma) unresolved (EREG–CAF r=0.38). |

Cards: [`CANDIDATE_MET.md`](CANDIDATE_MET.md), [`CANDIDATE_DDR2.md`](CANDIDATE_DDR2.md),
[`CANDIDATE_EGFR.md`](CANDIDATE_EGFR.md). Raw checks: `validation_{MET,DDR2,EGFR}.json`.
Machine-readable: `class3_verdicts.csv`.

## Classification-rule accounting

Provisional Class 3 requires all of: (1) no canonical driver, (2) robust protein-state outlier, (3) not
explained by total protein, (4) not explained by RNA, (5) coordinated downstream, (6) QC pass.

| Rule | MET | DDR2 | EGFR |
|---|---|---|---|
| 1 no canonical driver | ✅ | ✅ | ✅ |
| 2 robust outlier | ✅ (z=3.67) | ⚠️ protein only | ✅ (z=3.49) |
| 3 not total-protein-driven | ✅ | ❌ (z=5.87) | ✅ |
| 4 not RNA-driven | ⚠️ RNA 89th pct | ✅ | ✅ (RNA −0.31) |
| 5 coordinated downstream | ✅ | ❌ | ✅ |
| 6 QC pass | ⚠️ 50% missing | ❌ 56% missing | ✅ 3.6% missing |
| **Provisional C3?** | **yes (weak)** | **no → Unknown** | **yes, but ligand → C5** |

No candidate meets an **upgrade** criterion (independent recurrence, negative-regulator mechanism proven,
ligand-independent evidence, functional dependency, drug sensitivity, or independent-cohort support) —
functional (check 7, DepMap) was **not executable in-environment** and is the top `required_next_test` for
all three; a second proteogenomic cohort (knowledge gap G2) is required for independent replication.

## Literature novelty (check 8, Asta Paper Finder)

- **EGFR** (`lit/EGFR.json`, 38 papers): closest prior art centres on **MIG6/ERRFI1 and CBL feedback loss**
  and EREG/AREG autocrine EGFR signalling (e.g. *Mig6/EGFR expression & TKI resistance*, 2013; *c-Cbl
  induction in lung cancer*, 2015; MIG6-loss tumorigenesis, 2024). Wild-type ligand-driven EGFR activation
  is described; the untested novelty is the specific oncogene-negative-LUAD proteogenomic pattern
  (pEGFR↑/EGFR-RNA↓/EREG↑ + PLCG1-AKT).
- **MET / DDR2**: searches `lit/MET.json`, `lit/DDR2.json` (see files). MET is canonically amp/exon-14
  driven; DDR2 is a squamous-mutation / CAF-biology gene — neither has an established "wild-type
  protein-state activation in oncogene-negative LUAD" precedent, which is exactly why these remain
  candidate and require functional support before any novelty claim.

## Bottom line
The Class-3 mechanism survives as a **real, DNA/RNA-invisible phenomenon** (MET, and EGFR's activation
read-out), but the specific candidates are individually fragile: one is an abundance artifact (DDR2), one
is more likely ligand-driven (EGFR), and the strongest (MET) is caveated by RNA and QC. The honest current
status of all three is **candidate unexplained protein-state activation**, pending pTyr-resolved phospho,
ligand-source resolution, and functional dependency.
