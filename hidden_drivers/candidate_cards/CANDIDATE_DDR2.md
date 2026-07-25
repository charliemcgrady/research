# Candidate card — DDR2 (CPTAC C3N-02587, AT2/TRU lineage)

> **Status after validation: REJECTED as Class 3 → reclassified UNKNOWN.** The signal is high **total
> DDR2 protein abundance**, not receptor activation. Multiple stop conditions are met.

Automated first-pass call: C3 protein-state (p=0.58). Validation overturns it. Raw data: `validation_DDR2.json`.

## Check 1 — phosphosite verification (fails)
| Field | Value | Read |
|---|---|---|
| Top DDR2 phosphosite | `DDR2:NP_001014796.1:S461s` | |
| Robust z (MAD) | **−0.07** (42.9th pct) | **NOT an outlier** — the phosphosite is unremarkable |
| Missingness | **0.56** | ⚠️ 56% missing |
| Total-protein robust z | **+5.87** | the entire signal is total protein, and it is a **lone extreme** (2nd-highest tumor z = 3.1) in a feature measured in only 48/109 tumors (56% missing) |
| Activation source | `protein_abundance_driven` | **fails rule #3** — activation follows total protein abundance |

The engine scored this tumor via `max(phospho, protein)` and was captured by the total-protein value. On
phosphorylation — the actual read-out of *activation* — DDR2 is at the cohort median.

## Check 2 — coordinated downstream output (fails)
Only the SHP2 module is (weakly) elevated; SRC-family, ERK, AKT, focal-adhesion modules are **not**.
No coordinated DDR2-consistent signalling. Fails the "coordinated downstream" requirement.

## Check 3 — genomic exclusion
No DDR2 mutation; focal CN = −0.07 (neutral); purity 0.71; no suppressor alt; CNA burden 0.02 (very
quiet genome). Genomically clean — but with no activation to explain, this is moot.

## Check 4 — RNA exclusion
Gene RNA robust z = +0.16 (54th pct) — flat. So DDR2 has **high total protein with flat RNA** (a genuine
post-transcriptional protein *accumulation*), but accumulation ≠ activation without phospho/downstream.

## Check 5 — tumor-intrinsic vs ecosystem (the required DDR2 stromal analysis)
DDR2 is a collagen receptor and can reflect fibroblast biology, so this check was decisive:
- In-tumor CAF signature z = **−2.02**, collagen z = **−3.05** — stroma is **low**, so the DDR2 protein is
  **not** coming from an obvious CAF compartment in this tumor.
- **Collagen ligands are absent**: COL1A1 z = −2.17, COL1A2 = −3.23, COL3A1 = −3.66. A collagen receptor
  with no collagen ligand cannot be collagen-activated.
- Across-cohort corr(DDR2 protein, CAF) = 0.36, corr(DDR2, collagen) = 0.31 — DDR2 protein does partly
  track stroma cohort-wide, consistent with mixed epithelial/mesenchymal origin.
→ Not clearly stromal *in this tumor*, but also not ligand-engaged. The high protein is unexplained
accumulation, not microenvironment-driven signalling.

## Check 6 — negative-regulator loss
Only PTPN1 flagged (weak); no coherent negative-regulator mechanism.

## Check 7 — functional support
Not executed (DepMap not fetchable). Moot given the rejection; a dependency test is not warranted without
activation evidence.

## Check 8 — literature novelty
`lit/DDR2.json` (Asta Paper Finder, 5 papers). Notably: *Phosphoproteomics of collagen receptor networks
reveals SHP-2 phosphorylation downstream of wild-type DDR2* (2013) — SHP2/PTPN11 is the documented DDR2
effector, and SHP2 was indeed the only (weak) module here, but with a **median receptor phosphosite and
absent collagen ligand** that axis is not evidence of activation. *Type I collagen aging impairs
DDR2-mediated tumor-cell growth **suppression*** (2016) — DDR2 can be **growth-suppressive**, so even high
DDR2 protein does not imply oncogenic activation. DDR2 in NSCLC is otherwise known for activating
*mutations* in squamous carcinoma (2007). Net: high wild-type DDR2 protein without phospho/downstream is
not a described LUAD driver — consistent with the abundance-artifact verdict.

## Evidence graph
```mermaid
graph LR
  T["Tumor C3N-02587<br/>AT2/TRU · purity 0.71"]:::t
  M["Candidate: Class-3 (DDR2)"]:::m
  U["Reclassified: UNKNOWN<br/>(unexplained protein accumulation)"]:::u
  PR["DDR2 total protein z=5.87<br/>(lone extreme; 56% missing)"]:::obs
  PH["DDR2 phosphosite z=-0.07 (median)"]:::obs
  R["DDR2 RNA z=+0.16 (flat)"]:::obs
  L["Collagen ligands ABSENT (z=-3)"]:::env
  D["Downstream: only weak SHP2"]:::path
  PR -->|derived_from str=strong ig=DDR2_protein| M
  PH -->|contradicts str=strong ig=DDR2_phospho| M
  D -->|contradicts str=strong ig=downstream| M
  L -->|rules_out ligand str=strong ig=ligand| M
  M -->|reclassified_by| U
  classDef t fill:#032629,color:#faf2e9; classDef m fill:#105257,color:#faf2e9;
  classDef u fill:#7a3b3b,color:#fff; classDef obs fill:#0fcb8c,color:#032629;
  classDef path fill:#f0529c,color:#fff; classDef env fill:#ccc,color:#000;
```

## Verdict
**Rejected as Class 3; reclassified Unknown.** Stop conditions met: activation follows total-protein
abundance; the phosphosite fails (median value, 56% missingness); downstream activation absent; no ligand.
The high DDR2 protein is a real but mechanistically uninterpreted post-transcriptional accumulation.
**Required next test:** pTyr DDR2 (Y740) with pTyr enrichment; if still flat, close as abundance artifact.
This is the clearest demonstration that the CPTAC C3 screen needs per-candidate phospho+downstream
validation — the automated `max(phospho, protein)` score over-called it.
