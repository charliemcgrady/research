# Candidate card — MET (CPTAC C3N-02422, Basal lineage)

> **Status after validation: PROVISIONAL Class 3 (weak, caveated).** Candidate unexplained
> protein-state activation — **not** established as a post-transcriptional driver. Do not label causal.

Automated first-pass call: C3 protein-state (p=0.42). The validation checks below both support and qualify
that call. Raw check data: `validation_MET.json`.

## Check 1 — phosphosite verification
| Field | Value | Read |
|---|---|---|
| Top phosphosite | `MET:NP_001120972.1:T995t` | **Ser/Thr, not the activation-loop tyrosine** (Y1234/Y1235). |
| Robust z (MAD, all 109 tumors) | **+3.67** (98.2nd pct) | strong outlier |
| Missingness | **0.50** | ⚠️ 50% of tumors missing this site — a real QC concern |
| Lineage-matched robust z | n/a | too few Basal peers to compute |
| Total-protein robust z | +1.42 | phospho exceeds protein → **phospho-driven**, not abundance-driven |
| Activation source | `phospho_driven` | passes "not explained by total protein" |

**Caveat (applies to all three candidates):** global phosphoproteomics under-samples phosphotyrosine, so
the canonical MET activation-loop pY (Y1234/Y1235) is absent from the panel. The outlier sits on a
juxtamembrane/C-terminal S/T site of uncertain activating function. The activation claim therefore leans
on coordinated downstream signalling, not on a validated activating residue.

## Check 2 — coordinated downstream output (MET sites excluded)
Active modules (robust z > 1.5): **adaptor (GAB1), SHP2 (PTPN11), AKT, SRC** → `coordinated = True`.
Caveat: GAB1/SHP2/SRC are shared by many RTKs, so this supports "an RTK program is active," not "MET
specifically." STAT3 and ERK modules were not elevated.

## Check 3 — genomic exclusion
No MET mutation; focal CN = +0.26 (sub-threshold, not amplified); purity **0.85** (high — not a low-purity
artifact); no RTK/RAS/RAF or suppressor alteration; total CNA burden 0.29; **arm-level CN, WGD, and
fusion/SV not available from open calls** (study exposes 0 SVs). "No amplification ≠ no genomic
explanation" — a MET fusion/rearrangement cannot be excluded here.

## Check 4 — RNA exclusion
Gene RNA robust z = **+1.54 (89th percentile)** — MET RNA is **modestly elevated**, so RNA discordance is
only partial; this weakens the pure protein-state reading (borderline against rule #4). HGF ligand RNA
z = +0.07 (no autocrine ligand).

## Check 5 — tumor-intrinsic vs ecosystem
Purity 0.85; CAF signature z = −0.70, collagen z = −0.93 (stroma **low**); corr(activation, CAF) = 0.10,
corr(activation, collagen) = 0.09. → **tumor-intrinsic, not stromal.** Passes.

## Check 6 — negative-regulator loss
Flagged low: **CBL**, PTPN1, SPRY2, DUSP6. **CBL is MET's E3 ubiquitin ligase** — CBL loss stabilises
active MET and would explain protein-state activation without DNA/RNA change. This is recorded as a
**competing / refining mechanism** (protein-state activation *via* negative-regulator loss).

## Check 7 — functional support (DepMap/GDSC/PRISM)
**Not executed** (DepMap CRISPR/drug matrices not fetchable in this environment). Specified as
`required_next_test`: MET dependency (CRISPR) and MET-inhibitor sensitivity in RTK/RAS/RAF-WT LUAD lines
with matched high-MET-protein / low-CBL state, controlling for lineage.

## Check 8 — literature novelty
`lit/MET.json` (Asta Paper Finder, 30 papers). Closest prior observations: **HGF–MET autocrine loop**
enhances tumorigenicity in LUAD (2000); **constitutive c-Met activation correlated with overexpression and
dependent on cell–matrix adhesion** (2007); KIF5B–MET fusion (2018). So non-mutational MET activation
(autocrine, adhesion-dependent, overexpression, fusion) is established in principle. The specific untested
claim here is a **oncogene-negative LUAD tumor with no amp / no exon-14 / only modestly elevated RNA, where
MET activation is read out proteomically with coordinated downstream and candidate CBL loss** — not an
established driver mechanism, hence "candidate."

## Evidence graph
```mermaid
graph LR
  T["Tumor C3N-02422<br/>Basal · purity 0.85"]:::t
  M["Candidate: Class-3 protein-state (MET)"]:::m
  T -->|may_explain str=weak| M
  P["MET pT995 robust z=3.67<br/>(S/T site; pY not captured; 50% missing)"]:::obs
  PR["MET total protein z=1.42"]:::obs
  R["MET RNA z=+1.54 (89th pct)"]:::obs
  CN["MET focal CN +0.26; no mut; no fusion(untested)"]:::gen
  D["Downstream GAB1/SHP2/AKT/SRC active"]:::path
  CBL["CBL low (E3 ligase) + SPRY2/DUSP6 low"]:::path
  S["Stroma low (CAF z=-0.7)"]:::env
  P -->|supports str=moderate ig=MET_phospho| M
  D -->|supports str=moderate ig=downstream| M
  CBL -->|may_explain str=weak ig=negreg| M
  R -->|contradicts str=weak ig=RNA| M
  CN -->|requires_test str=strong ig=genomic| M
  S -->|rules_out stromal str=moderate ig=microenv| M
  classDef t fill:#032629,color:#faf2e9; classDef m fill:#105257,color:#faf2e9;
  classDef obs fill:#0fcb8c,color:#032629; classDef path fill:#f0529c,color:#fff;
  classDef gen fill:#888,color:#fff; classDef env fill:#ccc,color:#000;
```
Edges are `may_explain`/`supports`/`contradicts`/`requires_test` with independence groups (ig); no causal
edge is asserted from association.

## Verdict
**Provisional Class 3, weak.** Meets the provisional criteria (no canonical/genomic driver detected;
phospho-driven, not protein-abundance-driven; coordinated downstream; QC-passable on purity; tumor-
intrinsic) **but** does not upgrade to supported Class 3: RNA is modestly elevated, the outlier site is a
non-canonical S/T with 50% missingness, and MET fusion is untested. CBL-loss is a plausible competing
refinement. **Required next tests:** pTyr-enriched phospho (Y1234/5), MET fusion screen (WGS/SV), CBL
protein + functional MET dependency.
