# Candidate card — EGFR (CPTAC C3N-02588, Ciliated lineage)

> **Status after validation: reassigned toward Class 5 (ligand / autocrine), with Class-3 protein-state
> activation as the read-out.** The receptor *is* activated at protein level with silent DNA/RNA — but a
> convincing ligand (EREG) provides the more parsimonious upstream cause. Candidate, not established.

Automated first-pass call: C3 protein-state (p=0.32; C5 ligand was already the 3rd competitor at 0.21).
Raw data: `validation_EGFR.json`.

## Check 1 — phosphosite verification (cleanest of the three)
| Field | Value | Read |
|---|---|---|
| Top phosphosite | `EGFR:NP_005219.2:S1042s` | C-terminal regulatory S/T site (canonical pY not in panel) |
| Robust z (MAD) | **+3.49** (99.1st pct) | strong outlier |
| Missingness | **0.036** | ✅ excellent — 3.6% missing |
| Total-protein robust z | **−0.07** | protein at median → **phospho-driven**, not abundance |
| Activation source | `phospho_driven` | passes rule #3 cleanly |

Best data quality of the three candidates (low missingness, phospho-specific, protein not elevated).

## Check 2 — coordinated downstream output (EGFR sites excluded)
Active modules (robust z > 1.5): **AKT, adaptor (SHC1/GAB1/GRB2), PLCG1** → `coordinated = True`. PLCG1 is
a **direct EGFR effector**, which makes the coordinated-output evidence here stronger and more receptor-
specific than for MET.

## Check 3 — genomic exclusion
No EGFR mutation; focal CN +0.15 (neutral); purity **0.89** (high); no suppressor/RAS-pathway alt; CNA
burden 0.03 (quiet genome). Arm-CN/WGD/fusion not available (untested). Clean at available resolution.

## Check 4 — RNA exclusion + the ligand finding
- EGFR gene RNA robust z = **−0.31 (36th percentile — below median).** Strong receptor-RNA discordance —
  the protein-state read-out is genuine.
- **Ligand outliers:** **EREG (epiregulin) z = +3.14**, EGF +1.77, AREG +1.66, BTC +1.40. Multiple EGFR
  ligands are elevated, EREG strongly.
→ Per the directive's rule ("elevated ligand expression may instead support Class 5 … reassign when
receptor activation is accompanied by convincing ligand evidence"), EREG is **convincing ligand evidence**.
The mechanism is most parsimoniously **ligand-driven (autocrine/paracrine) EGFR activation**, with the
high pEGFR/low-EGFR-RNA pattern as the expected read-out of a ligand-activated, internalised/recycled
receptor.

## Check 5 — tumor-intrinsic vs ecosystem (ligand source)
Purity 0.89; CAF z = −0.18, collagen z = −0.20 (stroma low); corr(EGFR activation, CAF) ≈ 0. So the
receptor activation is tumor-intrinsic. **Ligand source is unresolved:** EREG–CAF corr = 0.38 and
EREG–EPCAM corr = −0.10 across the cohort, so EREG is **not confirmed tumor-derived** vs stromal. This is
the key open question separating autocrine (tumor EREG → Class 5 intrinsic) from paracrine
(stromal EREG → Class 6 microenvironment).

## Check 6 — negative-regulator loss
Flagged: PTPN1, SPRY2, SPRY4, DUSP6. SPRY/DUSP are ERK-feedback regulators normally *induced* by active
ERK, so low values here are ambiguous (could reflect blunted feedback or noise); not treated as a primary
mechanism. No strong CBL/ERRFI1 loss.

## Check 7 — functional support
Not executed (DepMap not fetchable). `required_next_test`: EGFR and EREG dependency (CRISPR + EGFR-TKI /
anti-EREG) in RTK/RAS/RAF-WT LUAD lines with high-EREG / high-pEGFR / low-EGFR-RNA state.

## Check 8 — literature novelty
`lit/EGFR.json`. EREG/AREG autocrine EGFR signalling is described in colorectal and some NSCLC contexts and
is a known cetuximab-response biomarker; wild-type EGFR activation without mutation/amplification is
reported. The specific claim to test for novelty is the proteogenomic pattern in **oncogene-negative
LUAD**: pEGFR-high / EGFR-RNA-low / EREG-high with coordinated PLCG1-AKT output.

## Evidence graph
```mermaid
graph LR
  T["Tumor C3N-02588<br/>Ciliated · purity 0.89"]:::t
  C3["Class-3 read-out<br/>(protein-state EGFR)"]:::m
  C5["Reassigned: Class-5<br/>ligand/autocrine (EREG)"]:::m5
  PH["EGFR pS1042 z=3.49 (99th pct, 3.6% missing)"]:::obs
  PRc["EGFR protein z=-0.07 · RNA z=-0.31"]:::obs
  D["Downstream AKT/adaptor/PLCG1 active"]:::path
  LIG["EREG z=3.14 (+EGF/AREG/BTC)"]:::env
  SRC["EREG source tumor-vs-stroma UNRESOLVED<br/>(EREG–CAF r=0.38)"]:::env
  PH -->|supports str=strong ig=EGFR_phospho| C3
  D -->|supports str=moderate ig=downstream| C3
  PRc -->|supports str=strong ig=EGFR_abundance| C3
  LIG -->|may_explain str=strong ig=ligand| C5
  C3 -->|refined_by| C5
  SRC -->|requires_test str=strong ig=ligand_source| C5
  classDef t fill:#032629,color:#faf2e9; classDef m fill:#105257,color:#faf2e9;
  classDef m5 fill:#f0529c,color:#fff; classDef obs fill:#0fcb8c,color:#032629;
  classDef path fill:#105257,color:#faf2e9; classDef env fill:#ccc,color:#000;
```

## Verdict
**Candidate Class 5 (ligand/autocrine EGFR), Class-3 protein-state as read-out.** The receptor activation
is real, high-quality, tumor-intrinsic, and DNA/RNA-silent — but EREG elevation is the more parsimonious
cause, so it does **not** qualify as an intrinsic protein-state driver. It is neither pure C3 nor confirmed
C5 until the **ligand source** (tumor vs stroma) is resolved. **Required next tests:** single-cell/spatial
to localise EREG; pEGFR-pY confirmation; EGFR/EREG functional dependency.
