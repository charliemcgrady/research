# Dataset-Access Report — Lineage-Conditioned Hidden Drivers in Oncogene-Negative LUAD

**Mission step 1 deliverable.** Assesses what is obtainable by a *credential-free* cloud environment
(no dbGaP/EGA) before any discovery. Access legend: **OPEN** = processed matrices downloadable now;
**PROC** = processed/supplement open, raw controlled; **CTRL** = controlled (dbGaP/EGA/national archive).
Status: 5/6 dataset families assessed; CPTAC multi-omic linkage pending.

---

## Executive access verdict (drives the feasibility gate)

1. **TCGA-LUAD processed data is OPEN and rich** (mutations, GISTIC/ASCAT CN, RNA, RNA-fusions,
   methylation, ABSOLUTE purity, clinical incl. smoking). **But WGS BAMs, DNA-breakpoint SVs,
   splice junctions, and allele-specific/non-coding resolution are CONTROLLED (dbGaP phs000178).**
2. **The TCGA oncogene-negative WGS "hidden driver" study already exists:** Chen et al. 2021,
   *Cell Reports* — WGS on 85 of ~118 TCGA pan-negatives; 28/85 reclassified as driver-positive on
   deeper data; residual hidden drivers = **STK11/KEAP1 promoter/TSS focal deletions, ILF2 promoter
   non-coding mutations, and complex amplicons (ecDNA/double-minutes, BFB)**. Processed calls are OPEN
   (supplement); raw WGS controlled. **This pre-empts much of the structural-driver aim and raises the
   novelty bar.**
3. **No open cohort has ≥100 LUAD with WGS-derived SV + matched RNA.** PCAWG has open consensus SV/CNV
   but only **~38** Lung-AdenoCA; every scaled WGS+RNA cohort (Sherlock-Lung ~871, NCC-Japan ~938,
   Korean n=99) is CTRL or just-under-threshold. **→ Phase 3 (Exp 11–20: local SVs, enhancer hijacking,
   cryptic splice, amplicon geometry) is NOT runnable at scale on open data.**
4. **Regulatory tracks are OPEN but contextual only** (A549 tumor line / IMR90 fibroblast, much hg19) —
   corroborative, never tumor-specific.
5. **DepMap + GDSC/CTRP/PRISM are OPEN** → functional dependency (Phase 6) is feasible, but only
   **~8–15 plausibly driver-negative LUAD lines** with full CRISPR+expression+drug (clears the ≥5 rule,
   thin).
6. **Single-cell lineage model is buildable** from LuCA (in hand) + HLCA/Travaglini normal refs, with
   Kim/Laughney/Bischoff as open validators — **except mucinous / gastric-like(HNF4A): no open human
   scRNA exists** (that lineage arm cannot be trained/validated on open data).
7. **CPTAC (pending)** determines feasibility of the protein/phospho-state outlier experiments
   (Phase 2, Exp 9–10) — the most promising *still-open* novel direction.

---

## Preliminary driver-negative census — TCGA-LUAD (local data)

| Definition | n | Basis |
|---|---:|---|
| RPA-classified total | 501 | published RTK/RAS/RAF call (Chen 2021 companion) |
| Driver-positive | 383 | |
| RPA-negative | 118 | |
| DN-Clinical (+ no clinical activating mut in MAF) | 86 | EGFR/KRAS/ALK/ROS1/RET/BRAF/MET-ex14/ERBB2/NTRK |
| DN-Extended (+ no extended RAS-pathway mut) | 75 | + RIT1/ARAF/MAP2K1/HRAS/NRAS/RAF1/FGFR/NRG1/SOS1 |
| **DN-Strict** (purity ≥0.3 + RNA + MAF) | **48** | 27 with WGS; median purity 0.55 |

Suppressor annotation within DN-Extended (n=75): TP53 48%, KEAP1 24%, STK11 17%, SMARCA4 9%, RB1 5%,
PTEN 1%, NF1/RASA1 0% (MAF-only; deep deletions not yet added). Consistent with Chen 2021's ~57 true
pan-negative after WGS reclassification. **Caveats:** fusions/SV beyond MET-ex14 rely on the published
RPA call (not re-called locally); smoking not yet joined (GDC clinical); purity floor 0.3 provisional.
Table: `driver_audit_preliminary.csv`.

---

## Per-family access summary

### TCGA-LUAD — OPEN (processed) / CTRL (WGS, SV, junctions)
Best one-stop OPEN matrix: cBioPortal `luad_tcga_pan_can_atlas_2018` (566; mut, GISTIC CN, RSEM RNA,
HM450, RPPA, fusion-SVs, ABSOLUTE clinical, TCGA-CDR survival, smoking). hg38 primaries via GDC; Xena
for pre-built matrices. Open fusions: GDC Arriba/STAR-Fusion, Gao 2018 PanCanAtlas, TumorFusions,
ChimerDB. Open junctions only via **recount3/Snaptron** (reprocessed). ~500 tumors meet DN-Strict
profiling; oncogene-negative subset ~48–57. **Blocks:** genome-wide SV/enhancer-hijacking, splice,
allele-specific/non-coding (all need controlled WGS/BAM).

### Large WGS LUAD cohorts — none OPEN at ≥100 WGS+SV+RNA
PCAWG Lung-AdenoCA ~38 (OPEN consensus SV/CNV/driver, no matched open RNA). OncoSG 305 & TRACERx ~248
are **WES** (RNA-fusion "SV" only; OPEN processed via cBioPortal/supplements). Sherlock-Lung ~871 and
NCC-Japan ~938 = **CTRL** (dbGaP/national archive). Korean EGFR/ALK-WT never-smoker **n=99** — oncogene-
negative *by design*, WGS-SV + RNA in OPEN supplements (but <100, no breakpoint matrix; PXD033360
proteome open). **For real structural work at scale → dbGaP phs001697/phs002346 (Sherlock) or EGA
TRACERx required.**

### Single-cell / spatial — OPEN (lineage substrate), mucinous gap
LuCA core/extended (in hand; CELLxGENE `edb893ee…`) — malignant cells annotated coarsely, **lineage
states NOT in obs → re-derive** (markers + NMF). Normal refs: **HLCA (Sikkema 2023), Travaglini 2020**
(OPEN CELLxGENE). Validators (OPEN): Kim GSE131907, Laughney GSE123904, Bischoff (Code Ocean). NKX2-1-
high→low axis explicit in Maynard PRJNA591860 + AAH→IA cohort. Spatial (OPEN): GSE189357/GSE189487
(Visium+scRNA), Nat Commun 2024 Visium+Xenium. **Gap: mucinous / gastric-like(HNF4A) — no open human
scRNA** (only mouse GEMM + one Visium IMA); this lineage arm is unvalidatable on open data.

### DepMap + drug response — OPEN
DepMap 24Q4 (Figshare CC-BY): `CRISPRGeneEffect.csv`, `OmicsExpressionProteinCodingGenesTPMLogp1.csv`,
`OmicsSomaticMutations.csv`, `OmicsFusionFiltered.csv`, `OmicsCNGene.csv`, `Model.csv` (filter
`OncotreeCode==LUAD` → ~80 lines). Drug: GDSC2 (`GDSC2_fitted_dose_response_27Oct23.xlsx` +
`screened_compounds_rel_8.5.csv` target/pathway), CTRP v2 (`CTRPv2.0_2015_ctd2_ExpandedDataset.zip`),
PRISM 24Q2 (Figshare). ~**8–15 driver-negative LUAD lines** with full CRISPR+expr+drug after audit
(clears ≥5-line rule; NF1/RIT1/fusion mis-calls are the main risk). Clinical-phase via Drug Repurposing
Hub.

### Regulatory annotations — OPEN, contextual only
ENCODE ATAC/DNase/H3K27ac (A549, SAEC, lung tissue), SCREEN cCREs (GRCh38), dbSUPER/SEdb super-enhancers,
Roadmap ChromHMM (E096 Lung, hg19), A549/IMR90 Hi-C (TADs/loops), GTEx lung expr+eQTL (open summary),
gnomAD-SV v4 + DGV (germline-SV filter), ENCODE blacklist v2 + Umap/Bismap (artifact masks). **Caveats:**
hg19↔GRCh38 liftover; A549(tumor)/IMR90(fibroblast) proxies, no adult-normal-AT2 chromatin panel;
strictly corroborative per mission rule.

*(CPTAC section pending — determines Phase 2 protein/phospho feasibility.)*
