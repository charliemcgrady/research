# Evidence Catalog — public data to test the per-cell LAG3 finding

**Finding to validate:** In oncogene/driver-negative (RTK/RAS/RAF wild-type) LUAD, chromosomal
instability (CIN) is associated with (a) **reduced** T-cell infiltration but (b) **increased LAG3
(not TIGIT) per infiltrating T cell**, with elevated per-cell IFN-γ. Derived from TCGA bulk RNA
with crude proxies; needs (1) proxy-hardening, (2) composition-free confirmation, (3) independent
replication, (4) protein/clinical relevance.

## Key structural conclusion from the search
**No single open dataset has LAG3 protein + T-cell localization + genomics-for-CIN in the same
NSCLC samples.** So the validation strategy is layered: harden in TCGA → confirm per-cell in
single-cell → replicate in open bulk cohorts → corroborate at protein/clinical level (some access-gated).
Notably, the closest protein evidence (**Datar & Schalper 2019, CCR**) independently found that
**baseline LAG-3 in TILs predicts anti–PD-1 non-response while PD-1/TIM-3 do not** — a protein-level
echo of our LAG3-specific, per-T-cell hypothesis.

---

## Tier A — Harden the TCGA analysis (open, one-click GDC downloads)
Replace our two proxies and add the critical purity confounder.

| Resource | What it gives | Download | Replaces / adds |
|---|---|---|---|
| **Thorsson 2018 leukocyte fraction** | Methylation-based immune infiltrate per TCGA sample (RNA-independent) | GDC panimmune `TCGA_all_leuk_estimate.masked.20170107.tsv` | our self-computed CD8 score |
| **Thorsson CIBERSORT fractions** | 22 immune subsets incl. CD8 (× leukocyte fraction = absolute CD8) | GDC panimmune `TCGA.Kallisto...cibersort.relative.tsv` | independent CD8 abundance |
| **Taylor 2018 Aneuploidy Score** | Purity-aware per-sample AS (0–39) + arm-level calls | GDC pancan-aneuploidy `PANCAN_ArmCallsAndAneuploidyScore_092817.txt` | our FGA proxy |
| **ABSOLUTE purity/ploidy/WGD** | Tumor purity, ploidy, genome-doublings | GDC pancanatlas `TCGA_mastercalls.abs_tables_JSedit.fixed.txt` | **purity confounder** (purity dilutes bulk immune signal) + WGD as 2nd CIN axis |
| **TIMER2.0 / immunedeconv** | 6-method CD8 estimates | timer.cistrome.org / omnideconv/immunedeconv | multi-method robustness |
| **CIBERSORTx high-res / CODEFACS** | **cell-type-specific LAG3 from bulk** | cibersortx.stanford.edu (token) / ruppinlab | turns "per-cell LAG3" from inferred → measured |
| **Saltz 2018 TIL maps** | H&E image-based TIL fraction (orthogonal) | GDC tilmap / TCIA | image confirmation of infiltration ↓ |

## Tier B — Single-cell LUAD (open; composition-free per-cell test via inferCNV/CopyKAT)

| Dataset | n | Access | Note |
|---|---|---|---|
| **Salcher 2022 LuCA** (best-powered) | 1.3M cells / 318 pts, LUAD-dominant | OPEN — CELLxGENE (h5ad/rds) + Zenodo | only atlas where driver-negative is a viable subgroup; genotype for a subset; batch/integration caveat |
| **Kim 2020 GSE131907** | 208k cells / 44 pts | OPEN — GEO | largest single consistent 10x LUAD; malignant + deep T; driver annotation partial |
| **Bischoff 2021** | 114k cells / 10 LUAD | OPEN — Code Ocean | most complete per specimen (deep CD8 + malignant + genotype); small n |
| **Maynard 2020 PRJNA591860** | 23k cells / 30 pts, Smart-seq2 | OPEN (processed, Code Ocean) | truest per-cell reads + clear driver status; treatment-exposed |
| **Laughney 2020 GSE123902 / Lambrechts E-MTAB-6149/6653** | small | OPEN | supporting replication |

## Tier C — Independent bulk replication (CNA+RNA+mutations, same tumors)

| Cohort | n LUAD | Access | CIN? | Driver-neg? |
|---|---|---|---|---|
| **OncoSG** `luad_oncosg_2020` | 305 (East-Asian) | **OPEN — cBioPortal** | ✅ GISTIC | ✅ |
| **CPTAC-LUAD** `luad_cptac_2020` | 110 (+ NAT) | **OPEN — cBioPortal**; +**proteomics** | ✅ | ✅ |
| **TRACERx421** | ~248 LUAD, multi-region | CIN/driver OPEN (suppl.); **RNA EGA-controlled** | ✅ (open tables: wGII/WGD) | ✅ | 
| **Korean EGFR/ALK-WT never-smoker** | 99, oncogene-neg *by design* | Controlled (KoNA); proteomics open (PXD033360) | ✅ | ✅ exact-match |
| **Sherlock-Lung** | 232 never-smoker | RNA open GSE171415; WGS dbGaP | ✅ (controlled) | ✅ |

## Tier D — Immunotherapy response & LAG3 protein/spatial

| Dataset | What | Access |
|---|---|---|
| **Ravi 2023 / SU2C-MARK** | 393 NSCLC ICI; WES+RNA+response+survival — supports the *entire* hypothesis | **Controlled — dbGaP phs002822** |
| **Rizvi 2015** `luad_mskcc_2015` / **Hellmann 2018** `nsclc_mskcc_2018` | WES + ICI response (CIN/TMB derivable; no RNA) | **OPEN — cBioPortal** |
| **GSE126044 / GSE135222** | anti-PD-1 NSCLC bulk RNA + response (no CNA) | **OPEN — GEO** |
| **Datar & Schalper 2019** | QIF LAG-3/PD-1/TIM-3 protein per TIL, >800 NSCLC + 90 ICI — *the* on-point protein biology | By request / MTA (Yale) |
| **Sorin 2023 Nature IMC** | 416 LUAD, 35-plex, 1.6M cells + survival — T-cell/tumor spatial (no LAG3 channel, no CIN) | **OPEN — Zenodo 7383627** |
| **Walsh/Quail 2025 IMC** | 102 NSCLC spatial (no LAG3, no CIN) | OPEN — Zenodo 14625562 |

**Clinical fact base:** RELATIVITY-104 (nivolumab+relatlimab 1L NSCLC; ORR 51.3% vs 43.7%, benefit
strongest in PD-L1≥1% non-squamous) and TACTI-002 (eftilagimod+pembrolizumab) establish LAG-3 as an
actionable co-checkpoint in NSCLC — raising the stakes of a CIN→LAG3 axis.
