# Feasibility Gate — Decision (mission step 6)

**Decision: NARROW SCOPE, proceed on open data.** (Autonomous call after the access report +
harmonized driver census; the gate question was left to the analyst.)

## Basis
- **Oncogene-negative n is small but real:** TCGA-LUAD DN-Strict **~48** (≈57 true pan-negative per
  Chen 2021); CPTAC-LUAD **~18**; DepMap **~8–15** lines. Enough for *exploratory / hypothesis-generating*
  multi-omic and functional work with honest power caveats; **not** enough for well-powered
  lineage×alteration interaction tests with 16-covariate control.
- **The mission's primary novel targets are blocked and pre-empted.** Structural rearrangement /
  enhancer-hijacking / promoter / cryptic-splice / amplicon-geometry (Phase 3, Exp 11–20) need
  genome-wide WGS SVs, which are **controlled (dbGaP/EGA)** — no open cohort has ≥100 LUAD WGS+SV+RNA —
  and were **already characterized in the same TCGA set by Chen et al. 2021 Cell Reports**
  (STK11/KEAP1 promoter deletions, ILF2 promoter mutations, ecDNA/BFB amplicons).

## Scope IN (feasible on open data; will execute)
- **Phase 1** cohort + lineage validation (Exp 1–5).
- **Phase 2** expression/protein/**phospho-state outliers** (Exp 6–10) — the most novel *open* angle
  (CPTAC protein/phospho, n≈18) + TCGA RNA outliers (n≈48–57).
- **Phase 4** lineage-conditioning of *expression/protein* classes (Exp 21–23, 25–30; not 24).
- **Phase 5** pathway-output validation (Exp 31–35).
- **Phase 6** functional dependency in DepMap (Exp 36–40).
- **Phase 7** lineage-state entropy (Exp 41–45).
- **Phase 8** ecosystem/spatial where open (Exp 46–48, partial).
- **Phase 9** integrated classifier + TCGA↔CPTAC replication (Exp 49–50), power-limited.

## Scope OUT — Failure ledger (blocked / not testable on open data)
| Experiment(s) | Status | Reason |
|---|---|---|
| 11, 12, 13, 14 (local SVs, enhancer proximity, super-enhancer overlap) | **not testable (open)** | genome-wide WGS SV controlled; no open ≥100-LUAD WGS+SV cohort |
| 17 (cryptic splice, genome-wide) | **degraded** | GDC junctions not open; recount3/Snaptron only → different pipeline, QC-limited |
| 18 (noncoding promoter mutations) | **not testable (open)** | needs WGS non-coding calls (controlled) |
| 19 (focal amplification geometry / ecDNA) | **not testable (open)** + pre-empted | WGS reconstruction controlled; Chen 2021 already reported ecDNA/BFB |
| 20, 24 (composite cis-regulatory; lineage×enhancer class) | **blocked** | depend on 11–19 |
| 15, 16 (allele-specific expression, alt promoter) | **partial** | ASE needs phased het SNPs (RNA-only approximation weak); alt-promoter via recount3 only |
| Mucinous / gastric-like (HNF4A) lineage arm | **unvalidatable (open)** | no open human mucinous-LUAD scRNA reference |

## Standing caveats applied to every downstream result
Small n (report exact n + CI + FDR family), tumor purity, smoking, TP53/STK11/KEAP1/NF1/PTEN,
proliferation, WGD, CNA burden, TMB, stage, cohort/platform; negative controls; TCGA=discovery,
CPTAC=validation (no threshold tuning on CPTAC); regulatory tracks are contextual, not tumor-specific.

## Immediate next steps (narrowed order)
1. Build + validate the lineage-state model (LuCA malignant cells + normal-lung refs → project to
   TCGA & CPTAC bulk; exclude candidate genes to prevent circularity).
2. Run Phase 1 Exp 1–5 (census already preliminary; formalize + false-negative + definition-stability
   + lineage reproducibility + lineage positive controls).
3. Proceed through the in-scope phases, logging every null/blocked result honestly.
