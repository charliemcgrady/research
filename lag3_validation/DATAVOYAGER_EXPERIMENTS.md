# Data Voyager experiment plan — understanding the LAG3 effect (10 experiments)

Each experiment = (dataset, research question) for the DataVoyager `analyze-data` agent, which runs
code on a local tabular file. "Prep" = what must be assembled into an analysis-ready table first;
"Runnable now" = data already local or one-click open.

Central claim under test: **in driver-negative LUAD, CIN → fewer T cells but higher LAG3 (not TIGIT)
& IFN-γ per T cell.** The plan moves: harden → confirm per-cell → replicate → mechanism → clinical.

---

## Tier 1 — Harden the TCGA finding (runnable now; kills the biggest confounds)

**E1. Gold-standard re-test.** *In driver-negative TCGA-LUAD, does per-T-cell LAG3 still rise with
aneuploidy after replacing both proxies and adjusting for tumor purity?*
- Data: our `master_sample_table` + Taylor **Aneuploidy Score**, Thorsson **leukocyte fraction** &
  CIBERSORT **absolute CD8**, ABSOLUTE **purity/ploidy/WGD**. Prep: join GDC tables on barcode.
- Discriminates: whether the effect survives purity adjustment (purity dilutes bulk immune signal —
  the untested confound most likely to explain it away) and better metrics. **Highest priority.**

**E2. Infiltration-method robustness.** *Is the per-cell LAG3~CIN signal stable across 6 deconvolution
methods (CIBERSORT-ABS, quanTIseq, xCell, MCP-counter, EPIC, TIMER)?*
- Data: TIMER2.0 per-sample CD8 estimates (all methods) joined to LAG3/AS. Prep: pull TIMER2.0 LUAD table.
- Discriminates: whether the earlier regression-vs-ratio discordance is a method artifact or robust.

**E3. Cell-type-specific LAG3 (measured, not inferred).** *Estimate LAG3 expression within the CD8/T-cell
compartment per sample (CIBERSORTx high-res / CODEFACS) and regress on aneuploidy in driver-negative LUAD.*
- Data: TCGA-LUAD bulk RNA → CIBERSORTx high-res GEPs. Prep: run CIBERSORTx (token) with a LUAD scRNA
  reference; produces per-sample T-cell-specific LAG3.
- Discriminates: converts "per-cell LAG3" from a covariate-adjusted inference into a direct measurement.

## Tier 2 — Single-cell confirmation (composition-free; the decisive layer)

**E4. Per-cell LAG3 vs tumor CNV in scRNA.** *Within CD8 T cells, is LAG3 higher in patients whose
malignant cells show high inferred aneuploidy?*
- Data: **Kim GSE131907** (and/or **LuCA**); inferCNV/CopyKAT on malignant cells → per-tumor CNV burden;
  pseudobulk/per-cell CD8 LAG3. Prep: scRNA download + inferCNV + build per-cell/per-tumor table.
- Discriminates: the definitive composition-free test — bulk cannot separate this.

**E5. Reproduce the dissociation at single-cell.** *Does CNV-high associate with (a) lower CD8 fraction
AND (b) higher per-CD8 LAG3 in the same cohort?*
- Data: same scRNA tables as E4.
- Discriminates: confirms the two-part dissociation (exclusion + per-cell up-regulation) directly.

**E6. Selectivity + activation-vs-exhaustion state.** *Is the effect LAG3-selective vs TIGIT/PD-1, and are
high-LAG3 CD8 cells activated effectors (IFNG⁺, GZMB⁺) or terminally exhausted (TOX⁺, TCF7⁻)?*
- Data: scRNA CD8 subset with checkpoint + exhaustion/effector programs vs per-tumor CNV.
- Discriminates: interprets the elevated per-cell IFN-γ — genuine activation vs early/terminal exhaustion;
  central to the biological meaning.

## Tier 3 — Independent bulk replication (open cBioPortal)

**E7. OncoSG replication.** *Does infiltration-adjusted LAG3(>TIGIT)~CIN replicate in OncoSG LUAD?*
- Data: `luad_oncosg_2020` CNA+RNA+mutations. Prep: cBioPortal download → build master-style table +
  derive AS + driver-negative + a CD8 signature.
- Discriminates: independent East-Asian replication.

**E8. CPTAC replication at RNA *and* protein.** *Does it hold in CPTAC-LUAD, including LAG3/CD8 at the
protein level?*
- Data: `luad_cptac_2020` CNA+RNA+**proteomics**. Prep: build table incl. protein LAG3/CD8.
- Discriminates: cross-platform + cross-ancestry, and RNA→protein consistency.

## Tier 4 — Mechanism & clinical relevance

**E9. Functional coherence (ligand + STK11).** *Does CIN track with LAG3's functional context — FGL1 and
MHC-II expression — and does STK11 loss (our run5 focal hit) modulate the LAG3 axis?*
- Data: TCGA + OncoSG (FGL1, HLA-DR/CIITA, CN_STK11, AS, LAG3).
- Discriminates: whether a CIN→LAG3 axis is functionally actionable (ligand present) vs epiphenomenal;
  connects to the STK11-cold result.

**E10. Clinical / ICI relevance.** *Is high T-cell-normalized LAG3 associated with anti-PD-1 non-response
in NSCLC — and does the CIN→LAG3 axis mark a PD-1-resistant, LAG3-co-blockade-candidate subset?*
- Data: **GSE126044 + GSE135222** (open RNA + response) for the LAG3~response test; **Rizvi 2015 /
  Hellmann 2018** (open WES) for CIN~response; **Ravi/SU2C-MARK** (dbGaP, if approved) for the full
  CIN×LAG3×response model. Prep: assemble per-patient tables.
- Discriminates: whether the finding has therapeutic traction (mirrors Datar/Schalper 2019 protein result).

---

## Sequencing & access summary
- **Run immediately (open + mostly local):** E1, E2, E7, E9, and the open half of E10 (GSE126044/135222,
  Rizvi/Hellmann). E3 open but needs a CIBERSORTx run.
- **Needs download + preprocessing (open):** E4, E5, E6 (scRNA + inferCNV), E8 (CPTAC protein).
- **Access-gated (parallel applications):** Ravi/SU2C-MARK dbGaP (E10 full model); TRACERx EGA (confirmatory);
  Datar/Schalper MTA (protein gold standard).
- **Recommended first batch:** E1 → E3 → E4 (harden → measure → confirm), then E7/E8 replication.
