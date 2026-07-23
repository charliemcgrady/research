# Single-cell LAG3 validation pipeline (cohort-agnostic)

Tests whether **tumor** chromosomal instability (CIN) associates with **CD8 T-cell LAG3
(vs TIGIT)** and whether LAG3-high CD8 cells are an activated vs exhausted state — ideally
in **driver-negative LUAD**. Built for a resourced environment (LuCA-scale atlases).

## Design guarantees (why this is statistically honest)
- **Independent unit = tumor, not cell.** Primary inference is tumor/patient-level pseudobulk
  (`tumor_level_tests`); a patient-clustered mixed model (`cell_level_mixedmodel`) is only a
  sensitivity check and its note states inference is still bounded by n_patients.
- **CIN is a malignant-cell property** inferred per malignant cell (infercnvpy) and aggregated
  to one value per tumor.
- **Cohort/platform adjusted** (covariate; per-dataset meta-analysis advised for >1 platform).
- **Driver-negative is gated:** without a populated driver column the script refuses to present
  results as validating the driver-negative claim and says so in `results.json` (`DRIVER_CLAIM`).
- **Audit-first:** reports the number of eligible independent tumors before any inferCNV; stops
  if <5.

## Files
- `sc_lag3_pipeline.py` — main pipeline (audit → per-tumor CIN → CD8 pseudobulk → tumor-level
  tests + mixed-model sensitivity → CD8-state description). Cohort-agnostic via `--col-*` args.
- `kim_adapter.py` — builds a Kim GSE131907 h5ad (smoke/descriptive cohort; ~13 tumors → exploratory).
- `smoke_test.py` — synthetic multi-tumor data with a *planted* CIN→LAG3 signal; verifies the
  pipeline recovers it (LAG3 coef>0 p<0.05, TIGIT null, n=tumors). **PASSES** (coef≈0.62, p≈1e-21).

## Prereqs
`pip install scanpy infercnvpy anndata statsmodels scipy` and a GENCODE GTF
(e.g. `gencode.v44.basic.annotation.gtf.gz`) for gene positions.

## Recommended cohorts (see ../EVIDENCE_CATALOG.md)
- **LuCA** (Salcher 2022, CELLxGENE) — properly powered (hundreds of patients) *if* driver
  annotation is adequately populated (run the feasibility audit first).
- **Kim GSE131907** — smoke test + descriptive CD8-state only (~13 tumors; not definitive).

## Run
```bash
# 1. audit only (cheap — decide feasibility before inferCNV)
python sc_lag3_pipeline.py --h5ad luca.h5ad --audit-only \
  --col-sample sample --col-patient patient --col-celltype cell_type_major \
  --malignant "Epithelial cell (malignant)" --cd8 "T cell CD8" \
  --col-histology condition --luad-values "LUAD" \
  --col-treatment treatment --naive-values "naive" \
  --col-driver driver_mutation --col-platform platform

# 2. full run (adds --gtf, --reference normal cell types, --driver-negative-values)
```

## Key caveat carried from bulk (E1b)
Bulk showed the per-cell LAG3 signal is purity-robust and LAG3-specific but **covariate-dependent**
(significant with an expression-based T-cell score, weaker with independent CIBERSORT-abs). Single
cell removes that deconvolution ambiguity — but only if an adequately powered, driver-annotated
cohort exists. If not, single-cell validation of the *driver-negative* claim is **infeasible**, and
this pipeline will say so rather than over-generalize.
