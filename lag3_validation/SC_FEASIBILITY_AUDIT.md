# Single-cell feasibility audit — can we validate the DRIVER-NEGATIVE per-cell LAG3 claim?

**Verdict: No — not at defensible power with current public single-cell data.** LuCA is the
only atlas at the right scale, and its driver annotation is too sparse and too soft to support a
per-tumor, single-cell test of a *driver-negative*-specific claim. (Counts parsed directly from the
atlas's own `patient_metadata_corrected.xlsx`, 303 patients.)

## Why (the gate is driver annotation)
- LuCA **does** carry a driver field (`driver_genes` + per-gene EGFR/KRAS/TP53/ALK/BRAF/ERBB2/ROS
  columns) and **does** separate "malignant cell" from "CD8+ T cell" — so the cell types are there.
- But driver annotation is **sparse**: **75 / 161 LUAD patients**; and the large *treatment-naive
  primary* LUAD cohorts (Kim 44, Wu 18, Guo 11, …) have **no** driver annotation at all. Annotated
  LUAD collapse to ~2 usable cohorts (Maier 24 naive-early; Laughney 14 primary).
- **Explicitly oncogene-negative LUAD ≈ 5–9 tumors**, from 2 cohorts — and "negative" is *inferred
  from a list*, not a documented wild-type test (absence ≠ tested-negative). Single digits, soft label.
- **No CNV in the atlas** — malignant cells are labeled by expression clustering, so we'd run
  inferCNV ourselves on expression-defined malignant cells (mild circularity) with no per-cell CN ground truth.
- **The LuCA authors themselves** ran their genotype↔CD8 associations on **TCGA bulk (1,026 pts) via
  Scissor**, not on their single-cell patients — an explicit admission that single-cell driver
  annotation is too thin for this exact question.

## Eligible independent LUAD tumors (tumor = the unit)
| Filter | Estimated N |
|---|---|
| LUAD + primary + naive + both malignant & CD8 (ignoring driver) | **~70–120** (upper band; exact needs an h5ad `sample × cell_type` crosstab) |
| …AND driver-annotated | ~35–40 |
| …AND **driver-negative** specifically | **~5–9** ← underpowered |

## What this means for the plan
1. **The driver-negative-specific per-cell claim cannot be validated in single cell right now.** Per
   directive, we state this plainly rather than pass off general-LUAD single-cell results as if they
   settled the driver-negative question.
2. **A single-cell test in LUAD *broadly* IS feasible** (~70–120 tumors) and would resolve the one
   real weakness of the bulk result — the deconvolution/covariate ambiguity from E1b (per-cell LAG3
   was purity-robust and LAG3-specific but sensitive to how infiltration was modeled). It just cannot
   be labeled driver-negative-specific.
3. Where the driver-negative claim therefore stands: on **bulk TCGA (n≈115 driver-negative)** with a
   per-cell inference that is purity-robust, LAG3-specific, but covariate-dependent — and **no public
   single-cell cohort can currently adjudicate it for the driver-negative subset.**

## Recommended next steps (ranked)
1. **Cheap, decisive first move:** download only the **LuCA core atlas** (or a cell-type subset) and
   run `sc_lag3_pipeline.py --audit-only` to get the *exact* per-tumor malignant+CD8 counts and the
   real driver-negative N — confirming the ~5–9 estimate before any compute.
2. **If you want the mechanism nailed (general LUAD):** run the full pipeline on LuCA ignoring driver
   status (~70–120 tumors) to settle per-cell vs composition for CIN→CD8-LAG3 broadly. Report it as
   all-LUAD, explicitly not driver-negative.
3. **For the driver-negative claim specifically:** it stays a bulk finding with the E1b caveat.
   Closing it would require either a driver-annotated scRNA LUAD cohort (not public at scale) or
   orthogonal validation (multiplex-IHC LAG3 protein per T cell in driver-stratified tumors — e.g. a
   Datar/Schalper-style TMA, access-gated).
4. Kim remains a **smoke-test/descriptive** cohort only (~13 tumors); the pipeline + smoke test are
   ready (`pipeline/`), and the Kim full-matrix download is capped in this container.
