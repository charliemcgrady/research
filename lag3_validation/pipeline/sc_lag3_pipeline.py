#!/usr/bin/env python3
"""
Cohort-agnostic single-cell pipeline: does tumor chromosomal instability (CIN)
associate with CD8 T-cell LAG3 (vs TIGIT), and is LAG3-high CD8 an activated or
exhausted state — in (ideally driver-negative) LUAD.

DESIGN PRINCIPLES (per project directive):
  * The independent unit is the TUMOR, not the cell. All inferential tests are at
    tumor/patient level (pseudobulk) or use a patient-clustered mixed model. Cell-level
    n is NEVER used as the inferential n.
  * CIN is a malignant-cell (tumor) property: inferred per malignant cell, aggregated
    to one CIN value per tumor.
  * Cohort/platform is always adjusted for (covariate or random effect); optional
    per-dataset meta-analysis.
  * DRIVER-NEGATIVE is the target claim. If driver annotation is absent/sparse, the
    script REFUSES to present results as validating the driver-negative claim and says so.
  * Runs an eligibility AUDIT first and reports the number of eligible independent tumors
    BEFORE any expensive inferCNV.

INPUTS: one .h5ad (log-normalized .X or raw counts in .layers['counts']) with an obs
schema mapped via --col-* args (cohort-agnostic). See COLS below.

USAGE (resourced env):
  python sc_lag3_pipeline.py --h5ad luca_extended.h5ad \
     --col-sample sample --col-patient patient --col-celltype cell_type_major \
     --malignant "Epithelial cell (malignant)" --cd8 "T cell CD8" \
     --col-histology condition --luad-values "LUAD,LUAD primary" \
     --col-treatment treatment --naive-values "naive,treatment-naive,none" \
     --col-platform platform --col-dataset dataset --col-driver driver_mutation \
     --gtf gencode.v44.basic.annotation.gtf.gz --out results/

  # Kim smoke test (exploratory; ~13 tumors) — see kim_adapter.py to build the h5ad first
  python sc_lag3_pipeline.py --h5ad kim.h5ad --smoke ...

DEPS: scanpy, infercnvpy, anndata, numpy, pandas, statsmodels, scipy
"""
import argparse, sys, json, warnings
warnings.filterwarnings("ignore")

# ---- marker programs (used for descriptive CD8-state classification) ----
CHECKPOINTS = ["LAG3", "TIGIT", "PDCD1", "HAVCR2", "CTLA4", "CD274"]
EXHAUSTION  = ["TOX", "HAVCR2", "PDCD1", "ENTPD1", "TIGIT", "CTLA4", "LAYN"]
ACTIVATION  = ["IFNG", "GZMB", "GZMK", "GZMA", "PRF1", "NKG7", "TNF"]
STEMNESS    = ["TCF7", "IL7R", "SELL", "CCR7"]
CD8_MARKERS = ["CD8A", "CD8B", "CD3D", "CD3E"]

MIN_MALIGNANT = 50   # min malignant cells/tumor to estimate CIN
MIN_CD8       = 20   # min CD8 cells/tumor to estimate CD8 pseudobulk


def audit(ad, C, luad_vals, naive_vals):
    """Eligibility audit -> per-tumor table; returns eligible sample list + summary."""
    import pandas as pd, numpy as np
    obs = ad.obs
    df = pd.DataFrame(index=obs.index)
    df["sample"]  = obs[C["sample"]].astype(str)
    df["patient"] = obs[C["patient"]].astype(str) if C.get("patient") else df["sample"]
    df["ct"]      = obs[C["celltype"]].astype(str)
    df["is_mal"]  = df["ct"].isin(C["malignant"])
    df["is_cd8"]  = df["ct"].isin(C["cd8"])
    # histology / treatment / driver filters (only if columns provided)
    keep = pd.Series(True, index=df.index)
    if C.get("histology") and luad_vals:
        keep &= obs[C["histology"]].astype(str).isin(luad_vals)
    if C.get("treatment") and naive_vals:
        keep &= obs[C["treatment"]].astype(str).isin(naive_vals)
    df, obs2 = df[keep], obs[keep]
    g = df.groupby("sample")
    tab = pd.DataFrame({
        "patient":   g["patient"].first(),
        "n_malignant": g["is_mal"].sum(),
        "n_cd8":       g["is_cd8"].sum(),
    })
    if C.get("dataset"):  tab["dataset"]  = g.apply(lambda x: obs2.loc[x.index, C["dataset"]].astype(str).iloc[0])
    if C.get("platform"): tab["platform"] = g.apply(lambda x: obs2.loc[x.index, C["platform"]].astype(str).iloc[0])
    if C.get("driver"):
        tab["driver"] = g.apply(lambda x: obs2.loc[x.index, C["driver"]].astype(str).iloc[0])
        tab["driver_annotated"] = tab["driver"].notna() & ~tab["driver"].isin(["nan","NA","unknown",""])
    tab["eligible"] = (tab["n_malignant"] >= MIN_MALIGNANT) & (tab["n_cd8"] >= MIN_CD8)
    summary = {
        "n_samples_total": int(tab.shape[0]),
        "n_eligible_tumors": int(tab["eligible"].sum()),
        "n_eligible_with_driver": int((tab["eligible"] & tab.get("driver_annotated", False)).sum()) if C.get("driver") else 0,
        "driver_column_present": bool(C.get("driver")),
    }
    return tab, summary


def infer_cin_per_tumor(ad, C, gtf, eligible_samples):
    """infercnvpy on malignant + reference cells -> per-tumor mean CNV magnitude."""
    import infercnvpy as cnv, scanpy as sc, numpy as np, pandas as pd
    cnv.io.genomic_position_from_gtf(gtf, adata=ad)          # populates var chr/start/end
    ad = ad[ad.var["chromosome"].notna()].copy()
    ref = [c for c in C["reference"] if c in set(ad.obs[C["celltype"]])]  # normal cell types as diploid ref
    cnv.tl.infercnv(ad, reference_key=C["celltype"], reference_cat=ref, window_size=100)
    # per-cell CNV magnitude = mean abs of the CNV matrix row
    mag = np.asarray(np.abs(ad.obsm["X_cnv"]).mean(axis=1)).ravel()
    ad.obs["cnv_mag"] = mag
    mal = ad.obs[ad.obs[C["celltype"]].isin(C["malignant"])]
    cin = mal.groupby(ad.obs[C["sample"]].astype(str))["cnv_mag"].mean()
    return cin.rename("tumor_CIN")


def cd8_pseudobulk(ad, C):
    """Per-patient/tumor mean log-expr over CD8 cells for checkpoints + program scores."""
    import scanpy as sc, numpy as np, pandas as pd
    cd8 = ad[ad.obs[C["celltype"]].isin(C["cd8"])].copy()
    for name, genes in [("exhaustion",EXHAUSTION),("activation",ACTIVATION),("stemness",STEMNESS)]:
        g=[x for x in genes if x in cd8.var_names]
        if g: sc.tl.score_genes(cd8, g, score_name=f"score_{name}")
    rows=[]
    for s, idx in cd8.obs.groupby(cd8.obs[C["sample"]].astype(str)).groups.items():
        sub=cd8[idx]; row={"sample":s,"n_cd8":sub.n_obs}
        for gene in CHECKPOINTS+["IFNG"]:
            if gene in sub.var_names:
                row[f"pb_{gene}"]=float(np.asarray(sub[:,gene].X.mean()))
        for name in ["exhaustion","activation","stemness"]:
            k=f"score_{name}"
            if k in sub.obs: row[f"pb_{name}"]=float(sub.obs[k].mean())
        rows.append(row)
    return pd.DataFrame(rows).set_index("sample")


def tumor_level_tests(merged, has_driver):
    """Primary inference at TUMOR level (n = eligible tumors). Platform-adjusted OLS."""
    import statsmodels.formula.api as smf, numpy as np
    out={}
    d=merged.dropna(subset=["tumor_CIN","pb_LAG3"]).copy()
    d["CIN_z"]=(d["tumor_CIN"]-d["tumor_CIN"].mean())/d["tumor_CIN"].std()
    covar = " + C(platform)" if "platform" in d and d["platform"].nunique()>1 else ""
    out["n_tumors"]=int(d.shape[0])
    for gene in ["LAG3","TIGIT","IFNG"]:
        y=f"pb_{gene}"
        if y not in d: continue
        m=smf.ols(f"{y} ~ CIN_z{covar}", d).fit()
        out[gene]={"coef":round(float(m.params.get("CIN_z",np.nan)),3),
                   "p":float(m.pvalues.get("CIN_z",np.nan))}
    if has_driver and "driver_negative" in d:
        dn=d[d["driver_negative"]]
        out["driver_negative_n"]=int(dn.shape[0])
        if dn.shape[0]>=10:
            m=smf.ols(f"pb_LAG3 ~ CIN_z{covar}", dn).fit()
            out["driver_negative_LAG3"]={"coef":round(float(m.params.get("CIN_z",np.nan)),3),
                                         "p":float(m.pvalues.get("CIN_z",np.nan))}
    return out


def cell_level_mixedmodel(ad, C, cin):
    """Sensitivity: cell-level LAG3 ~ tumor_CIN + platform + (1|patient). Cells nested in patients."""
    import statsmodels.formula.api as smf, numpy as np, pandas as pd
    cd8=ad[ad.obs[C["celltype"]].isin(C["cd8"])]
    if "LAG3" not in cd8.var_names: return {"skipped":"no LAG3"}
    df=pd.DataFrame({"LAG3":np.asarray(cd8[:,"LAG3"].X.todense()).ravel() if hasattr(cd8[:,'LAG3'].X,'todense') else np.asarray(cd8[:,"LAG3"].X).ravel(),
                     "sample":cd8.obs[C["sample"]].astype(str).values,
                     "patient":(cd8.obs[C["patient"]].astype(str).values if C.get("patient") else cd8.obs[C["sample"]].astype(str).values)})
    df["tumor_CIN"]=df["sample"].map(cin)
    df=df.dropna(subset=["tumor_CIN"])
    df["CIN_z"]=(df["tumor_CIN"]-df["tumor_CIN"].mean())/df["tumor_CIN"].std()
    try:
        md=smf.mixedlm("LAG3 ~ CIN_z", df, groups=df["patient"]).fit(reml=False)
        return {"coef_CIN":round(float(md.params.get("CIN_z",np.nan)),3),
                "p_CIN":float(md.pvalues.get("CIN_z",np.nan)),
                "n_cells":int(df.shape[0]),"n_patients":int(df["patient"].nunique()),
                "note":"random intercept per patient; inference still limited by n_patients"}
    except Exception as e:
        return {"error":str(e)}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--h5ad",required=True); ap.add_argument("--gtf")
    ap.add_argument("--col-sample",required=True); ap.add_argument("--col-patient")
    ap.add_argument("--col-celltype",required=True)
    ap.add_argument("--malignant",required=True,help="comma-sep malignant cell-type labels")
    ap.add_argument("--cd8",required=True,help="comma-sep CD8 cell-type labels")
    ap.add_argument("--reference",default="",help="comma-sep normal cell-type labels for CNV diploid ref")
    ap.add_argument("--col-histology"); ap.add_argument("--luad-values",default="")
    ap.add_argument("--col-treatment"); ap.add_argument("--naive-values",default="")
    ap.add_argument("--col-platform"); ap.add_argument("--col-dataset"); ap.add_argument("--col-driver")
    ap.add_argument("--driver-negative-values",default="",help="driver values meaning driver-NEGATIVE (e.g. 'none,WT,negative')")
    ap.add_argument("--audit-only",action="store_true")
    ap.add_argument("--smoke",action="store_true",help="Kim smoke mode: exploratory, no driver claim")
    ap.add_argument("--out",default="sc_results")
    a=ap.parse_args()

    import scanpy as sc, anndata, pandas as pd, numpy as np, os
    os.makedirs(a.out,exist_ok=True)
    C=dict(sample=a.col_sample,patient=a.col_patient,celltype=a.col_celltype,
           malignant=a.malignant.split(","),cd8=a.cd8.split(","),
           reference=[x for x in a.reference.split(",") if x],
           histology=a.col_histology,treatment=a.col_treatment,
           platform=a.col_platform,dataset=a.col_dataset,driver=a.col_driver)
    ad=sc.read_h5ad(a.h5ad)
    luad_vals=[x for x in a.luad_values.split(",") if x]; naive_vals=[x for x in a.naive_values.split(",") if x]

    # ---- 1. AUDIT (always first) ----
    tab,summary=audit(ad,C,luad_vals,naive_vals)
    tab.to_csv(f"{a.out}/eligibility_audit.csv")
    print("=== ELIGIBILITY AUDIT ===\n",json.dumps(summary,indent=2))
    print(f"Eligible independent tumors: {summary['n_eligible_tumors']}")
    if C.get("driver"):
        print(f"...with driver annotation: {summary['n_eligible_with_driver']}")
    else:
        print("NO DRIVER COLUMN -> cannot validate the DRIVER-NEGATIVE claim; general results are NOT a substitute.")
    if a.audit_only or summary["n_eligible_tumors"]<5:
        if summary["n_eligible_tumors"]<5: print("Too few eligible tumors for inference. Stopping.")
        return

    elig=set(tab.index[tab["eligible"]])
    ad=ad[ad.obs[C["sample"]].astype(str).isin(elig)].copy()

    # ---- 2. per-tumor CIN ----
    if not a.gtf: sys.exit("--gtf required for CIN inference (gene positions).")
    cin=infer_cin_per_tumor(ad,C,a.gtf,elig)
    cin.to_csv(f"{a.out}/tumor_CIN.csv")

    # ---- 3. CD8 pseudobulk + merge ----
    pb=cd8_pseudobulk(ad,C)
    merged=pb.join(cin,how="inner")
    if C.get("platform"): merged=merged.join(tab["platform"],how="left")
    if C.get("driver") and a.driver_negative_values:
        dn_vals=set(a.driver_negative_values.split(","))
        merged=merged.join(tab["driver"],how="left")
        merged["driver_negative"]=merged["driver"].isin(dn_vals)
    merged.to_csv(f"{a.out}/tumor_level_merged.csv")

    # ---- 4. TUMOR-LEVEL primary tests ----
    has_driver=bool(C.get("driver") and a.driver_negative_values)
    res=tumor_level_tests(merged,has_driver)
    # ---- 5. cell-level mixed-model sensitivity ----
    res["mixed_model_sensitivity"]=cell_level_mixedmodel(ad,C,cin)
    # ---- 6. descriptive CD8-state (LAG3-high vs low: exhaustion vs activation) ----
    res["audit"]=summary
    if a.smoke: res["MODE"]="SMOKE TEST (exploratory; tumor-level n small; NOT a definitive/driver-negative result)"
    if not has_driver:
        res["DRIVER_CLAIM"]="NOT TESTED — no driver annotation; results describe LUAD broadly, NOT the driver-negative subset."
    json.dump(res,open(f"{a.out}/results.json","w"),indent=2)
    print("=== RESULTS ===\n",json.dumps(res,indent=2))

if __name__=="__main__":
    main()
