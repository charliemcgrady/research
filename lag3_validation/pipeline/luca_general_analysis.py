#!/usr/bin/env python3
"""
General all-LUAD LuCA analysis (memory-frugal, checkpointed, restart-safe).
Tests tumor-level CIN -> CD8 LAG3 (vs TIGIT comparator) across eligible primary LUAD tumors.
NOT driver-negative-specific (n=1 there); framed as a GENERAL LUAD mechanism.

Stages (each writes a checkpoint; rerun skips finished stages):
  1 extract  -> OUT/subset.h5ad         (malignant + CD8 + normal-epithelial reference cells)
  2 cin      -> OUT/tumor_CIN.csv        (infercnvpy per-cell CNV -> per-tumor mean over malignant)
  3 analyze  -> OUT/results.json + tumor_level.csv
"""
import os, json, sys, numpy as np, pandas as pd, h5py
from scipy import sparse

H5="/home/user/research/lag3_validation/scrna/luca/luca_core.h5ad"
AUDIT="/home/user/research/lag3_validation/scrna/luca/luca_luad_audit.csv"
GTF="/home/user/research/lag3_validation/scrna/luca/gencode.basic.gtf.gz"
OUT="/home/user/research/lag3_validation/scrna/luca/general_out"; os.makedirs(OUT,exist_ok=True)

MALIGNANT=["Tumor cells"]; CD8=["T cell CD8"]
REF_TYPES=["Alveolar cell type 1","Alveolar cell type 2","Club","Ciliated","transitional club/AT2"]
REF_CAP=15000   # cap reference cells (fixed stride) to control memory
CHECKPOINTS=["LAG3","TIGIT","PDCD1","HAVCR2","CTLA4","IFNG"]
EXH=["TOX","HAVCR2","PDCD1","ENTPD1","TIGIT","CTLA4","LAYN"]
ACT=["IFNG","GZMB","GZMK","GZMA","PRF1","NKG7"]
ENSG={"LAG3":"ENSG00000089692","TIGIT":"ENSG00000181847","IFNG":"ENSG00000111537",
 "PDCD1":"ENSG00000188389","HAVCR2":"ENSG00000135077","CTLA4":"ENSG00000163599",
 "TOX":"ENSG00000198846","ENTPD1":"ENSG00000138185","LAYN":"ENSG00000079215",
 "GZMB":"ENSG00000100453","GZMK":"ENSG00000113088","GZMA":"ENSG00000145649",
 "PRF1":"ENSG00000180644","NKG7":"ENSG00000105374","CD8A":"ENSG00000153563"}

def catcol(o,name):
    g=o[name]; cats=np.array([c.decode() if isinstance(c,bytes) else str(c) for c in g["categories"][:]],dtype=object)
    codes=g["codes"][:]; out=np.where(codes<0,None,cats[codes]); return out

def stage1_extract():
    if os.path.exists(f"{OUT}/subset.h5ad"): print("stage1: cached"); return
    import anndata as ad
    elig=set(pd.read_csv(AUDIT,index_col=0).query("eligible").index.astype(str))
    f=h5py.File(H5,"r"); o=f["obs"]
    sample=catcol(o,"sample").astype(str); ct=catcol(o,"cell_type_major").astype(str); origin=catcol(o,"origin").astype(str)
    in_elig=np.isin(sample,list(elig))
    is_mal=in_elig&np.isin(ct,MALIGNANT); is_cd8=in_elig&np.isin(ct,CD8)
    is_ref=np.isin(ct,REF_TYPES)&np.isin(origin,["normal","normal_adjacent"])
    ref_idx=np.where(is_ref)[0]
    if len(ref_idx)>REF_CAP: ref_idx=ref_idx[::max(1,len(ref_idx)//REF_CAP)][:REF_CAP]
    sel=np.sort(np.concatenate([np.where(is_mal)[0],np.where(is_cd8)[0],ref_idx]))
    print(f"stage1: malignant={int(is_mal.sum())} cd8={int(is_cd8.sum())} ref={len(ref_idx)} total={len(sel)}")
    # CSR row extraction
    Xg=f["X"]; indptr=Xg["indptr"][:]; n_genes=Xg.attrs["shape"][1]
    data_ds=Xg["data"]; ind_ds=Xg["indices"]
    rows_data=[]; rows_ind=[]; new_indptr=[0]
    for i in sel:
        s,e=indptr[i],indptr[i+1]
        rows_data.append(data_ds[s:e]); rows_ind.append(ind_ds[s:e]); new_indptr.append(new_indptr[-1]+(e-s))
    X=sparse.csr_matrix((np.concatenate(rows_data),np.concatenate(rows_ind),np.array(new_indptr)),
                        shape=(len(sel),n_genes))
    var=np.array([x.decode() if isinstance(x,bytes) else x for x in f["var"]["_index"][:]])
    var=np.array([v.split('.')[0] for v in var])
    obs=pd.DataFrame({"sample":sample[sel],"cell_type":ct[sel],"origin":origin[sel],
                      "dataset":catcol(o,"dataset").astype(str)[sel],
                      "platform":catcol(o,"platform").astype(str)[sel]})
    A=ad.AnnData(X,obs=obs,var=pd.DataFrame(index=var))
    A.write_h5ad(f"{OUT}/subset.h5ad"); print("stage1: wrote subset.h5ad",A.shape)

def stage2_cin():
    if os.path.exists(f"{OUT}/tumor_CIN.csv"): print("stage2: cached"); return
    import anndata as ad, infercnvpy as cnv, scanpy as sc, gzip
    A=ad.read_h5ad(f"{OUT}/subset.h5ad")
    sc.pp.normalize_total(A,target_sum=1e4); sc.pp.log1p(A)   # ensure log-norm for CNV
    # build ENSG(version-stripped) -> (chr,start,end) from GTF directly (robust to id versions)
    pos={}
    with gzip.open(GTF,"rt") as fh:
        for line in fh:
            if line.startswith("#"): continue
            p=line.split("\t")
            if len(p)>8 and p[2]=="gene":
                gid=None
                for fld in p[8].split(";"):
                    fld=fld.strip()
                    if fld.startswith("gene_id"): gid=fld.split('"')[1].split(".")[0]; break
                if gid: pos[gid]=(p[0].replace("chr",""),int(p[3]),int(p[4]))
    A.var["chromosome"]=[("chr"+pos[g][0]) if g in pos else None for g in A.var_names]
    A.var["start"]=[pos[g][1] if g in pos else None for g in A.var_names]
    A.var["end"]=[pos[g][2] if g in pos else None for g in A.var_names]
    A=A[:,A.var["chromosome"].notna()].copy()
    A.var["start"]=A.var["start"].astype(int); A.var["end"]=A.var["end"].astype(int)
    A.obs["cnv_ref"]=np.where(A.obs["cell_type"].isin(REF_TYPES),"reference","tumor_or_T")
    cnv.tl.infercnv(A,reference_key="cnv_ref",reference_cat=["reference"],window_size=100)
    A.obs["cnv_mag"]=np.asarray(np.abs(A.obsm["X_cnv"]).mean(axis=1)).ravel()
    mal=A.obs[A.obs["cell_type"].isin(MALIGNANT)]
    cin=mal.groupby("sample")["cnv_mag"].agg(["mean","count"]).rename(columns={"mean":"tumor_CIN","count":"n_malignant"})
    ds=A.obs.groupby("sample")[["dataset","platform"]].first()
    cin.join(ds).to_csv(f"{OUT}/tumor_CIN.csv"); print("stage2: wrote tumor_CIN.csv",cin.shape)

def stage3_analyze():
    import anndata as ad, scanpy as sc, statsmodels.formula.api as smf
    A=ad.read_h5ad(f"{OUT}/subset.h5ad")
    sc.pp.normalize_total(A,target_sum=1e4); sc.pp.log1p(A)
    cd8=A[A.obs["cell_type"].isin(CD8)].copy()
    def score(genes,name):
        ids=[ENSG[g] for g in genes if g in ENSG and ENSG[g] in cd8.var_names]
        if ids: sc.tl.score_genes(cd8,ids,score_name=name)
    score(EXH,"exhaustion"); score(ACT,"activation")
    rows=[]
    for s,idx in cd8.obs.groupby("sample").groups.items():
        sub=cd8[idx]; r={"sample":s,"n_cd8":sub.n_obs}
        for g in CHECKPOINTS:
            e=ENSG.get(g)
            if e in sub.var_names: r[f"pb_{g}"]=float(np.asarray(sub[:,e].X.mean()))
        for nm in ["exhaustion","activation"]:
            if nm in sub.obs: r[f"pb_{nm}"]=float(sub.obs[nm].mean())
        rows.append(r)
    pb=pd.DataFrame(rows).set_index("sample")
    cin=pd.read_csv(f"{OUT}/tumor_CIN.csv",index_col=0).dropna(subset=["tumor_CIN"])
    m_all=pb.join(cin,how="inner")
    # restrict to PRIMARY tumors (audit carries origin)
    aud=pd.read_csv(AUDIT,index_col=0)
    prim=set(aud.index[aud["origin"].astype(str).str.contains("primary",case=False,na=False)].astype(str))
    m=m_all[m_all.index.astype(str).isin(prim)].copy()
    m["CIN_z"]=(m.tumor_CIN-m.tumor_CIN.mean())/m.tumor_CIN.std()
    m.to_csv(f"{OUT}/tumor_level.csv")
    cov=" + C(dataset)" if m["dataset"].nunique()>1 else ""
    res={"n_tumors_primary":int(m.shape[0]),"n_tumors_all_eligible":int(m_all.shape[0]),
         "n_datasets":int(m["dataset"].nunique())}
    for g in ["LAG3","TIGIT","IFNG"]:
        y=f"pb_{g}"
        if y in m:
            mod=smf.ols(f"{y} ~ CIN_z{cov}",m).fit()
            res[g]={"coef":round(float(mod.params.get("CIN_z",np.nan)),3),"p":float(mod.pvalues.get("CIN_z",np.nan))}
            # no-cohort-adjust sensitivity
            mod2=smf.ols(f"{y} ~ CIN_z",m).fit()
            res[g]["p_no_cohort_adj"]=float(mod2.pvalues.get("CIN_z",np.nan))
    res["FRAMING"]="GENERAL all-LUAD mechanism (tumor-level, n tumors above). NOT driver-negative-specific."
    json.dump(res,open(f"{OUT}/results.json","w"),indent=2); print("stage3 RESULTS:",json.dumps(res,indent=2))

if __name__=="__main__":
    stage=sys.argv[1] if len(sys.argv)>1 else "all"
    if stage in ("1","all","extract"): stage1_extract()
    if stage in ("2","all","cin"): stage2_cin()
    if stage in ("3","all","analyze"): stage3_analyze()
