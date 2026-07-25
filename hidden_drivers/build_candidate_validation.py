#!/usr/bin/env python3
"""
Class-3 candidate validation — MET, DDR2, EGFR (CPTAC-LUAD).
Runs checks 1-6 from the directive using data on hand (proteome, phosphoproteome,
RNA, CNV, mutations, clinical purity/TMB). Robust outlier scoring uses median/MAD
(not standard z) because the oncogene-negative subset is small. Downstream pathway
modules are built WITHOUT the initiating receptor's own sites.

Checks 7 (functional/DepMap) and 8 (literature) are handled separately
(DepMap not fetchable in-env -> required_next_test; literature via Asta).

Writes: candidate_cards/validation_<GENE>.json  +  console summary.
"""
import csv, json, os
import numpy as np, pandas as pd
RAW="cptac_luad/raw/"; CB="cptac_luad/"; OUT="hidden_drivers/candidate_cards"
os.makedirs(OUT,exist_ok=True)

def norm_id(s): return s.lstrip("X").replace(".","-")
def load_cct(fn):
    with open(RAW+fn) as f:
        r=csv.reader(f,delimiter="\t"); hdr=next(r); cols=[norm_id(x) for x in hdr[1:]]
        idx=[]; data=[]
        for row in r:
            idx.append(row[0]); data.append([np.nan if v in ("NA","","NaN") else float(v) for v in row[1:]])
    return pd.DataFrame(data,index=idx,columns=cols)

prot=load_cct("HS_CPTAC_LUAD_proteome_ratio_NArm_TUMOR.cct")
pho =load_cct("HS_CPTAC_LUAD_phosphoproteome_ratio_norm_NArm_TUMOR.cct")
rna =load_cct("HS_CPTAC_LUAD_rnaseq_uq_rpkm_log2_NArm_TUMOR.cct")
cnv =load_cct("HS_CPTAC_LUAD_cnv_gene_LR.cct")
mutm=load_cct("HS_CPTAC_LUAD_somatic_mutation_gene.cbt")
clin=json.load(open(CB+"cptac_clinical_all.json"))
purity={norm_id(k):float(v) for k,v in clin.get("TSNET PURITY",{}).items() if v not in ("NA","","[Not Available]")}
tmb={norm_id(k):float(v) for k,v in clin.get("TMB_NONSYNONYMOUS",{}).items() if v not in ("NA","",None)}
onco_mut=json.load(open(CB+"cptac_onco_mutations.json"))
supp_mut=json.load(open(CB+"cptac_suppressor_mutations.json"))
cna_burden=(cnv.abs()>0.3).sum(axis=0)/cnv.shape[0]

def rz_series(s):  # robust z across tumors (median/MAD)
    x=s.astype(float); med=np.nanmedian(x); mad=np.nanmedian(np.abs(x-med))
    sc=1.4826*mad if mad>0 else np.nanstd(x); sc=sc if sc>0 else np.nan
    return (x-med)/sc
def robust_z(df,key,samp):
    if key not in df.index or samp not in df.columns: return np.nan
    return float(rz_series(df.loc[key]).get(samp,np.nan))
def percentile(df,key,samp):
    if key not in df.index or samp not in df.columns: return np.nan
    row=df.loc[key].astype(float); v=row.get(samp)
    if pd.isna(v): return np.nan
    return round(100*float((row<v).sum())/row.notna().sum(),1)
def missing(df,key):
    if key not in df.index: return np.nan
    return round(float(df.loc[key].isna().mean()),3)

# lineage (for lineage-matched robust z) — reuse dominant lineage from CPTAC taxonomy table
Tc=pd.read_csv("hidden_drivers/taxonomy/tumor_evidence_CPTAC.csv").set_index("pid")
def lineage_matched_rz(df,key,samp,lin):
    peers=[p for p in Tc.index[Tc.dominant_lineage==lin] if p in df.columns]
    if key not in df.index or len(peers)<4: return np.nan
    sub=df.loc[key,[p for p in peers if p in df.columns]].astype(float)
    med=np.nanmedian(sub); mad=np.nanmedian(np.abs(sub-med)); sc=1.4826*mad if mad>0 else np.nanstd(sub)
    return float((df.loc[key].get(samp)-med)/sc) if sc and sc>0 else np.nan

# signatures for stromal-source check (robust-z mean of member RNA)
CAF=["FAP","COL1A1","COL1A2","COL3A1","PDGFRB","ACTA2","THY1","LUM","DCN","POSTN"]
MAC=["CD68","CD163","CSF1R","MRC1","MSR1","LYZ"]
COLLAGEN=["COL1A1","COL1A2","COL3A1","COL5A1","COL6A1","COL6A3"]
def sig_rz(genes,samp):
    vs=[robust_z(rna,g,samp) for g in genes if g in rna.index]
    vs=[v for v in vs if pd.notna(v)]; return float(np.mean(vs)) if vs else np.nan
def sig_corr(genes, actseries):  # corr of activation across tumors with signature
    sig=pd.DataFrame({g:rz_series(rna.loc[g]) for g in genes if g in rna.index}).mean(axis=1)
    j=pd.concat([actseries.rename("act"),sig.rename("sig")],axis=1).dropna()
    return round(float(j.act.corr(j.sig)),2) if len(j)>5 else np.nan

# downstream modules (activation-relevant phosphosites), receptor EXCLUDED
DOWNSTREAM={
 "EGFR":{"ERK":["MAPK1","MAPK3"],"MEK":["MAP2K1","MAP2K2"],"AKT":["AKT1","AKT2"],
         "adaptor":["SHC1","GAB1","GRB2"],"SRC":["SRC"],"PLCG1":["PLCG1"]},
 "MET":{"adaptor":["GAB1"],"SHP2":["PTPN11"],"AKT":["AKT1","AKT2","PIK3CA"],
        "ERK":["MAPK1","MAPK3"],"STAT3":["STAT3"],"SRC":["SRC"]},
 "DDR2":{"SRC_family":["SRC","FYN","LYN","YES1"],"SHP2":["PTPN11"],"ERK":["MAPK1","MAPK3"],
         "AKT":["AKT1","AKT2"],"focal_adhesion":["PTK2","PXN","BCAR1"]},
}
NEGREG=["PTPN1","PTPN2","PTPRJ","CBL","ERRFI1","LRIG1","SPRY1","SPRY2","SPRY4","DUSP1","DUSP4","DUSP6"]

CANDIDATES={"MET":"C3N-02422","DDR2":"C3N-02587","EGFR":"C3N-02588"}
LIG_OF={"MET":["HGF"],"EGFR":["EGF","TGFA","AREG","EREG","HBEGF","BTC"],"DDR2":["COL1A1","COL1A2","COL3A1"]}

def module_score(gene,samp):
    out={}
    for mod,members in DOWNSTREAM[gene].items():
        sites=[s for s in pho.index for m in members if s.split(":")[0]==m]
        sites=[s for s in sites if s.split(":")[0]!=gene]   # exclude receptor
        # prefer activation-loop / known-active sites: keep max robust-z across the module's sites
        vals=[robust_z(pho,s,samp) for s in sites]; vals=[v for v in vals if pd.notna(v)]
        out[mod]={"n_sites":len(sites),"max_rz":round(max(vals),2) if vals else None,
                  "mean_rz":round(float(np.mean(vals)),2) if vals else None}
    return out

report={}
for gene,samp in CANDIDATES.items():
    lin=Tc.dominant_lineage.get(samp,"")
    # --- Check 1: phosphosites ---
    sites=[s for s in pho.index if s.split(":")[0]==gene]
    site_rows=[]
    act_series=None
    for s in sites:
        rzv=robust_z(pho,s,samp)
        site_rows.append({"site":s,"robust_z":round(rzv,2) if pd.notna(rzv) else None,
                          "percentile":percentile(pho,s,samp),"missingness":missing(pho,s),
                          "lineage_matched_rz":round(lineage_matched_rz(pho,s,samp,lin),2) if pd.notna(lineage_matched_rz(pho,s,samp,lin)) else None})
    site_rows=sorted(site_rows,key=lambda r:-(r["robust_z"] or -9))
    top_site=site_rows[0] if site_rows else None
    if top_site: act_series=rz_series(pho.loc[top_site["site"]])
    prot_rz=robust_z(prot,gene,samp); rna_rz=robust_z(rna,gene,samp)
    # --- Check 2: downstream modules (receptor excluded) ---
    modules=module_score(gene,samp)
    modules_active=[m for m,d in modules.items() if d["max_rz"] and d["max_rz"]>1.5]
    # --- Check 3: genomic exclusion ---
    gene_mut=[ (g,pc) for (sid,g,pc,mt,vt) in onco_mut if norm_id(sid)==samp and g==gene ] + \
             [ (g,pc) for (sid,g,pc,mt) in supp_mut if norm_id(sid)==samp and g==gene ]
    focal_cn=float(cnv.at[gene,samp]) if gene in cnv.index and samp in cnv.columns else np.nan
    supp_alts=[g for (sid,g,pc,mt) in supp_mut if norm_id(sid)==samp and g in ("NF1","RASA1","PTEN","STK11","KEAP1")]
    genomic={"gene_mutation":gene_mut,"focal_CN_LR":round(focal_cn,3) if pd.notna(focal_cn) else None,
             "arm_level_CN":"not_available_open_calls","fusion_SV":"not_available_open_calls (0 SVs in study)",
             "pathway_suppressor_alts":sorted(set(supp_alts)),
             "purity":round(purity.get(samp,np.nan),3) if samp in purity else None,
             "WGD":"not_available_open_calls","total_CNA_burden":round(float(cna_burden.get(samp,np.nan)),3),
             "TMB_nonsyn":tmb.get(samp)}
    # --- Check 4: RNA exclusion ---
    lig_rna={g:round(robust_z(rna,g,samp),2) for g in LIG_OF[gene] if g in rna.index}
    rna_block={"gene_RNA_robust_z":round(rna_rz,2) if pd.notna(rna_rz) else None,
               "gene_RNA_percentile":percentile(rna,gene,samp),
               "ligand_RNA_robust_z":lig_rna,
               "note":"RNA file is tumor abundance (log2 UQ-RPKM); T/N ratio not computed (normal RNA not downloaded)"}
    # --- Check 5: stromal source ---
    act_for_corr=act_series if act_series is not None else rz_series(prot.loc[gene]) if gene in prot.index else None
    stromal={"tumor_purity":round(purity.get(samp,np.nan),3) if samp in purity else None,
             "CAF_sig_rz_in_tumor":round(sig_rz(CAF,samp),2) if pd.notna(sig_rz(CAF,samp)) else None,
             "macrophage_sig_rz_in_tumor":round(sig_rz(MAC,samp),2) if pd.notna(sig_rz(MAC,samp)) else None,
             "collagen_sig_rz_in_tumor":round(sig_rz(COLLAGEN,samp),2) if pd.notna(sig_rz(COLLAGEN,samp)) else None,
             "corr_activation_vs_CAF_acrossTumors":sig_corr(CAF,act_for_corr) if act_for_corr is not None else None,
             "corr_activation_vs_collagen_acrossTumors":sig_corr(COLLAGEN,act_for_corr) if act_for_corr is not None else None,
             "receptor_RNA_vs_purity_corr":round(float(pd.concat([rz_series(rna.loc[gene]).rename("g"),pd.Series(purity).rename("p")],axis=1).dropna().corr().iloc[0,1]),2) if gene in rna.index and len(purity)>5 else None}
    # --- Check 6: negative-regulator loss ---
    negreg={}
    for nr in NEGREG:
        p_rz=robust_z(prot,nr,samp); r_rz=robust_z(rna,nr,samp); c=float(cnv.at[nr,samp]) if nr in cnv.index and samp in cnv.columns else np.nan
        lost = (pd.notna(p_rz) and p_rz<-1) or (pd.notna(r_rz) and r_rz<-1) or (pd.notna(c) and c<-0.3)
        if lost or (pd.notna(p_rz) and abs(p_rz)>1):
            negreg[nr]={"protein_rz":round(p_rz,2) if pd.notna(p_rz) else None,"rna_rz":round(r_rz,2) if pd.notna(r_rz) else None,
                        "CN_LR":round(c,2) if pd.notna(c) else None,"flagged_loss":bool(lost)}
    # --- recurrence: other tumors with this receptor unexplained-activated ---
    recur=[]
    if act_series is not None:
        for p in pho.columns:
            az=act_series.get(p)
            if pd.notna(az) and az>2 and p!=samp:
                rz_r=robust_z(rna,gene,p); amp=float(cnv.at[gene,p]) if gene in cnv.index and p in cnv.columns else np.nan
                if not((pd.notna(rz_r) and rz_r>1.5) or (pd.notna(amp) and amp>0.3)):
                    recur.append(p)
    # --- activation source: phospho-driven vs protein-abundance-driven ---
    top_rz=top_site["robust_z"] if top_site else None
    source = "phospho_driven" if (top_rz and (prot_rz is None or (pd.notna(prot_rz) and top_rz>prot_rz+0.5))) else \
             ("protein_abundance_driven" if (pd.notna(prot_rz) and prot_rz>2) else "ambiguous")
    report[gene]={"tumor":samp,"dominant_lineage":lin,
        "check1_phosphosite":{"sites":site_rows,"top_site":top_site,"total_protein_robust_z":round(prot_rz,2) if pd.notna(prot_rz) else None,
                              "activation_source_vs_total_protein":source},
        "check2_downstream":{"modules":modules,"active_modules":modules_active,
                             "coordinated_output": len(modules_active)>=1},
        "check3_genomic":genomic,"check4_rna":rna_block,"check5_stromal":stromal,"check6_negreg":negreg,
        "recurrence_other_tumors":recur}
    json.dump(report[gene],open(f"{OUT}/validation_{gene}.json","w"),indent=1)

# ---- console summary ----
for gene,r in report.items():
    print(f"\n{'='*70}\n{gene}  (tumor {r['tumor']}, {r['dominant_lineage']})")
    ts=r["check1_phosphosite"]["top_site"]
    print(f"  [1] top phosphosite: {ts['site'] if ts else None}  robust_z={ts['robust_z'] if ts else None}  "
          f"pct={ts['percentile'] if ts else None}  missing={ts['missingness'] if ts else None}  "
          f"lineage_matched_rz={ts['lineage_matched_rz'] if ts else None}")
    print(f"      total_protein_rz={r['check1_phosphosite']['total_protein_robust_z']}  source={r['check1_phosphosite']['activation_source_vs_total_protein']}")
    print(f"  [2] active downstream modules (receptor excluded): {r['check2_downstream']['active_modules']}  coordinated={r['check2_downstream']['coordinated_output']}")
    print(f"  [3] genomic: mut={r['check3_genomic']['gene_mutation']} focalCN={r['check3_genomic']['focal_CN_LR']} "
          f"purity={r['check3_genomic']['purity']} CNAburden={r['check3_genomic']['total_CNA_burden']} suppAlts={r['check3_genomic']['pathway_suppressor_alts']}")
    print(f"  [4] RNA gene_rz={r['check4_rna']['gene_RNA_robust_z']} (pct {r['check4_rna']['gene_RNA_percentile']})  ligand_rna={r['check4_rna']['ligand_RNA_robust_z']}")
    print(f"  [5] stromal: purity={r['check5_stromal']['tumor_purity']} CAF_rz={r['check5_stromal']['CAF_sig_rz_in_tumor']} "
          f"collagen_rz={r['check5_stromal']['collagen_sig_rz_in_tumor']} corr(act,CAF)={r['check5_stromal']['corr_activation_vs_CAF_acrossTumors']} "
          f"corr(act,collagen)={r['check5_stromal']['corr_activation_vs_collagen_acrossTumors']}")
    print(f"  [6] neg-reg flagged: {[k for k,v in r['check6_negreg'].items() if v['flagged_loss']]}")
    print(f"  recurrence (other unexplained-activated tumors): {r['recurrence_other_tumors']}")
print("\nwrote candidate_cards/validation_{MET,DDR2,EGFR}.json")
