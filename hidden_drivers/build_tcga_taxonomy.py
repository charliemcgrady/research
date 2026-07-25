#!/usr/bin/env python3
"""
Mechanistic taxonomy engine — TCGA-LUAD cohort (bulk multi-omics).
Recovered/cleaned from the interactive engine. Emits per-tumor posteriors AND
explicit evidence items (for the per-tumor evidence-graph deliverable).

Class 3 (protein-state activation) is NOT scorable here: TCGA in this project
has no proteome/phosphoproteome. C3 is a CPTAC-only column, reported per-cohort.

Shared with CPTAC engine: gene sets, soft evidence combiner, lineage model.
Pathway output here is an ssGSEA proxy (HALLMARK_KRAS_UP / PI3K_AKT_MTOR),
a known weaker surrogate than CPTAC's phospho readout.
"""
import gzip
import numpy as np, pandas as pd
A="asta_autods_data/"

# ---- symbol->ENSG ----
sym2ensg={}
with gzip.open("lag3_validation/scrna/luca/gencode.basic.gtf.gz","rt") as f:
    for line in f:
        if line.startswith("#"):continue
        p=line.split("\t")
        if len(p)>8 and p[2]=="gene":
            gid=nm=None
            for fld in p[8].split(";"):
                fld=fld.strip()
                if fld.startswith("gene_id"):gid=fld.split('"')[1].split(".")[0]
                elif fld.startswith("gene_name"):nm=fld.split('"')[1]
            if gid and nm:sym2ensg.setdefault(nm,gid)

RTK=["EGFR","ERBB2","ERBB3","MET","ALK","ROS1","RET","FGFR1","FGFR2","FGFR3","FGFR4","NTRK1","NTRK2","NTRK3","AXL","IGF1R","DDR2","KIT","PDGFRA","PDGFRB"]
LIG=["NRG1","NRG2","HGF","EGF","TGFA","AREG","EREG","HBEGF","FGF2","FGF9","FGF10","IGF1","IGF2","PDGFA","PDGFB","BTC"]
SUPP=["NF1","RASA1","PTEN","STK11","KEAP1"]
need=set(RTK)|set(LIG); needensg={sym2ensg[g] for g in need if g in sym2ensg}

with gzip.open(A+"cleaned_gene_expression.csv.gz","rt") as f:
    hdr=f.readline().rstrip("\n").split(","); samples=hdr[1:]; rows={}
    for line in f:
        b=line.split(",",1)[0].split(".")[0]
        if b in needensg: rows[b]=pd.to_numeric(pd.Series(line.rstrip("\n").split(",")[1:]),errors="coerce").values
e=pd.DataFrame(rows,index=samples); e.columns=[{v:k for k,v in sym2ensg.items()}.get(c,c) for c in e.columns]
ez=(e-e.mean())/e.std()
def p12(x):return str(x)[:12]
ez.index=[p12(s) for s in ez.index]; ez=ez[~ez.index.duplicated()]

cna=pd.read_csv(A+"cleaned_cna (1).csv").set_index("Gene Symbol"); cna.columns=[p12(c) for c in cna.columns]
def cn(gene): return cna.loc[gene] if gene in cna.index else pd.Series(dtype=float)

maf=pd.read_csv(A+"cleaned_mutations (1).csv"); maf["pid"]=maf.submitter_id.map(p12)
nonsyn={"Missense_Mutation","Nonsense_Mutation","Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Splice_Site","Translation_Start_Site","Nonstop_Mutation"}
mns=maf[maf.Variant_Classification.isin(nonsyn)]
def mut(genes): return set(mns[mns.Hugo_Symbol.isin(genes)].pid)

ss=pd.read_csv(A+"ssgsea_pathway_scores_hallmark_complete.csv"); ss["pid"]=ss.submitter_id.map(p12); ss=ss.set_index("pid")
def zc(col): s=ss[col]; return ((s-s.mean())/s.std())
mapk=zc("HALLMARK_KRAS_SIGNALING_UP"); pi3k=zc("HALLMARK_PI3K_AKT_MTOR_SIGNALING") if "HALLMARK_PI3K_AKT_MTOR_SIGNALING" in ss else zc("HALLMARK_MTORC1_SIGNALING")

lin=pd.read_csv("hidden_drivers/lineage/lineage_state_table_TCGA.csv"); lin["pid"]=lin.submitter_id.map(p12); lin=lin.set_index("pid")
ent_hi=lin.lineage_entropy.quantile(0.75)

lf=pd.read_csv("lag3_validation/refs/leukocyte_fraction.tsv",sep="\t",header=None,names=["ct","bc","leuk"]); lf["pid"]=lf.bc.map(p12); leuk=lf.groupby("pid").leuk.mean()

on=pd.read_excel(A+"ON classification TCGA.xlsx",sheet_name="Table S1",header=1)
on=on.rename(columns={on.columns[0]:"SampleID"}); on["pid"]=on.SampleID.astype(str).str.replace("^LUAD-","",regex=True).map(p12)
on["RPA_pos"]=pd.to_numeric(on.Previous_RPA_Positive,errors="coerce")
ab=pd.read_csv("lag3_validation/refs/absolute_purity.txt",sep="\t"); ab["pid"]=ab.array.astype(str).str.replace(".","-",regex=False).map(p12); purity=ab.groupby("pid").purity.mean()

def norm(x,lo,hi): return float(np.clip((x-lo)/(hi-lo),0,1)) if pd.notna(x) else 0.0
def score(pid, sign, items=None):
    def add(cls,src,feat,val,w):
        if items is not None and w>0.05:
            items.append({"pid":pid,"class":cls,"source":src,"feature":feat,"value":round(float(val),3),"weight":round(float(w),3)})
    rtk_out=max([ez.loc[pid,g] for g in RTK if g in ez.columns and pid in ez.index], default=np.nan)
    rtk_amp=max([cn(g).get(pid,np.nan) for g in ["EGFR","ERBB2","MET"]], default=np.nan)
    c1=max(norm(rtk_out,2,4), norm(rtk_amp,0.3,1.0))
    add("C1_hidden_canonical","RNA","RTK_expr_outlier",rtk_out if pd.notna(rtk_out) else 0,norm(rtk_out,2,4))
    add("C1_hidden_canonical","CNV","RTK_amp",rtk_amp if pd.notna(rtk_amp) else 0,norm(rtk_amp,0.3,1.0))
    conv=0; hit=[]
    for g in SUPP:
        m=pid in mut([g]); l=cn(g).get(pid,0)<=-0.3
        if m or l: conv+=1; hit.append(g+("/mut" if m else "/loss"))
    out=max(norm(mapk.get(pid,np.nan),1,2.5), norm(pi3k.get(pid,np.nan),1,2.5))
    c2=min(1.0,0.35*conv)*0.7+0.3*out if conv>0 else 0.0
    if conv>0:
        add("C2_suppressor_conv","MUT/CNV","suppressor_hits({})".format("+".join(hit)),conv,min(1.0,0.35*conv)*0.7)
        add("C2_suppressor_conv","ssGSEA","pathway_output",max(mapk.get(pid,0),pi3k.get(pid,0)),0.3*out)
    dom=lin.dominant_lineage.get(pid,"")
    lin_str=norm(lin.filter(like="_score").loc[pid].max() if pid in lin.index else np.nan,0.5,1.5) if pid in lin.index else 0
    c4=lin_str*out if dom in ("Mucinous_gastric","Basal","EMT_mesenchymal","Proliferative") else 0.3*lin_str
    add("C4_lineage_cond","RNA",f"dom={dom}xoutput",lin_str,c4)
    lig_out=max([ez.loc[pid,g] for g in LIG if g in ez.columns and pid in ez.index], default=np.nan)
    c5=norm(lig_out,2,4); add("C5_ligand","RNA","ligand_outlier",lig_out if pd.notna(lig_out) else 0,c5)
    c6=norm(leuk.get(pid,np.nan),0.25,0.6); add("C6_microenv","GDC","leukocyte_fraction",leuk.get(pid,np.nan) if pd.notna(leuk.get(pid,np.nan)) else 0,c6)
    ent=lin.lineage_entropy.get(pid,np.nan)
    lesion=max(c1,c2,c4,c5)
    c7=(1.0 if (pid in lin.index and ent>=ent_hi) else norm(ent,ent_hi*0.9,ent_hi))*(1-lesion)
    add("C7_entropy","RNA","lineage_entropy",ent if pd.notna(ent) else 0,c7)
    raw={"C1_hidden_canonical":c1,"C2_suppressor_conv":c2,"C4_lineage_cond":c4,"C5_ligand":c5,"C6_microenv":c6,"C7_entropy":c7}
    tot=sum(raw.values()); unknown=max(0.0,1.0-tot) if tot<1.0 else 0.0
    Z=tot+unknown if (tot+unknown)>0 else 1
    post={k:v/Z for k,v in raw.items()}; post["C8_unknown"]=unknown/Z
    best=max(post,key=post.get)
    return {"pid":pid,"best_class":best,"best_p":round(post[best],2),"dominant_lineage":dom,
            "lineage_entropy":round(ent,3) if pd.notna(ent) else np.nan,
            "onco_status":"neg" if sign==0 else "pos",
            "purity":round(purity.get(pid,np.nan),2) if pd.notna(purity.get(pid,np.nan)) else np.nan,
            **{k:round(v,3) for k,v in post.items()}}

DN=[pid for pid in on[on.RPA_pos==0].pid if pid in ez.index or pid in ss.index]
POS=[pid for pid in on[on.RPA_pos==1].pid if pid in ez.index or pid in ss.index]
items=[]
Trows=[score(p,0,items) for p in DN]; Prows=[score(p,1,None) for p in POS]
T=pd.DataFrame(Trows).set_index("pid"); P=pd.DataFrame(Prows).set_index("pid")
T.reset_index().to_csv("hidden_drivers/taxonomy/tumor_evidence_TCGA.csv",index=False)
pd.DataFrame(items).to_csv("hidden_drivers/taxonomy/evidence_items_TCGA.csv",index=False)
print("TCGA oncogene-negative n=%d"%len(T))
print(T.best_class.value_counts().to_string())
print("\nmean posterior mass:"); print(T[[c for c in T.columns if c.startswith("C")]].mean().round(3).sort_values(ascending=False).to_string())
print("\nPOSITIVE CONTROL oncogene-positive n=%d (expect C1 lead):"%len(P))
print(P.best_class.value_counts().to_string())
print("mean C1 pos=%.3f neg=%.3f"%(P.C1_hidden_canonical.mean(),T.C1_hidden_canonical.mean()))
print("wrote tumor_evidence_TCGA.csv, evidence_items_TCGA.csv")
