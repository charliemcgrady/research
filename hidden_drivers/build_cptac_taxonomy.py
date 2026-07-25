#!/usr/bin/env python3
"""
Mechanistic taxonomy engine — CPTAC-LUAD cohort (protein/phospho-enabled).

Adds Class 3 (protein-state activation), which is only measurable where
proteome+phosphoproteome exist. Shared classes (C1,C2,C4,C5,C6,C7,C8) mirror
the TCGA engine, but pathway *output* here is a real phospho-signaling readout
(ERK activation-loop, AKT, RPS6) rather than an ssGSEA proxy.

Inputs (open CPTAC-LUAD, LinkedOmics + cBioPortal):
  raw/HS_CPTAC_LUAD_{proteome,phosphoproteome,rnaseq,cnv}_...cct
  raw/HS_CPTAC_LUAD_somatic_mutation_gene.cbt
  cptac_onco_mutations.json, cptac_suppressor_mutations.json  (cBioPortal, protein-change annotated)

Philosophy (per mission): NOT p-value optimization. Soft, evidence-weighted
per-tumor posteriors; multiple weak convergent signals beat one strong RNA hit.
"""
import csv, json, math
import numpy as np, pandas as pd

RAW="cptac_luad/raw/"
CB ="cptac_luad/"

# ---------------- IO helpers ----------------
def norm_id(s): return s.lstrip("X").replace(".","-")
def load_cct(fn, is_phospho=False):
    with open(RAW+fn) as f:
        r=csv.reader(f,delimiter="\t"); hdr=next(r); samples=[norm_id(x) for x in hdr[1:]]
        idx=[]; data=[]
        for row in r:
            idx.append(row[0])
            data.append([np.nan if v in ("NA","","NaN") else float(v) for v in row[1:]])
    df=pd.DataFrame(data,index=idx,columns=samples)
    return df
def zrow(df):  # z-score each feature (row) across samples
    return df.sub(df.mean(axis=1),axis=0).div(df.std(axis=1).replace(0,np.nan),axis=0)

# ---------------- load layers ----------------
prot=load_cct("HS_CPTAC_LUAD_proteome_ratio_NArm_TUMOR.cct")
pho =load_cct("HS_CPTAC_LUAD_phosphoproteome_ratio_norm_NArm_TUMOR.cct")
rna =load_cct("HS_CPTAC_LUAD_rnaseq_uq_rpkm_log2_NArm_TUMOR.cct")
cnv =load_cct("HS_CPTAC_LUAD_cnv_gene_LR.cct")
# gene-level binary mutation matrix (suppressor presence fallback)
mutm=load_cct("HS_CPTAC_LUAD_somatic_mutation_gene.cbt")

samples=sorted(set(prot.columns)&set(pho.columns)&set(rna.columns)&set(cnv.columns))
protZ=zrow(prot); rnaZ=zrow(rna); phoZ=zrow(pho)

# ---------------- gene sets ----------------
RTK=["EGFR","ERBB2","ERBB3","MET","ALK","ROS1","RET","FGFR1","FGFR2","FGFR3","FGFR4",
     "NTRK1","NTRK2","NTRK3","AXL","IGF1R","DDR2","KIT","PDGFRA","PDGFRB"]
FUSION_RTK=["ALK","ROS1","RET","NTRK1","NTRK2","NTRK3"]     # activated by fusion (no point-mut/amp signal)
LIG=["NRG1","NRG2","HGF","EGF","TGFA","AREG","EREG","HBEGF","FGF2","FGF9","FGF10","IGF1","IGF2","PDGFA","PDGFB","BTC"]
SUPP=["NF1","RASA1","PTEN","STK11","KEAP1"]
# activation-loop / autophospho readouts
ERK_SITES=["MAPK1:NP_002736.3:T185tY187y","MAPK3:NP_002737.2:T202tY204y"]
AKT_SITES=["AKT1:NP_001014431.1:S477s","AKT2:NP_001617.1:S474s"]
S6_SITES =["RPS6:NP_001001.2:S235sS236s","RPS6:NP_001001.2:S240s"]
# epithelial lineage signatures (candidate-driver genes excluded; same as TCGA model)
SIG={
 "AT2_TRU":["SFTPC","SFTPB","SFTPA1","SFTPA2","NAPSA","LAMP3","SFTPD","PGC","SLC34A2","ETV5","LPCAT1","ABCA3"],
 "AT1":["AGER","PDPN","CAV1","CAV2","RTKN2","CLIC5","AKAP5","CLDN18","SPOCK2"],
 "Club_secretory":["SCGB1A1","SCGB3A2","SCGB3A1","MUC5B","BPIFB1","CYP2F1","WFDC2"],
 "Ciliated":["FOXJ1","TPPP3","PIFO","SNTN","CAPS","DNAI1","CCDC78"],
 "Basal":["KRT5","KRT14","TP63","KRT6A","DAPL1","KRT17","SERPINB3","KRT15"],
 "Mucinous_gastric":["HNF4A","FOXA3","TFF1","TFF2","TFF3","MUC5AC","GKN1","GKN2","CDX2","MUC2","SPINK1"],
 "Proliferative":["MKI67","TOP2A","CCNB1","CDK1","BIRC5","UBE2C","CENPF","PCNA"],
 "EMT_mesenchymal":["VIM","ZEB1","ZEB2","SNAI2","TWIST1","FN1","CDH2","SPARC"],
}
IMMUNE=["PTPRC","CD3D","CD3E","CD8A","CD2","CD68","CD163","MS4A1","IL2RG","LCK","GZMB","NKG7"]
STROMAL=["FAP","COL1A1","COL1A2","COL3A1","PDGFRB","THY1","ACTA2","LUM","DCN"]

# ---------------- mutation calls (activating) ----------------
onco_mut=json.load(open(CB+"cptac_onco_mutations.json"))       # [sampleId, gene, protChange, mutType, varType]
supp_mut=json.load(open(CB+"cptac_suppressor_mutations.json")) # [sampleId, gene, protChange, mutType]
def _codon(pc):
    import re
    m=re.match(r"[A-Z](\d+)",pc or ""); return int(m.group(1)) if m else None
# LUAD driver genes where any mutation is treated activating; fusion RTKs handled via RNA outlier
ACT_GENES={"EGFR","KRAS","NRAS","HRAS","BRAF","ERBB2","MAP2K1","RIT1"}
onco_pos=set()          # sample -> oncogene-positive
onco_reason={}
for s,g,pc,mt,vt in onco_mut:
    s=norm_id(s)
    if g in ACT_GENES:
        onco_pos.add(s); onco_reason.setdefault(s,[]).append(f"{g}:{pc}")
# focal amplification positivity (LR>0.6) for classic amp drivers
for g in ["EGFR","ERBB2","MET","FGFR1"]:
    if g in cnv.index:
        for s in samples:
            v=cnv.at[g,s]
            if pd.notna(v) and v>0.6:
                onco_pos.add(s); onco_reason.setdefault(s,[]).append(f"{g}:amp(LR={v:.2f})")
# fusion-RTK RNA outlier (z>3) -> treat as likely fusion-driven -> exclude from clean negative set
for g in FUSION_RTK:
    if g in rnaZ.index:
        for s in samples:
            v=rnaZ.at[g,s]
            if pd.notna(v) and v>3:
                onco_pos.add(s); onco_reason.setdefault(s,[]).append(f"{g}:RNA_outlier(z={v:.1f})")
# MET exon14 / high MET from binary matrix + MET amp already covered; also count binary matrix MET splice-like? skip (no splice detail)

NEG=[s for s in samples if s not in onco_pos]
print(f"CPTAC-LUAD: {len(samples)} tumors with full multi-omics; oncogene-POSITIVE={len(onco_pos)}, "
      f"comprehensively oncogene-NEGATIVE={len(NEG)}")

# suppressor alteration set
supp_mut_by={g:set() for g in SUPP}
for s,g,pc,mt in supp_mut:
    if g in SUPP: supp_mut_by[g].add(norm_id(s))
# also binary matrix
for g in SUPP:
    if g in mutm.index:
        for s in samples:
            if mutm.at[g,s]==1: supp_mut_by[g].add(s)

# ---------------- feature helpers ----------------
def gz(df, gene, s):
    return df.at[gene,s] if gene in df.index and s in df.columns else np.nan
def maxz(df, genes, s):
    vs=[df.at[g,s] for g in genes if g in df.index and pd.notna(df.at[g,s])]
    return max(vs) if vs else np.nan
def phos_act(sites, s):
    vs=[phoZ.at[k,s] for k in sites if k in phoZ.index and pd.notna(phoZ.at[k,s])]
    return max(vs) if vs else np.nan
def sig_score(genes, s, Z):
    vs=[Z.at[g,s] for g in genes if g in Z.index and pd.notna(Z.at[g,s])]
    return float(np.mean(vs)) if vs else np.nan
def nz(x,lo,hi):
    return float(np.clip((x-lo)/(hi-lo),0,1)) if pd.notna(x) else 0.0

# pathway phospho output per sample
erk={s:phos_act(ERK_SITES,s) for s in samples}
akt={s:phos_act(AKT_SITES,s) for s in samples}
s6 ={s:phos_act(S6_SITES,s)  for s in samples}
def mapk_out(s): return erk[s]
def pi3k_out(s): return np.nanmax([akt[s], s6[s]]) if not (pd.isna(akt[s]) and pd.isna(s6[s])) else np.nan

# lineage scoring from RNA (z of signature means, softmax)
lin_scores={s:{k:sig_score(v,s,rnaZ) for k,v in SIG.items()} for s in samples}
def lineage(s):
    d=lin_scores[s]; ks=list(d); arr=np.array([d[k] if pd.notna(d[k]) else 0 for k in ks])
    ex=np.exp(arr-arr.max()); p=ex/ex.sum()
    dom=ks[int(p.argmax())]; ent=float(-(p*np.log(p+1e-12)).sum()); strength=float(arr.max())
    return dom,strength,ent
ent_all=[lineage(s)[2] for s in samples]; ent_hi=np.quantile(ent_all,0.75)

# ---------------- score one tumor ----------------
def score(s, record_items=None):
    ev={}; items=[]
    def add(cls,src,feat,val,w):
        if record_items is not None and w>0.05:
            items.append({"pid":s,"class":cls,"source":src,"feature":feat,"value":round(float(val),3),"weight":round(float(w),3)})
    # ---- pathway output (phospho) ----
    mo=mapk_out(s); po=pi3k_out(s); out=max(nz(mo,1,2.5),nz(po,1,2.5))
    # ---- C1 hidden canonical: sub-threshold RTK amp OR RTK RNA/protein outlier ----
    rtk_rna=maxz(rnaZ,RTK,s); rtk_amp=max([gz(cnv,g,s) for g in ["EGFR","ERBB2","MET","FGFR1"] if pd.notna(gz(cnv,g,s))],default=np.nan)
    c1=max(nz(rtk_rna,2,4), nz(rtk_amp,0.3,0.6))
    add("C1_hidden_canonical","RNA","RTK_expr_outlier",rtk_rna if pd.notna(rtk_rna) else 0,nz(rtk_rna,2,4))
    add("C1_hidden_canonical","CNV","RTK_subthreshold_amp",rtk_amp if pd.notna(rtk_amp) else 0,nz(rtk_amp,0.3,0.6))
    # ---- C2 suppressor convergence ----
    conv=0; hit=[]
    for g in SUPP:
        m = s in supp_mut_by[g]; l = gz(cnv,g,s)<=-0.3 if pd.notna(gz(cnv,g,s)) else False
        if m or l: conv+=1; hit.append(g+("/mut" if m else "/loss"))
    c2=(min(1.0,0.35*conv)*0.7 + 0.3*out) if conv>0 else 0.0
    if conv>0:
        add("C2_suppressor_conv","MUT/CNV","suppressor_hits({})".format("+".join(hit)),conv,min(1.0,0.35*conv)*0.7)
        add("C2_suppressor_conv","PHOSPHO","pathway_output(ERK/AKT/S6)",max(mo if pd.notna(mo) else 0,po if pd.notna(po) else 0),0.3*out)
    # ---- C3 protein-state activation (CPTAC-only) ----
    # (a) RTK activated at protein/phospho level, unexplained by DNA(amp/mut) or RNA
    c3=0.0; c3feat=""
    for g in RTK:
        pa=phos_act([k for k in phoZ.index if k.startswith(g+":")],s) if any(k.startswith(g+":") for k in phoZ.index) else np.nan
        pr_z=gz(protZ,g,s); rna_z=gz(rnaZ,g,s); amp=gz(cnv,g,s)
        act=np.nanmax([pa,pr_z]) if not(pd.isna(pa) and pd.isna(pr_z)) else np.nan
        dna_expl = (pd.notna(amp) and amp>0.3)
        rna_expl = (pd.notna(rna_z) and rna_z>1.5)
        if pd.notna(act) and act>2 and not dna_expl and not rna_expl:
            v=nz(act,2,4)
            if v>c3: c3=v; c3feat=f"{g}:activation_z={act:.1f}(protein/phospho),RNA_z={rna_z if pd.notna(rna_z) else float('nan'):.1f},amp={amp if pd.notna(amp) else float('nan'):.2f}"
    # (b) downstream ERK/AKT/S6 output high with NO upstream DNA/RNA lesion at all
    upstream_lesion = (pd.notna(rtk_amp) and rtk_amp>0.3) or (pd.notna(rtk_rna) and rtk_rna>1.5) or conv>0
    out_z=max(mo if pd.notna(mo) else -9, po if pd.notna(po) else -9)
    if out_z>2 and not upstream_lesion:
        v=nz(out_z,2,4)
        if v>c3: c3=v; c3feat=f"downstream_output_z={out_z:.1f}_no_upstream_lesion"
    if c3>0.05: add("C3_protein_state","PHOSPHO/PROTEIN",c3feat,c3,c3)
    # ---- C4 lineage conditioned ----
    dom,strength,ent=lineage(s)
    lin_str=nz(strength,0.5,1.5)
    c4=lin_str*out if dom in ("Mucinous_gastric","Basal","EMT_mesenchymal","Proliferative") else 0.3*lin_str
    add("C4_lineage_cond","RNA",f"dom={dom}(strength={strength:.2f})xoutput",strength,c4)
    # ---- C5 ligand ----
    lig_rna=maxz(rnaZ,LIG,s); lig_prot=maxz(protZ,LIG,s)
    lig=np.nanmax([lig_rna,lig_prot]) if not(pd.isna(lig_rna) and pd.isna(lig_prot)) else np.nan
    c5=nz(lig,2,4); add("C5_ligand","RNA/PROTEIN","ligand_outlier",lig if pd.notna(lig) else 0,c5)
    # ---- C6 microenvironment ----
    imm=sig_score(IMMUNE,s,rnaZ); strom=sig_score(STROMAL,s,rnaZ)
    micro=np.nanmax([imm,strom]) if not(pd.isna(imm) and pd.isna(strom)) else np.nan
    c6=nz(micro,0.5,1.5); add("C6_microenv","RNA","immune/stromal_signature",micro if pd.notna(micro) else 0,c6)
    # ---- C7 entropy (residual) ----
    lesion=max(c1,c2,c3,c4,c5)
    c7=(1.0 if ent>=ent_hi else nz(ent,ent_hi*0.9,ent_hi))*(1-lesion)
    add("C7_entropy","RNA","lineage_entropy",ent,c7)
    raw={"C1_hidden_canonical":c1,"C2_suppressor_conv":c2,"C3_protein_state":c3,
         "C4_lineage_cond":c4,"C5_ligand":c5,"C6_microenv":c6,"C7_entropy":c7}
    tot=sum(raw.values()); unknown=max(0.0,1.0-tot) if tot<1.0 else 0.0
    Z=tot+unknown if (tot+unknown)>0 else 1
    post={k:v/Z for k,v in raw.items()}; post["C8_unknown"]=unknown/Z
    best=max(post,key=post.get)
    row={"pid":s,"best_class":best,"best_p":round(post[best],2),"dominant_lineage":dom,
         "lineage_entropy":round(ent,3),"onco_status":"neg" if s in NEG else "pos",
         "erk_pTpY_z":round(mo,2) if pd.notna(mo) else np.nan,
         **{k:round(v,3) for k,v in post.items()}}
    if record_items is not None: record_items.extend(items)
    return row

# ---------------- run on oncogene-negative + positive control ----------------
items=[]
neg_rows=[score(s,items) for s in NEG]
pos_rows=[score(s,None) for s in sorted(onco_pos)]
T=pd.DataFrame(neg_rows).set_index("pid")
P=pd.DataFrame(pos_rows).set_index("pid")

import os
os.makedirs("hidden_drivers/taxonomy",exist_ok=True)
T.reset_index().to_csv("hidden_drivers/taxonomy/tumor_evidence_CPTAC.csv",index=False)
pd.DataFrame(items).to_csv("hidden_drivers/taxonomy/evidence_items_CPTAC.csv",index=False)
# reasons for onco-positive calls (audit trail)
pd.DataFrame([{"pid":norm_id(s),"onco_reason":"; ".join(sorted(set(r)))} for s,r in onco_reason.items()]
             ).to_csv("hidden_drivers/taxonomy/cptac_onco_positive_calls.csv",index=False)

print("\n=== CPTAC oncogene-NEGATIVE mechanistic taxonomy (n=%d) ==="%len(T))
print(T.best_class.value_counts().to_string())
cols=[c for c in T.columns if c.startswith("C")]
print("\nmean posterior mass per class:")
print(T[cols].mean().round(3).sort_values(ascending=False).to_string())
print("\n=== POSITIVE CONTROL: oncogene-POSITIVE tumors (n=%d) — expect C1 to lead ==="%len(P))
print(P.best_class.value_counts().to_string())
print("mean C1 (pos)=%.3f vs (neg)=%.3f | mean C3 (pos)=%.3f vs (neg)=%.3f"%(
    P.C1_hidden_canonical.mean(),T.C1_hidden_canonical.mean(),P.C3_protein_state.mean(),T.C3_protein_state.mean()))
print("\nwrote tumor_evidence_CPTAC.csv, evidence_items_CPTAC.csv, cptac_onco_positive_calls.csv")
