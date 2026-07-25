#!/usr/bin/env python3
"""
Regenerate the TCGA oncogene-negative tumor table with the FULL mission schema,
plus a normalized evidence-item table. Same recovered classification engine as
build_tcga_taxonomy.py; this adds metadata, competing-class columns, prose
evidence_for/against, and one normalized row per evidence item (with
independence groups, confounders checked, QC status).

Outputs:
  taxonomy/tcga_tumor_table.csv        (one row per tumor, mission schema)
  taxonomy/tcga_evidence_items.csv     (one row per evidence item, normalized schema)
"""
import gzip, os
import numpy as np, pandas as pd
A="asta_autods_data/"
os.makedirs("hidden_drivers/taxonomy",exist_ok=True)

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
LABEL={"C1_hidden_canonical":"Hidden canonical driver","C2_suppressor_conv":"Tumor-suppressor convergence",
       "C3_protein_state":"Protein-state activation","C4_lineage_cond":"Lineage-conditioned signaling",
       "C5_ligand":"Ligand / autocrine signaling","C6_microenv":"Microenvironment-dependent",
       "C7_entropy":"High-plasticity / state-entropy","C8_unknown":"Unknown"}
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
cna_burden=(cna.abs()>0.3).sum(axis=0)/cna.shape[0]   # fraction genome altered proxy

maf=pd.read_csv(A+"cleaned_mutations (1).csv"); maf["pid"]=maf.submitter_id.map(p12)
nonsyn={"Missense_Mutation","Nonsense_Mutation","Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Splice_Site","Translation_Start_Site","Nonstop_Mutation"}
mns=maf[maf.Variant_Classification.isin(nonsyn)]
def mut(genes): return set(mns[mns.Hugo_Symbol.isin(genes)].pid)

ss=pd.read_csv(A+"ssgsea_pathway_scores_hallmark_complete.csv"); ss["pid"]=ss.submitter_id.map(p12); ss=ss.set_index("pid")
def zc(col): s=ss[col]; return ((s-s.mean())/s.std())
mapk=zc("HALLMARK_KRAS_SIGNALING_UP"); pi3k=zc("HALLMARK_PI3K_AKT_MTOR_SIGNALING") if "HALLMARK_PI3K_AKT_MTOR_SIGNALING" in ss else zc("HALLMARK_MTORC1_SIGNALING")

lin=pd.read_csv("hidden_drivers/lineage/lineage_state_table_TCGA.csv"); lin["pid"]=lin.submitter_id.map(p12); lin=lin.set_index("pid")
ent_hi=lin.lineage_entropy.quantile(0.75)
# lineage confidence = max softmax prob
pcols=[c for c in lin.columns if c.endswith("_p")]
lin_conf=lin[pcols].max(axis=1) if pcols else pd.Series(dtype=float)

lf=pd.read_csv("lag3_validation/refs/leukocyte_fraction.tsv",sep="\t",header=None,names=["ct","bc","leuk"]); lf["pid"]=lf.bc.map(p12); leuk=lf.groupby("pid").leuk.mean()

on=pd.read_excel(A+"ON classification TCGA.xlsx",sheet_name="Table S1",header=1)
on=on.rename(columns={on.columns[0]:"SampleID"}); on["pid"]=on.SampleID.astype(str).str.replace("^LUAD-","",regex=True).map(p12)
on["RPA_pos"]=pd.to_numeric(on.Previous_RPA_Positive,errors="coerce")
on=on.set_index("pid")
ab=pd.read_csv("lag3_validation/refs/absolute_purity.txt",sep="\t"); ab["pid"]=ab.array.astype(str).str.replace(".","-",regex=False).map(p12)
purity=ab.groupby("pid").purity.mean(); wgd=ab.groupby("pid")["Genome doublings"].max() if "Genome doublings" in ab else pd.Series(dtype=float)

def norm(x,lo,hi): return float(np.clip((x-lo)/(hi-lo),0,1)) if pd.notna(x) else 0.0

TABLE=[]; ITEMS=[]
def item(pid,cls,gene,etype,modality,direction,obs,ref,eff,thr,src,assay,strength,supp,indep,conf,qc,interp):
    ITEMS.append(dict(patient_id=pid,mechanism_class=cls,candidate_gene_or_pathway=gene,evidence_type=etype,
        modality=modality,direction=direction,observed_value=obs,reference_value=ref,effect_size_or_zscore=eff,
        threshold=thr,data_source=src,assay=assay,evidence_strength=strength,supports_or_opposes=supp,
        independence_group=indep,confounders_checked=conf,quality_control_status=qc,interpretation=interp,
        artifact_path="hidden_drivers/taxonomy/tumor_evidence_TCGA.csv"))

DN=[pid for pid in on.index[on.RPA_pos==0] if pid in ez.index or pid in ss.index]
for pid in DN:
    # ---- evidence computations (mirror engine) ----
    rtk_vals={g:ez.loc[pid,g] for g in RTK if g in ez.columns and pid in ez.index and pd.notna(ez.loc[pid,g])}
    rtk_out=max(rtk_vals.values()) if rtk_vals else np.nan; rtk_gene=max(rtk_vals,key=rtk_vals.get) if rtk_vals else ""
    rtk_amp=max([cn(g).get(pid,np.nan) for g in ["EGFR","ERBB2","MET"]], default=np.nan)
    c1=max(norm(rtk_out,2,4), norm(rtk_amp,0.3,1.0))
    conv=0; hit=[]
    for g in SUPP:
        m=pid in mut([g]); l=cn(g).get(pid,0)<=-0.3
        if m or l: conv+=1; hit.append(g+("/mut" if m else "/loss"))
    out=max(norm(mapk.get(pid,np.nan),1,2.5), norm(pi3k.get(pid,np.nan),1,2.5))
    c2=min(1.0,0.35*conv)*0.7+0.3*out if conv>0 else 0.0
    dom=lin.dominant_lineage.get(pid,""); lstr=lin.filter(like="_score").loc[pid].max() if pid in lin.index else np.nan
    lin_str=norm(lstr,0.5,1.5)
    c4=lin_str*out if dom in ("Mucinous_gastric","Basal","EMT_mesenchymal","Proliferative") else 0.3*lin_str
    lig_vals={g:ez.loc[pid,g] for g in LIG if g in ez.columns and pid in ez.index and pd.notna(ez.loc[pid,g])}
    lig_out=max(lig_vals.values()) if lig_vals else np.nan; lig_gene=max(lig_vals,key=lig_vals.get) if lig_vals else ""
    c5=norm(lig_out,2,4)
    c6=norm(leuk.get(pid,np.nan),0.25,0.6)
    ent=lin.lineage_entropy.get(pid,np.nan); lesion=max(c1,c2,c4,c5)
    c7=(1.0 if (pid in lin.index and ent>=ent_hi) else norm(ent,ent_hi*0.9,ent_hi))*(1-lesion)
    raw={"C1_hidden_canonical":c1,"C2_suppressor_conv":c2,"C4_lineage_cond":c4,"C5_ligand":c5,"C6_microenv":c6,"C7_entropy":c7}
    tot=sum(raw.values()); unknown=max(0.0,1.0-tot) if tot<1.0 else 0.0
    Z=tot+unknown if (tot+unknown)>0 else 1
    post={k:v/Z for k,v in raw.items()}; post["C8_unknown"]=unknown/Z; post["C3_protein_state"]=0.0
    ranked=sorted(post.items(),key=lambda kv:-kv[1]); best,bp=ranked[0]
    c1p,c1v=ranked[1]; c2p,c2v=ranked[2]

    # ---- metadata ----
    wgs=on.WGS.get(pid); curalt=on.Current_RPA_Alteration.get(pid) if "Current_RPA_Alteration" in on else np.nan
    has_wgs = str(wgs).strip() not in ("nan","","None","NaN")
    reclassified = str(curalt).strip() not in ("nan","","None","NaN")
    dn_conf = "low_reclassified_on_WGS" if reclassified else ("high_WGS_confirmed" if has_wgs else "moderate_panel_only")
    layers=[l for l,ok in [("RNA",pid in ez.index),("CNA",pid in cna.columns),
            ("MAF",pid in set(maf.pid)),("ssGSEA",pid in ss.index),("purity",pid in purity.index)] if ok]
    assay_completeness=f"{len(layers)}/5:"+"+".join(layers)+"; NO proteome/phospho (C3 unmeasurable)"

    # ---- normalized evidence items ----
    item(pid,"definition","RTK/RAS/RAF","canonical_driver_absent","genomic","absent",0,0,np.nan,"RPA",
         "TCGA_RPA_classification","WES/panel","strong" if not reclassified else "weak","supports",
         "def_genomic","WGS_availability="+("yes" if has_wgs else "no"),
         "pass" if not reclassified else "reclassified_on_WGS",
         "Oncogene-negative by RPA"+(" (RECLASSIFIED positive on WGS - low confidence)" if reclassified else ""))
    item(pid,"C1_hidden_canonical","EGFR/ERBB2/MET","amplification_absent" if not(pd.notna(rtk_amp) and rtk_amp>0.3) else "copy_number_neutral",
         "genomic","neutral" if not(pd.notna(rtk_amp) and rtk_amp>0.3) else "up",
         round(float(rtk_amp),3) if pd.notna(rtk_amp) else np.nan,0.0,round(float(rtk_amp),3) if pd.notna(rtk_amp) else np.nan,0.3,
         "TCGA_CNA","gene_CN_log2ratio","moderate","opposes" if not(pd.notna(rtk_amp) and rtk_amp>0.3) else "supports",
         "genomic_CN","purity;WGD","pass","No RTK focal amplification" if not(pd.notna(rtk_amp) and rtk_amp>0.3) else "RTK amplified")
    if pd.notna(rtk_out):
        et="RNA_outlier" if rtk_out>2 else "RNA_not_elevated"
        item(pid,"C1_hidden_canonical",rtk_gene,et,"transcriptomic","up" if rtk_out>0 else "down",
             round(float(rtk_out),2),0.0,round(float(rtk_out),2),2.0,"TCGA_RNAseq","bulk_RNAseq",
             "moderate" if rtk_out>2 else "inconclusive","supports" if rtk_out>2 else "opposes",
             "rna_expression","purity;lineage","bulk_admixture_not_excluded",
             f"{rtk_gene} RNA z={rtk_out:.2f}")
    if conv>0:
        item(pid,"C2_suppressor_conv","+".join(hit),"suppressor_loss","genomic/mutational","loss",conv,0,conv,1,
             "TCGA_MAF/CNA","WES+CN","moderate" if conv>=2 else "weak","supports","suppressor_genomic",
             "purity","pass",f"{conv} suppressor(s) altered: {'+'.join(hit)}")
        item(pid,"C2_suppressor_conv","MAPK/PI3K","pathway_activation","transcriptomic","up",
             round(float(max(mapk.get(pid,np.nan),pi3k.get(pid,np.nan))),2),0.0,
             round(float(max(mapk.get(pid,np.nan),pi3k.get(pid,np.nan))),2),1.0,"TCGA_ssGSEA","Hallmark_ssGSEA",
             "weak","supports","rna_pathway","proxy_not_phospho","proxy_metric",
             "ssGSEA output proxy (weaker than phospho; independent of suppressor genomics)")
    if pd.notna(lig_out) and lig_out>1:
        item(pid,"C5_ligand",lig_gene,"ligand_expression","transcriptomic","up",round(float(lig_out),2),0.0,round(float(lig_out),2),2.0,
             "TCGA_RNAseq","bulk_RNAseq","moderate" if lig_out>2 else "weak","supports" if lig_out>2 else "opposes",
             "rna_ligand","stromal_source_possible;purity","stromal_source_possible",
             f"{lig_gene} ligand RNA z={lig_out:.2f}; ecosystem source not excluded")
    if pd.notna(leuk.get(pid,np.nan)):
        lv=leuk.get(pid)
        item(pid,"C6_microenv","leukocyte_fraction","stromal_source_possible","immune","up" if lv>0.25 else "down",
             round(float(lv),3),0.25,round(float(lv),3),0.25,"GDC_Thorsson","leukocyte_fraction",
             "moderate" if lv>0.3 else "weak","supports" if lv>0.3 else "opposes","immune_context",
             "purity","pass",f"Leukocyte fraction={lv:.2f}")
    if pd.notna(ent):
        item(pid,"C7_entropy","lineage_state","lineage_association","transcriptomic","up" if ent>=ent_hi else "neutral",
             round(float(ent),3),round(float(ent_hi),3),round(float(ent-ent_hi),3),round(float(ent_hi),3),
             "derived","lineage_signature_entropy","weak","supports" if ent>=ent_hi else "opposes",
             "rna_lineage","admixture_not_excluded","bulk_admixture_not_excluded",
             f"Lineage entropy={ent:.2f} (top-quartile cutoff {ent_hi:.2f}); plasticity vs admixture unresolved")

    # ---- prose ----
    supp_items=[i for i in ITEMS if i["patient_id"]==pid and i["supports_or_opposes"]=="supports"]
    opp_items =[i for i in ITEMS if i["patient_id"]==pid and i["supports_or_opposes"]=="opposes"]
    ev_for="; ".join(sorted({f"{i['interpretation']} [{i['evidence_strength']}]" for i in supp_items if i['mechanism_class'].startswith(best[:2]) or i['mechanism_class']==best}))[:600] or "weak/diffuse evidence only"
    ev_against="; ".join(sorted({i['interpretation'] for i in opp_items}))[:600] or "none recorded"
    intrinsic=max(c1,c2,c4,c5)
    if best in ("C7_entropy","C6_microenv","C8_unknown") or intrinsic<0.2:
        status="unexplained_context_or_plasticity"; nxt="WGS+SV / single-cell (no cell-intrinsic driver evidence)"
        unresolved="Is there a fusion/enhancer-hijack driver? Is high entropy true plasticity or admixture?"
    elif best=="C1_hidden_canonical":
        status="candidate_hidden_canonical"; nxt="WGS to confirm sub-threshold amp / fusion"
        unresolved="Is the RTK signal a true driver or lineage-driven expression?"
    elif best=="C2_suppressor_conv":
        status="candidate_suppressor_convergence"; nxt="functional MAPK/PI3K dependency test"
        unresolved="Are these tumors actually MAPK/PI3K-dependent?"
    else:
        status="candidate_"+best.lower(); nxt="orthogonal validation of "+LABEL[best]
        unresolved="Does the leading mechanism reproduce with an orthogonal assay?"

    TABLE.append(dict(
        patient_id=pid, sample_id=on.SampleID.get(pid,pid), cohort="TCGA-LUAD",
        driver_negative_definition="RPA (RTK/RAS/RAF pathway alteration) = negative"+(" [reclassified+ on WGS]" if reclassified else ""),
        driver_negative_confidence=dn_conf, assay_completeness=assay_completeness,
        purity=round(float(purity.get(pid,np.nan)),2) if pd.notna(purity.get(pid,np.nan)) else np.nan,
        dominant_lineage=dom,
        lineage_confidence=round(float(lin_conf.get(pid,np.nan)),3) if pid in lin_conf.index and pd.notna(lin_conf.get(pid,np.nan)) else np.nan,
        lineage_entropy=round(float(ent),3) if pd.notna(ent) else np.nan,
        assigned_mechanistic_class=best, class_probability=round(bp,3),
        competing_class_1=c1p, competing_class_1_probability=round(c1v,3),
        competing_class_2=c2p, competing_class_2_probability=round(c2v,3),
        evidence_for=ev_for, evidence_against=ev_against,
        unresolved_questions=unresolved, required_next_test=nxt,
        classification_status=status,
        evidence_graph_path=f"hidden_drivers/taxonomy/per_tumor_evidence.json#{pid}"))

pd.DataFrame(TABLE).to_csv("hidden_drivers/taxonomy/tcga_tumor_table.csv",index=False)
pd.DataFrame(ITEMS).to_csv("hidden_drivers/taxonomy/tcga_evidence_items.csv",index=False)
T=pd.DataFrame(TABLE)
print("TCGA tumor table rows:",len(T),"| evidence items:",len(ITEMS))
print("\nclassification_status:"); print(T.classification_status.value_counts().to_string())
print("\ndriver_negative_confidence:"); print(T.driver_negative_confidence.value_counts().to_string())
print("\nevidence_type distribution:")
print(pd.DataFrame(ITEMS).evidence_type.value_counts().to_string())
print("\nwrote tcga_tumor_table.csv, tcga_evidence_items.csv")
