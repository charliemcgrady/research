#!/usr/bin/env python3
"""
Assemble Deliverables #1-#4 data products from the two per-cohort taxonomy
tables + evidence items. Narrative deliverables (#5 failure report, #6 knowledge
gaps) are authored as markdown separately.

Outputs:
  taxonomy/mechanism_prevalence.csv        (#1/#3 cross-cohort prevalence)
  taxonomy/per_tumor_evidence.json         (#2 per-tumor evidence graph)
  taxonomy/unknown_ranking.csv             (#4 ranked unexplained tumors)
"""
import json
import numpy as np, pandas as pd

CLASSES=["C1_hidden_canonical","C2_suppressor_conv","C3_protein_state","C4_lineage_cond",
         "C5_ligand","C6_microenv","C7_entropy","C8_unknown"]
LABEL={"C1_hidden_canonical":"Hidden canonical driver","C2_suppressor_conv":"Tumor-suppressor convergence",
       "C3_protein_state":"Protein-state activation","C4_lineage_cond":"Lineage-conditioned signaling",
       "C5_ligand":"Ligand / autocrine signaling","C6_microenv":"Microenvironment-dependent",
       "C7_entropy":"High-plasticity / state-entropy","C8_unknown":"Unknown"}
INTRINSIC=["C1_hidden_canonical","C2_suppressor_conv","C3_protein_state","C4_lineage_cond","C5_ligand"]

def load(coh):
    T=pd.read_csv(f"hidden_drivers/taxonomy/tumor_evidence_{coh}.csv")
    for c in CLASSES:
        if c not in T: T[c]=0.0
    E=pd.read_csv(f"hidden_drivers/taxonomy/evidence_items_{coh}.csv")
    T["cohort"]=coh
    return T,E
Tt,Et=load("TCGA"); Tc,Ec=load("CPTAC")
T=pd.concat([Tt,Tc],ignore_index=True); E=pd.concat([Et,Ec],ignore_index=True)

# ---------- #1/#3 mechanism prevalence (per cohort) ----------
rows=[]
for c in CLASSES:
    r={"class":c,"label":LABEL[c]}
    for coh,Tx in [("TCGA",Tt),("CPTAC",Tc)]:
        n=len(Tx); best=(Tx.best_class==c).sum()
        r[f"{coh}_n_best"]=int(best); r[f"{coh}_pct_best"]=round(100*best/n,1)
        r[f"{coh}_mean_mass"]=round(Tx[c].mean(),3)
        r[f"{coh}_n_support"]=int((Tx[c]>0.1).sum())   # tumors with material evidence for this class
    r["reproducible_both"]= (r["TCGA_n_best"]>0 and r["CPTAC_n_best"]>0) or (c=="C3_protein_state")
    rows.append(r)
M=pd.DataFrame(rows)
M.to_csv("hidden_drivers/taxonomy/mechanism_prevalence.csv",index=False)

# ---------- #2 per-tumor evidence graph (JSON) ----------
graph={}
for _,t in T.iterrows():
    pid=t.pid; coh=t.cohort
    items=E[(E.pid==pid)]
    # competing classes ranked by posterior
    posts={c:float(t[c]) for c in CLASSES}
    ranked=sorted(posts.items(),key=lambda kv:-kv[1])
    evid=[]
    for _,it in items.iterrows():
        evid.append({"class":it["class"],"source":it["source"],"feature":it["feature"],
                     "value":float(it["value"]),"weight":float(it["weight"])})
    graph[pid]={"cohort":coh,"best_class":t.best_class,"best_p":float(t.best_p),
                "dominant_lineage":t.dominant_lineage,"lineage_entropy":float(t.lineage_entropy) if pd.notna(t.lineage_entropy) else None,
                "posteriors":{k:round(v,3) for k,v in posts.items()},
                "competing":[{"class":k,"label":LABEL[k],"p":round(v,3)} for k,v in ranked if v>0.02],
                "evidence":sorted(evid,key=lambda e:-e["weight"])}
json.dump(graph,open("hidden_drivers/taxonomy/per_tumor_evidence.json","w"),indent=1)

# ---------- #4 unknown / unexplained ranking ----------
def intrinsic_expl(row): return max(row[c] for c in INTRINSIC)
recs=[]
for _,t in T.iterrows():
    ie=intrinsic_expl(t)
    # "unexplained" = little cell-intrinsic-driver evidence; best explanation is context/plasticity/unknown
    unexplained = 1.0 - ie
    recs.append({"pid":t.pid,"cohort":t.cohort,"best_class":t.best_class,"best_p":t.best_p,
                 "intrinsic_explanation":round(ie,3),"unexplained_score":round(unexplained,3),
                 "dominant_lineage":t.dominant_lineage,
                 "lineage_entropy":t.lineage_entropy,
                 "C6_microenv":t.C6_microenv,"C7_entropy":t.C7_entropy,"C8_unknown":t.C8_unknown})
U=pd.DataFrame(recs).sort_values(["unexplained_score","best_p"],ascending=[False,True])
U.to_csv("hidden_drivers/taxonomy/unknown_ranking.csv",index=False)

# ---------- console summary ----------
print("=== Deliverable #1/#3: mechanism prevalence ===")
print(M[["label","TCGA_n_best","TCGA_pct_best","TCGA_mean_mass","CPTAC_n_best","CPTAC_pct_best","CPTAC_mean_mass","reproducible_both"]].to_string(index=False))
print("\n=== Deliverable #4: top-15 most-unexplained tumors ===")
print(U.head(15)[["pid","cohort","best_class","best_p","intrinsic_explanation","dominant_lineage"]].to_string(index=False))
print("\ncohort totals: TCGA neg n=%d, CPTAC neg n=%d"%(len(Tt),len(Tc)))
print("tumors with cell-intrinsic explanation >=0.35: TCGA %d/%d, CPTAC %d/%d"%(
    (Tt[INTRINSIC].max(axis=1)>=0.35).sum(),len(Tt),(Tc[INTRINSIC].max(axis=1)>=0.35).sum(),len(Tc)))
print("wrote mechanism_prevalence.csv, per_tumor_evidence.json, unknown_ranking.csv")
