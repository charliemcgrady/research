#!/usr/bin/env python3
"""
CPTAC oncogene-negative tumor table in the full mission schema (symmetric with
tcga_tumor_table.csv), integrating the validated Class-3 verdicts for the three
candidates. Reuses the CPTAC engine outputs; adds metadata + competing classes +
prose + validated classification_status.
Writes: taxonomy/cptac_tumor_table.csv
"""
import json, csv
import numpy as np, pandas as pd
CB="cptac_luad/"
LABEL={"C1_hidden_canonical":"Hidden canonical driver","C2_suppressor_conv":"Tumor-suppressor convergence",
       "C3_protein_state":"Protein-state activation","C4_lineage_cond":"Lineage-conditioned signaling",
       "C5_ligand":"Ligand / autocrine signaling","C6_microenv":"Microenvironment-dependent",
       "C7_entropy":"High-plasticity / state-entropy","C8_unknown":"Unknown"}
CLASSES=list(LABEL)
def norm_id(s): return s.lstrip("X").replace(".","-")

T=pd.read_csv("hidden_drivers/taxonomy/tumor_evidence_CPTAC.csv").set_index("pid")
E=pd.read_csv("hidden_drivers/taxonomy/evidence_items_CPTAC.csv")
clin=json.load(open(CB+"cptac_clinical_all.json"))
purity={norm_id(k):v for k,v in clin.get("TSNET PURITY",{}).items() if v not in ("NA","","[Not Available]")}
tmb={norm_id(k):v for k,v in clin.get("TMB_NONSYNONYMOUS",{}).items() if v not in ("NA","",None)}
onco_calls=pd.read_csv("hidden_drivers/taxonomy/cptac_onco_positive_calls.csv") if False else None
verd=pd.read_csv("hidden_drivers/candidate_cards/class3_verdicts.csv").set_index("tumor")

rows=[]
for pid,t in T.iterrows():
    posts={c:float(t.get(c,0.0)) for c in CLASSES}
    ranked=sorted(posts.items(),key=lambda kv:-kv[1])
    best,bp=ranked[0]; c1p,c1v=ranked[1]; c2p,c2v=ranked[2]
    items=E[E.pid==pid]
    supp=items[items.weight>0.1]
    ev_for="; ".join(f"{r['source']}:{r['feature']} (w={r['weight']})" for _,r in supp.sort_values('weight',ascending=False).head(4).iterrows()) or "diffuse/weak"
    # validated status override for the 3 candidates
    if pid in verd.index:
        status=verd.loc[pid,"validated_status"]; nxt=verd.loc[pid,"required_next_test"]
        ev_against="see candidate_cards/CANDIDATE_%s.md"%verd.loc[pid,"candidate_gene"]
        unresolved="Class-3 candidate validated: "+status
        graph="hidden_drivers/candidate_cards/CANDIDATE_%s.md"%verd.loc[pid,"candidate_gene"]
    else:
        intrinsic=max(posts["C1_hidden_canonical"],posts["C2_suppressor_conv"],posts["C3_protein_state"],posts["C4_lineage_cond"],posts["C5_ligand"])
        if best in ("C7_entropy","C6_microenv","C8_unknown") or intrinsic<0.2:
            status="unexplained_context_or_plasticity"; nxt="WGS+SV / single-cell"
        else: status="candidate_"+best.lower(); nxt="orthogonal validation of "+LABEL[best]
        ev_against="none recorded"; unresolved="Does the leading mechanism reproduce orthogonally?"
        graph=f"hidden_drivers/taxonomy/per_tumor_evidence.json#{pid}"
    layers="5/5:RNA+CNV+MAF+proteome+phospho (C3 measurable)"
    rows.append(dict(
        patient_id=pid, sample_id=pid, cohort="CPTAC-LUAD",
        driver_negative_definition="RTK/RAS/RAF WT (no activating mut/amp; fusion-RTK RNA-outlier proxy)",
        driver_negative_confidence="moderate_no_SV_layer (fusions untestable; 0 SVs in study)",
        assay_completeness=layers,
        purity=round(float(purity[pid]),2) if pid in purity else np.nan,
        dominant_lineage=t.dominant_lineage,
        lineage_confidence=np.nan,
        lineage_entropy=round(float(t.lineage_entropy),3) if pd.notna(t.lineage_entropy) else np.nan,
        assigned_mechanistic_class=best, class_probability=round(bp,3),
        competing_class_1=c1p, competing_class_1_probability=round(c1v,3),
        competing_class_2=c2p, competing_class_2_probability=round(c2v,3),
        evidence_for=ev_for, evidence_against=ev_against,
        unresolved_questions=unresolved, required_next_test=nxt,
        classification_status=status, evidence_graph_path=graph))
out=pd.DataFrame(rows)
out.to_csv("hidden_drivers/taxonomy/cptac_tumor_table.csv",index=False)
print("wrote cptac_tumor_table.csv (n=%d)"%len(out))
print(out.classification_status.value_counts().to_string())
print("\nvalidated Class-3 candidate rows:")
print(out[out.patient_id.isin(verd.index)][["patient_id","assigned_mechanistic_class","classification_status"]].to_string(index=False))
