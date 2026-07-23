#!/usr/bin/env python3
"""
Pipeline smoke test: synthetic multi-tumor scRNA with a KNOWN injected tumor-level
CIN->CD8-LAG3 signal. Verifies the pipeline's audit + pseudobulk + tumor-level
inference recover the planted effect and respect the tumor-as-unit design.
(inferCNV step is exercised separately in a resourced env; here CIN is injected so the
downstream statistics can be validated deterministically.)
"""
import numpy as np, pandas as pd, anndata, scanpy as sc, sys, os
sys.path.insert(0, os.path.dirname(__file__))
import sc_lag3_pipeline as P

rng = np.random.RandomState(0)  # deterministic
N_TUMORS = 24
genes = list(dict.fromkeys(P.CD8_MARKERS + P.CHECKPOINTS + P.EXHAUSTION +
                           P.ACTIVATION + P.STEMNESS + [f"BG{i}" for i in range(40)]))
gi = {g: j for j, g in enumerate(genes)}

# planted per-tumor CIN and a POSITIVE CIN->LAG3 effect in CD8 cells (effect size 0.8)
tumor_cin = rng.uniform(0, 3, N_TUMORS)
BETA_LAG3 = 0.8

X, obs = [], []
for t in range(N_TUMORS):
    plat = ["10x3p", "10x5p"][t % 2]
    drv  = ["none", "KRAS", "EGFR"][t % 3]          # "none" = driver-negative
    for kind, n in [("Epithelial cell (malignant)", 200), ("T cell CD8", 60), ("normal_ref", 80)]:
        for _ in range(n):
            v = rng.normal(0, 0.5, len(genes))
            if kind == "T cell CD8":
                for m in P.CD8_MARKERS: v[gi[m]] += 3
                v[gi["LAG3"]] += BETA_LAG3 * tumor_cin[t] + rng.normal(0, 0.3)  # planted signal
                v[gi["TIGIT"]] += rng.normal(0, 0.3)                             # no CIN effect (specificity)
                v[gi["IFNG"]] += 0.3 * tumor_cin[t]
            X.append(np.clip(v, 0, None))
            obs.append(dict(sample=f"T{t:02d}", patient=f"P{t:02d}", cell_type=kind,
                            platform=plat, condition="LUAD", treatment="naive", driver=drv))
ad = anndata.AnnData(np.array(X, dtype=np.float32),
                     obs=pd.DataFrame(obs), var=pd.DataFrame(index=genes))
ad.obs_names = [f"cell{i}" for i in range(ad.n_obs)]
print(f"synthetic AnnData: {ad.shape}, tumors={N_TUMORS}")

C = dict(sample="sample", patient="patient", celltype="cell_type",
         malignant=["Epithelial cell (malignant)"], cd8=["T cell CD8"],
         reference=["normal_ref"], histology="condition", treatment="treatment",
         platform="platform", dataset=None, driver="driver")

# 1) audit
tab, summary = P.audit(ad, C, ["LUAD"], ["naive"])
print("AUDIT:", summary)
assert summary["n_eligible_tumors"] == N_TUMORS, "audit should find all tumors eligible"

# 2) pseudobulk + inject known CIN + tumor-level test
pb = P.cd8_pseudobulk(ad, C)
cin = pd.Series(tumor_cin, index=[f"T{t:02d}" for t in range(N_TUMORS)], name="tumor_CIN")
merged = pb.join(cin, how="inner").join(tab["platform"], how="left")
merged["driver"] = tab["driver"]; merged["driver_negative"] = merged["driver"] == "none"
res = P.tumor_level_tests(merged, has_driver=True)
print("TUMOR-LEVEL RESULTS:", res)

# 3) assertions: recover planted positive LAG3 effect, TIGIT null, n=tumors
assert res["n_tumors"] == N_TUMORS
assert res["LAG3"]["coef"] > 0 and res["LAG3"]["p"] < 0.05, "should recover planted LAG3 signal"
assert res["TIGIT"]["p"] > 0.05, "TIGIT should be null (specificity)"
print("\nSMOKE TEST PASSED: pipeline recovers planted tumor-level CIN->LAG3 (coef=%.2f, p=%.1e), "
      "TIGIT null (p=%.2f), n_tumors=%d (unit=tumor)."
      % (res["LAG3"]["coef"], res["LAG3"]["p"], res["TIGIT"]["p"], res["n_tumors"]))
