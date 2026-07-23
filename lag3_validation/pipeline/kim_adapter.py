#!/usr/bin/env python3
"""
Kim GSE131907 -> AnnData adapter (run in a resourced environment).

Builds an h5ad the main pipeline can consume, from the GEO files:
  GSE131907_Lung_Cancer_normalized_log2TPM_matrix.txt.gz   (genes x cells, dense TSV)
  GSE131907_Lung_Cancer_cell_annotation.txt.gz

NOTE: the dense normalized matrix is ~4-5 GB and 208k cells; the GEO FTP endpoint
capped/stalled the download in our ephemeral container (hard stall at ~2.9 GB), so this
adapter is intended for a resourced environment with a reliable connection. Alternatively
load the sparse `raw_UMI_matrix.rds.gz` in R/Seurat and export to h5ad (smaller, sparse).

Kim is a SMOKE-TEST / descriptive cohort only: ~11-15 LUAD tumors (tLung/tL/B), so the
independent-unit count is ~13 and any CIN->LAG3 correlation is EXPLORATORY. Driver
annotation is not per-cell in GEO; treat driver-negative as untested here.
"""
import gzip, numpy as np, pandas as pd, anndata, scanpy as sc, argparse

# cell types worth carrying (streaming keeps only cells in tumor+normal-lung samples)
CD8_LABEL = "T cell CD8"          # assigned below from CD8A+ within 'T lymphocytes'
MALIGNANT = "Epithelial cell (malignant)"

def build(matrix_gz, annot_gz, out_h5ad, keep_origins=("tLung","tL/B","nLung")):
    ann = pd.read_csv(annot_gz, sep="\t")
    ann = ann[ann["Sample_Origin"].isin(keep_origins)].copy()
    keep_cells = set(ann["Index"])
    # stream dense matrix; keep only columns (cells) in keep_cells
    with gzip.open(matrix_gz, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        cell_ids = header[1:]
        col_keep = [i for i, c in enumerate(cell_ids) if c in keep_cells]
        kept_cells = [cell_ids[i] for i in col_keep]
        genes, rows = [], []
        for line in f:
            parts = line.rstrip("\n").split("\t")
            genes.append(parts[0])
            vals = np.array(parts[1:], dtype=np.float32)
            rows.append(vals[col_keep])
    X = np.vstack(rows).T            # cells x genes
    ad = anndata.AnnData(X, obs=pd.DataFrame(index=kept_cells), var=pd.DataFrame(index=genes))
    a = ann.set_index("Index").loc[kept_cells]
    ad.obs["sample"] = a["Sample"].values
    ad.obs["patient"] = a["Sample"].values           # 1 sample/patient in Kim tumors
    ad.obs["origin"] = a["Sample_Origin"].values
    ad.obs["condition"] = "LUAD"
    ad.obs["treatment"] = "naive"
    ad.obs["platform"] = "10x3p"
    # cell types: malignant = tumor-origin Epithelial; CD8 = T lymphocytes w/ CD8A>0
    ct = a["Cell_type"].astype(str).values
    ad.obs["cell_type"] = ct
    is_epi = np.char.startswith(ct.astype(str), "Epithelial")
    is_tum = np.isin(ad.obs["origin"].values, ["tLung","tL/B"])
    ad.obs.loc[is_epi & is_tum, "cell_type"] = MALIGNANT
    is_t = ct == "T lymphocytes"
    cd8a = np.asarray(ad[:, "CD8A"].X).ravel() if "CD8A" in ad.var_names else np.zeros(ad.n_obs)
    ad.obs.loc[is_t & (cd8a > 0), "cell_type"] = CD8_LABEL
    ad.write_h5ad(out_h5ad)
    print(f"wrote {out_h5ad}: {ad.shape}; malignant={int((ad.obs.cell_type==MALIGNANT).sum())}, "
          f"CD8={int((ad.obs.cell_type==CD8_LABEL).sum())}, tumors={ad.obs[is_tum].sample.nunique()}")

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", required=True); ap.add_argument("--annot", required=True)
    ap.add_argument("--out", default="kim.h5ad")
    a = ap.parse_args(); build(a.matrix, a.annot, a.out)
