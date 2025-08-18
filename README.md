# 🧬 10x Visium HD (8 µm) — End-to-End Analysis in Python

Welcome! 👋  
This guide walks you through an **analysis pipeline** for 10x Visium HD (8 µm bin) data using Python, [Scanpy](https://scanpy.readthedocs.io/), and [Squidpy](https://squidpy.readthedocs.io/).  

- ✅ Import Visium HD data directly from 10x outputs  
- ✅ Perform QC, filtering, normalization, HVG selection  
- ✅ PCA, UMAP, Leiden clustering  
- ✅ Spatial neighborhood statistics (Squidpy)  
- ✅ Plot genes + clusters on tissue  
- ✅ Rank marker genes, visualize expression panels  
- ✅ Annotate clusters using **literature-derived marker modules**  
- ✅ Export results and generate a **Methods text snippet**  


> **Tip for beginners**: Each section starts with an explanation of *why* we’re doing it, then gives runnable Python code.

---

## 📂 Data Structure

After processing with the 10x Genomics pipeline, your folder should look like this:

```

Visium\_HD\_data/
├── binned\_outputs/
│   └── square\_008um/
│       ├── filtered\_feature\_bc\_matrix.h5
│       └── spatial/
│           ├── tissue\_positions.parquet
│           ├── scalefactors\_json.json
│           ├── tissue\_hires\_image.png
│           └── tissue\_lowres\_image.png

````

- **`filtered_feature_bc_matrix.h5`** → expression matrix (spots × genes)  
- **`tissue_positions.parquet`** → where each spot is located on the tissue  
- **`scalefactors_json.json`** → scaling factors for plotting on histology images  
- **`tissue_hires_image.png` / `tissue_lowres_image.png`** → high/low-resolution histology images  

---

## 🚀 Getting Started

First, install the main Python packages (run this only once):

```python
# !pip install scanpy squidpy anndata h5py pandas numpy matplotlib pillow pyarrow leidenalg
````

Now import them and set paths to your data:

```python
import os, re, json
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt
from PIL import Image

RUN_DIR = "Visium_HD_data"   # change to your folder
BIN = "008"                  # '008' = 8 µm bin
BIN_DIR = os.path.join(RUN_DIR, f"binned_outputs/square_{BIN}um")
SPAT_DIR = os.path.join(BIN_DIR, "spatial")
MATRIX_H5 = os.path.join(BIN_DIR, "filtered_feature_bc_matrix.h5")
lib_id = f"square_{BIN}um"
```

---

## 📊 Loading the Data

Let’s load the raw counts into an `AnnData` object (Scanpy’s main data structure).
Think of `AnnData` like a “container” that holds your expression matrix, spot metadata, and spatial info.

```python
adata = sc.read_10x_h5(MATRIX_H5)
adata.var_names_make_unique()
print(adata)
```

---

## 🖼️ Adding Spatial Coordinates & Histology Images

Visium HD is special because every spot has **spatial coordinates**.
Here we connect the expression data to the histology images, so we can plot results directly on the tissue.

```python
pos = pd.read_parquet(os.path.join(SPAT_DIR, "tissue_positions.parquet"))

if "barcode" not in pos.columns:
    pos = pos.rename(columns={pos.columns[0]: "barcode"})
pos = pos.set_index("barcode")

_base = lambda s: s.to_series().str.replace(r"-\d+$", "", regex=True)
if not pos.index.isin(adata.obs_names).all():
    pos.index = _base(pos.index)
    adata.obs["__base"] = _base(adata.obs_names)
    pos = pos.reindex(adata.obs["__base"])
else:
    pos = pos.reindex(adata.obs_names)

adata.obsm["spatial"] = np.c_[pos["pxl_col_in_fullres"], pos["pxl_row_in_fullres"]]

with open(os.path.join(SPAT_DIR, "scalefactors_json.json")) as fh:
    scales = json.load(fh)
img_hires = np.array(Image.open(os.path.join(SPAT_DIR, "tissue_hires_image.png")))
img_low   = np.array(Image.open(os.path.join(SPAT_DIR, "tissue_lowres_image.png")))

adata.uns["spatial"] = {
    lib_id: {"images": {"hires": img_hires, "lowres": img_low}, "scalefactors": scales, "metadata": {}},
    "library_id": [lib_id],
}
```

---

## 🧹 Quality Control (QC)

Every dataset contains noisy spots or bad reads. We tag **mitochondrial** and **ribosomal** genes,
calculate QC metrics, and filter low-quality spots using **adaptive thresholds**.

```python
adata.var["mt"]   = adata.var_names.str.startswith(("mt-","mt."))
adata.var["ribo"] = adata.var_names.str.match(r"^(Rps|Rpl)\d", na=False)

sc.pp.calculate_qc_metrics(adata, qc_vars=["mt","ribo"], inplace=True)

sc.pl.violin(adata, ["n_genes_by_counts","total_counts","pct_counts_mt","pct_counts_ribo"],
             jitter=0.3, multi_panel=True)
```

Filter adaptively:

```python
low, high = adata.obs['total_counts'].quantile([0.02, 0.99])
min_counts = max(1000, int(low))
max_counts = max(int(high), 15000)

sc.pp.filter_cells(adata, min_counts=min_counts)
sc.pp.filter_cells(adata, max_counts=max_counts)
adata = adata[adata.obs['pct_counts_mt'] < 25].copy()
adata = adata[adata.obs['pct_counts_ribo'] < 30].copy()
sc.pp.filter_genes(adata, min_counts=10)
```

---

## ⚖️ Normalization, HVGs, PCA & UMAP

Now we normalize (to library size), log-transform, find **highly variable genes (HVGs)**,
and reduce dimensions using **PCA** and **UMAP**. Finally, we cluster spots with Leiden.

```python
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat_v3")
adata = adata[:, adata.var.highly_variable].copy()

sc.pp.pca(adata, n_comps=50)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata)

import leidenalg
sc.tl.leiden(adata, resolution=0.8, key_added="leiden_bin")
adata.obs["cluster"] = adata.obs["leiden_bin"].astype("category")

sc.pl.umap(adata, color=["cluster","total_counts","pct_counts_mt"], wspace=0.4)
sc.pl.spatial(adata, color="cluster", library_id=lib_id, spot_size=1.2)
```

---

## 🌐 Spatial Neighborhood Analysis

With Squidpy, we can test whether certain clusters prefer being neighbors or avoid each other.
This reveals spatial structure beyond just expression patterns.

```python
sq.gr.spatial_neighbors(adata, coord_type="generic", n_neighs=8)
sq.gr.nhood_enrichment(adata, cluster_key="cluster")
sq.pl.nhood_enrichment(adata, cluster_key="cluster")
```

---

## 🔬 Marker Genes

We identify marker genes for each cluster and visualize them with dot plots and heatmaps.

```python
sc.tl.rank_genes_groups(adata, "cluster", method="t-test")
sc.pl.rank_genes_groups(adata, n_genes=10, sharey=False)

rg_df = sc.get.rank_genes_groups_df(adata, group=None).sort_values(["group","pvals_adj","scores"])
top_df = rg_df.groupby("group").head(6)
ordered_genes = top_df["names"].unique().tolist()

sc.pl.dotplot(adata, var_names=ordered_genes, groupby="cluster", standard_scale="var")
```

---

## 🧭 Cluster Annotation

Finally, we assign biological meaning to clusters using curated marker sets
from the **Allen Brain Atlas** and the **Linnarson lab Mouse Brain Atlas**.

```python
marker_sets = {
    "Excitatory neurons": ["Slc17a7","Tbr1","Cux1","Rorb","Foxp2","Reln"],
    "Inhibitory neurons": ["Gad1","Gad2","Pvalb","Sst","Vip"],
    "Oligodendrocytes":  ["Mbp","Plp1","Mog","Cnp"],
    "Astrocytes":        ["Aqp4","Gfap","Aldh1l1","Slc1a3"],
    "Microglia":         ["C1qa","Cx3cr1","P2ry12"],
    "Endothelial":       ["Kdr","Pecam1","Cldn5"],
    "Hippocampal":       ["Prox1","Itpka"],
}

for name, genes in marker_sets.items():
    present = [g for g in genes if g in adata.var_names]
    if present:
        sc.tl.score_genes(adata, present, score_name=f"score_{name}")

avg_scores = adata.obs.groupby("cluster")[[c for c in adata.obs if c.startswith("score_")]].mean()
best_class = avg_scores.idxmax(axis=1).str.replace("score_", "", regex=False)
mapping = dict(zip(best_class.index.astype(str), best_class.values))
adata.obs["cluster_annot"] = adata.obs["cluster"].map(best_class)

print("Cluster → annotation:\n", mapping)
sc.pl.umap(adata, color=["cluster","cluster_annot"], wspace=0.4)
sc.pl.spatial(adata, color="cluster_annot", library_id=lib_id, spot_size=1.2)
```

---

## 💾 Save Your Work

```python
OUT = os.path.join(RUN_DIR, f"visium_hd_{BIN}um_processed.h5ad")
adata.write(OUT, compression="gzip")
print("Saved:", OUT)
```

---

## 📑 Methods (Auto-generated)

Copy-paste this into your manuscript/report:

```
Preprocessing followed the Scanpy spatial tutorial workflow:
https://scanpy-tutorials.readthedocs.io/en/latest/spatial/basic-analysis.html
We computed QC metrics (including mitochondrial and ribosomal fractions),
performed adaptive filtering, library-size normalization, log1p transform,
highly variable gene selection, PCA, neighborhood graph construction,
UMAP, and Leiden clustering.

Cluster annotation was guided by:
- Allen Brain Atlas (https://mouse.brain-map.org)
- Mouse Brain Gene Expression Atlas (Linnarson lab, http://mousebrain.org/)
- A recent preprint on spatial transcriptomics (https://www.biorxiv.org/content/10.1101/2020.07.24.219758v1)
```

---

## 📚 References

* [Scanpy spatial tutorial](https://scanpy-tutorials.readthedocs.io/en/latest/spatial/basic-analysis.html)
* [Allen Brain Atlas](https://mouse.brain-map.org/)
* [Mouse Brain Atlas — Linnarson Lab](http://mousebrain.org/)
* [Spatial transcriptomics preprint](https://www.biorxiv.org/content/10.1101/2020.07.24.219758v1)
* Wolf et al., *Genome Biology*, 2018 (Scanpy)
* Palla et al., *Nature Methods*, 2022 (Squidpy)

