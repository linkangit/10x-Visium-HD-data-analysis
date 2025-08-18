# 🧬 10x Visium HD (Square 8 µm) — End-to-End Python Analysis

A **beginner-friendly walkthrough** of a complete analysis pipeline for 10x **Visium HD** spatial transcriptomics data in Python.  
This guide mirrors the full analysis script, but breaks it into clear blocks with explanations so you can understand **what’s happening at each step** and run it confidently.

---

## 🎯 What You’ll Do

- ✅ Perform quality control (QC): filtering, normalization, HVG selection  
- ✅ Run PCA/UMAP and cluster with Leiden  
- ✅ Build spatial neighbors + run neighborhood enrichment  
- ✅ Discover cluster marker genes  
- ✅ Save a processed `.h5ad` object ready for downstream analysis  

---

## 📂 Data File Structure

After running the 10x Genomics pipeline for Visium HD, your folder should look like this:

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

**What these files mean:**
- 📊 `filtered_feature_bc_matrix.h5`: expression counts (spots × genes)  
- 📐 `tissue_positions.parquet`: spatial coordinates of each spot/barcode  
- 📏 `scalefactors_json.json`: image scaling info  
- 🖼️ `tissue_hires_image.png` / `tissue_lowres_image.png`: histology images  

---

## ⚙️ Environment & Install

> 📝 Best practice: create a fresh conda environment with Python ≥ 3.9

```bash
# install (first time):
pip install scanpy squidpy anndata h5py pandas numpy matplotlib pillow pyarrow leidenalg
````

---

## 📥 Imports & Setup

```python
import os, re, json
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt
from PIL import Image
```

---

## 📁 Paths

```python
RUN_DIR = "Visium_HD_data"   # <-- change this to your run folder
BIN = "008"                  # '008' = 8 µm bins
BIN_DIR = os.path.join(RUN_DIR, f"binned_outputs/square_{BIN}um")
SPAT_DIR = os.path.join(BIN_DIR, "spatial")
MATRIX_H5 = os.path.join(BIN_DIR, "filtered_feature_bc_matrix.h5")
lib_id = f"square_{BIN}um"
```

---

## 📊 Load the Count Matrix

```python
adata = sc.read_10x_h5(MATRIX_H5)
adata.var_names_make_unique()
print(adata)
```

---

## 🖼️ Attach Spatial Coordinates & Images

```python
# read positions
pos = pd.read_parquet(os.path.join(SPAT_DIR, "tissue_positions.parquet"))
if "barcode" not in pos.columns:
    pos = pos.rename(columns={pos.columns[0]: "barcode"})
pos = pos.set_index("barcode")

# align
_base = lambda s: s.to_series().str.replace(r"-\d+$", "", regex=True)
if not pos.index.isin(adata.obs_names).all():
    pos.index = _base(pos.index)
    adata.obs["__base"] = _base(adata.obs_names)
    pos = pos.reindex(adata.obs["__base"])
else:
    pos = pos.reindex(adata.obs_names)

# pixel coords
adata.obsm["spatial"] = np.c_[pos["pxl_col_in_fullres"], pos["pxl_row_in_fullres"]]

# images + scalefactors
with open(os.path.join(SPAT_DIR, "scalefactors_json.json")) as fh:
    scales = json.load(fh)
img_hires = np.array(Image.open(os.path.join(SPAT_DIR, "tissue_hires_image.png")))
img_low   = np.array(Image.open(os.path.join(SPAT_DIR, "tissue_lowres_image.png")))

adata.uns["spatial"] = {
    lib_id: {"images": {"hires": img_hires, "lowres": img_low}, "scalefactors": scales, "metadata": {}},
    "library_id": [lib_id],
}
```

> 💡 **Tip:** Scanpy expects `[x=column, y=row]` order for spatial coordinates.

---

## 🧩 (Optional) Parse Grid Coordinates

```python
coords = adata.obs_names.to_series().str.extract(r"s_(\d{3})um_(\d+)_(\d+)(?:-\d+)?$")
coords.columns = ["um","array_row","array_col"]
if coords.notna().all().all():
    adata.obs[["array_row","array_col"]] = coords[["array_row","array_col"]].astype(float).values
```

---

## 🔍 QC Features (Mito & Ribosomal)

```python
adata.var["mt"]   = adata.var_names.str.startswith(("mt-","mt."))
adata.var["ribo"] = adata.var_names.str.match(r"^(Rps|Rpl)\d", na=False)
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt","ribo"], inplace=True)
```

---

## 📈 Initial QC Visualization

```python
qc_keys = [k for k in ["n_genes_by_counts","total_counts","pct_counts_mt","pct_counts_ribo"] if k in adata.obs.columns]
if qc_keys:
    sc.pl.violin(adata, qc_keys, jitter=0.3, multi_panel=True)
```

---

## 🧹 Data-Adaptive Filtering

```python
low, high = adata.obs['total_counts'].quantile([0.02, 0.99])
min_counts = max(200, int(low))
max_counts = max(int(high), 15000)

sc.pp.filter_cells(adata, min_counts=min_counts)
sc.pp.filter_cells(adata, max_counts=max_counts)

adata = adata[adata.obs['pct_counts_mt'] < 25].copy()
adata = adata[adata.obs['pct_counts_ribo'] < 30].copy()

sc.pp.filter_genes(adata, min_counts=10)
```

---

## 🔄 Recompute QC

```python
adata.var["mt"]   = adata.var_names.str.startswith(("mt-","mt."))
adata.var["ribo"] = adata.var_names.str.match(r"^(Rps|Rpl)\d", na=False)
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt","ribo"], inplace=True)
```

---

## ⚖️ Normalize & Log-Transform

```python
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
```

---

## 📉 HVG, PCA, Neighbors, UMAP

```python
sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat_v3")
adata = adata[:, adata.var.highly_variable].copy()

sc.pp.pca(adata, n_comps=50)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata)
```

> ⚠️ **Scaling note:** This script does *not* scale before PCA.
> Add `sc.pp.scale(adata, max_value=10)` if you want each gene to contribute equally.

---

## 🧭 Leiden Clustering & UMAP

```python
import leidenalg
sc.tl.leiden(adata, resolution=0.8, key_added="leiden_bin")
adata.obs["cluster"] = adata.obs["leiden_bin"].astype("category")
sc.pl.umap(adata, color=["cluster","total_counts","pct_counts_mt"], wspace=0.4)
```

---

## 🌐 Spatial Neighbors & Enrichment

```python
sq.gr.spatial_neighbors(adata)
sq.gr.nhood_enrichment(adata, cluster_key="cluster")
sq.pl.nhood_enrichment(adata, cluster_key="cluster")
```

---

## 🗺️ Spatial Gene & Cluster Maps

```python
for g in ["Olfm1", "Plp1", "Mbp"]:
    print(g, g in adata.var_names)

sq.pl.spatial_scatter(
    adata,
    color=["Olfm1", "Plp1", "Mbp", "cluster"],
    library_id="square_008um",
    size=2
)
```

---

## 🧪 Marker Discovery & Plots

```python
sc.tl.rank_genes_groups(adata, "cluster", method="t-test")

sc.pl.rank_genes_groups(adata, n_genes=10, sharey=False)
sc.pl.rank_genes_groups_heatmap(adata, n_genes=10, groupby="cluster", show_gene_labels=True)

rg_df = sc.get.rank_genes_groups_df(adata, group=None).sort_values(["group","pvals_adj","scores"], ascending=[True, True, False])
topN = 6
top_df = rg_df.groupby("group", as_index=False, sort=False).head(topN)
ordered_genes = top_df.groupby("group")["names"].apply(list).explode().drop_duplicates().tolist()

sc.pl.dotplot(adata, var_names=ordered_genes, groupby="cluster", standard_scale="var", dendrogram=False)
sc.pl.heatmap(adata, var_names=ordered_genes, groupby="cluster", swap_axes=True, vmin=-2, vmax=2, cmap="viridis", show_gene_labels=True)
```

---

## 💾 Save Processed Object

```python
OUT = os.path.join(RUN_DIR, f"visium_hd_{BIN}um_processed.h5ad")
adata.write(OUT, compression="gzip")
print("Saved:", OUT)
```

---

## 📝 Notes & Tips

* 🐭 **Species gene prefixes**: Mouse = `mt-` / `mt.`, human = `MT-`.
* ⚖️ **Scaling before PCA**: Optional — add `sc.pp.scale` if you want all genes weighted equally.
* 📏 **QC thresholds**: Percentile-based defaults here; tweak for your dataset.
* 🧩 **Grid neighbors**: Use `coord_type="grid"` for true lattice adjacency in HD.
* 🖥️ **Plotting**: If running on a server, call `plt.show()` to render plots.

---

## 📚 References

* **Scanpy**: Wolf et al., *Genome Biology*, 2018
* **Squidpy**: Palla et al., *Nature Methods*, 2022
* **10x Genomics Visium HD**: [official documentation](https://www.10xgenomics.com/)
* Example dataset: [Visium HD 3' mouse brain](https://www.10xgenomics.com/datasets/visium-hd-three-prime-mouse-brain-fresh-frozen)

```

---

✨ This version is:
- **Visually appealing** with emojis and clear headers  
- **Readable for beginners** (plain-English explanations, tips inline)  
- **GitHub-friendly** (code blocks, quotes, bullets)

Would you like me to also generate a **short “Quickstart” snippet** at the top (like ~15 lines of code for impatient users who don’t want explanations), so your README works both for beginners *and* advanced readers?
```
