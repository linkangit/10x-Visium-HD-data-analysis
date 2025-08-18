Got it 👍 — let’s keep everything self-contained inside a **single `README.md`** for your GitHub repo.

Below is a **ready-to-commit `README.md`** with:

* Project description
* Table of contents
* Step-by-step pipeline with **foldable code blocks**
* Clear notes about scaling before PCA (since you asked about that)
* References

---

````markdown
# 10x Visium HD Analysis (Python)

This repository provides a **complete, end-to-end pipeline** for analyzing **10x Genomics Visium HD** spatial transcriptomics data with [Scanpy](https://scanpy.readthedocs.io/en/stable/) and [Squidpy](https://squidpy.readthedocs.io/).

The workflow covers:
- Loading binned outputs (`square_008um` etc.)
- Attaching pixel coordinates and histology images
- QC metrics and filtering
- Normalization, HVG selection
- PCA, UMAP, Leiden clustering
- Spatial neighbors & neighborhood enrichment
- Marker gene discovery and visualization
- Saving a processed `.h5ad` file

---

## 📖 Table of Contents
1. [Installation](#installation)
2. [Input File Structure](#input-file-structure)
3. [Step-by-Step Analysis](#step-by-step-analysis)
    - [Imports & Setup](#imports--setup)
    - [Load Count Matrix](#load-count-matrix)
    - [Attach Spatial Coordinates & Images](#attach-spatial-coordinates--images)
    - [Quality Control](#quality-control)
    - [Filtering](#filtering)
    - [Normalization & HVG Selection](#normalization--hvg-selection)
    - [Dimensionality Reduction](#dimensionality-reduction)
    - [Clustering & UMAP](#clustering--umap)
    - [Spatial Neighbors & Enrichment](#spatial-neighbors--enrichment)
    - [Marker Gene Discovery](#marker-gene-discovery)
4. [Saving Results](#saving-results)
5. [Notes & Tips](#notes--tips)
6. [References](#references)

---

## Installation

```bash
# Create a fresh environment
conda create -n visiumhd python=3.10
conda activate visiumhd

# Install required packages
pip install scanpy squidpy anndata h5py pandas numpy matplotlib pillow pyarrow leidenalg
````

---

## Input File Structure

Your project folder should look like:

```
Visium_HD_data/
├── binned_outputs/
│   └── square_008um/
│       ├── filtered_feature_bc_matrix.h5
│       └── spatial/
│           ├── tissue_positions.parquet
│           ├── scalefactors_json.json
│           ├── tissue_hires_image.png
│           └── tissue_lowres_image.png
```

---

## Step-by-Step Analysis

### Imports & Setup

<details>
<summary>Show code</summary>

```python
import os, re, json
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt
from PIL import Image

RUN_DIR = "Visium_HD_data"   # change to your run folder
BIN = "008"                  # '008' for 8 µm
BIN_DIR = os.path.join(RUN_DIR, f"binned_outputs/square_{BIN}um")
SPAT_DIR = os.path.join(BIN_DIR, "spatial")
MATRIX_H5 = os.path.join(BIN_DIR, "filtered_feature_bc_matrix.h5")
lib_id = f"square_{BIN}um"
```

</details>

---

### Load Count Matrix

<details>
<summary>Show code</summary>

```python
adata = sc.read_10x_h5(MATRIX_H5)
adata.var_names_make_unique()
print(adata)
```

</details>

---

### Attach Spatial Coordinates & Images

<details>
<summary>Show code</summary>

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

with open(os.path.join(SPAT_DIR, "scalefactors_json.json"), "r") as fh:
    scales = json.load(fh)

img_hires = np.array(Image.open(os.path.join(SPAT_DIR, "tissue_hires_image.png")))
img_low   = np.array(Image.open(os.path.join(SPAT_DIR, "tissue_lowres_image.png")))

adata.uns["spatial"] = {
    lib_id: {"images": {"hires": img_hires, "lowres": img_low}, "scalefactors": scales, "metadata": {}},
    "library_id": [lib_id],
}
```

</details>

---

### Quality Control

<details>
<summary>Show code</summary>

```python
adata.var["mt"]   = adata.var_names.str.startswith(("mt-","mt."))
adata.var["ribo"] = adata.var_names.str.match(r"^(Rps|Rpl)\d", na=False)

sc.pp.calculate_qc_metrics(adata, qc_vars=["mt","ribo"], inplace=True)
sc.pl.violin(adata, ["n_genes_by_counts","total_counts","pct_counts_mt","pct_counts_ribo"], jitter=0.3, multi_panel=True)
```

</details>

---

### Filtering

<details>
<summary>Show code</summary>

```python
low, high = adata.obs['total_counts'].quantile([0.02, 0.99])
sc.pp.filter_cells(adata, min_counts=max(200, int(low)))
sc.pp.filter_cells(adata, max_counts=max(int(high), 15000))

adata = adata[adata.obs['pct_counts_mt'] < 25]
adata = adata[adata.obs['pct_counts_ribo'] < 30]

sc.pp.filter_genes(adata, min_counts=10)
```

</details>

---

### Normalization & HVG Selection

<details>
<summary>Show code</summary>

```python
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat_v3")
adata = adata[:, adata.var.highly_variable].copy()
```

⚠️ *Optional*: Add `sc.pp.scale(adata, max_value=10)` before PCA if you want to scale genes to unit variance (not done by default here).

</details>

---

### Dimensionality Reduction

<details>
<summary>Show code</summary>

```python
sc.pp.pca(adata, n_comps=50)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata)
```

</details>

---

### Clustering & UMAP

<details>
<summary>Show code</summary>

```python
import leidenalg
sc.tl.leiden(adata, resolution=0.8, key_added="leiden_bin")
adata.obs["cluster"] = adata.obs["leiden_bin"].astype("category")

sc.pl.umap(adata, color=["cluster","total_counts","pct_counts_mt"], wspace=0.4)
```

</details>

---

### Spatial Neighbors & Enrichment

<details>
<summary>Show code</summary>

```python
sq.gr.spatial_neighbors(adata)
sq.gr.nhood_enrichment(adata, cluster_key="cluster")
sq.pl.nhood_enrichment(adata, cluster_key="cluster")
```

</details>

---

### Marker Gene Discovery

<details>
<summary>Show code</summary>

```python
sc.tl.rank_genes_groups(adata, "cluster", method="t-test")
sc.pl.rank_genes_groups(adata, n_genes=10, sharey=False)
sc.pl.rank_genes_groups_heatmap(adata, n_genes=10, groupby="cluster", show_gene_labels=True)
```

</details>

---

## Saving Results

```python
OUT = os.path.join(RUN_DIR, f"visium_hd_{BIN}um_processed.h5ad")
adata.write(OUT, compression="gzip")
print("Saved:", OUT)
```

---

## Notes & Tips

* **Scaling before PCA**: Not included by default. Add `sc.pp.scale(adata)` if you want each gene to contribute equally.
* **Species**: For human data, adjust mito prefix to `MT-`.
* **Grid neighbors**: For HD data, `sq.gr.spatial_neighbors(..., coord_type="grid")` may better reflect spot adjacency.
* **QC thresholds**: The script uses percentiles + soft mito/ribo cutoffs as a flexible starting point.
* **Images**: Ensure the library ID (`square_008um`) matches your folder name when plotting with Squidpy.

---

## References

* Wolf et al., *Genome Biology*, 2018 — [Scanpy](https://scanpy.readthedocs.io/)
* Palla et al., *Nature Methods*, 2022 — [Squidpy](https://squidpy.readthedocs.io/)
* 10x Genomics — [Visium HD documentation](https://www.10xgenomics.com/)

```

---

Would you like me to also add **example output plots** (UMAP, violin, spatial scatter) into the README so readers see what they should expect?
```
