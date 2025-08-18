Perfect ✅ Let’s do both options so you can choose what fits your repo best:

---

# **Option 1: `README.md` (ready to commit)**

Here’s a full GitHub-ready Markdown with:

* Project title + description
* Table of contents
* Walkthrough (with collapsible code blocks so readers aren’t overwhelmed)
* Notes & citation

````markdown
# 10x Visium HD Analysis (Python)

A complete **end-to-end pipeline** for processing and analyzing **10x Visium HD binned data** with [Scanpy](https://scanpy.readthedocs.io/en/stable/) and [Squidpy](https://squidpy.readthedocs.io/).

---

## 📖 Table of Contents
1. [Overview](#overview)
2. [Installation](#installation)
3. [File Structure](#file-structure)
4. [Step-by-Step Analysis](#step-by-step-analysis)
    - [Imports & Setup](#1-imports--setup)
    - [Load Data](#2-load-data)
    - [Attach Spatial Info](#3-attach-spatial-info)
    - [Quality Control](#4-quality-control)
    - [Filtering](#5-filtering)
    - [Normalization & HVG](#6-normalization--hvg)
    - [Dimensionality Reduction](#7-dimensionality-reduction)
    - [Clustering & UMAP](#8-clustering--umap)
    - [Spatial Neighbors & Enrichment](#9-spatial-neighbors--enrichment)
    - [Marker Discovery](#10-marker-discovery)
6. [Saving Results](#saving-results)
7. [Notes & Tips](#notes--tips)
8. [References](#references)

---

## Overview
This workflow:
- Loads **Visium HD binned outputs** (e.g., `square_008um`)
- Performs **QC, filtering, normalization, HVG selection**
- Computes **PCA/UMAP and Leiden clustering**
- Builds **spatial neighbors & neighborhood enrichment**
- Identifies **marker genes** and generates plots
- Saves a compact `.h5ad` object

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

## File Structure

Your folder should look like:

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

### 1. Imports & Setup

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
BIN = "008"                  # bin size: '008' for 8 µm
BIN_DIR = os.path.join(RUN_DIR, f"binned_outputs/square_{BIN}um")
SPAT_DIR = os.path.join(BIN_DIR, "spatial")
MATRIX_H5 = os.path.join(BIN_DIR, "filtered_feature_bc_matrix.h5")
lib_id = f"square_{BIN}um"
```

</details>

---

### 2. Load Data

<details>
<summary>Show code</summary>

```python
adata = sc.read_10x_h5(MATRIX_H5)
adata.var_names_make_unique()
print(adata)
```

</details>

---

### 3. Attach Spatial Info

<details>
<summary>Show code</summary>

```python
# read parquet, align barcodes, add pixel coords, images, scalefactors...
```

</details>

---

### 4. Quality Control

<details>
<summary>Show code</summary>

```python
# mark mito/ribo, calculate QC, violin plots
```

</details>

---

### 5. Filtering

<details>
<summary>Show code</summary>

```python
# percentile thresholds, mito/ribo caps, filter genes
```

</details>

---

### 6. Normalization & HVG

<details>
<summary>Show code</summary>

```python
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat_v3")
```

⚠️ *Optional*: Add `sc.pp.scale(adata, max_value=10)` before PCA if you prefer scaled input.

</details>

---

### 7. Dimensionality Reduction

<details>
<summary>Show code</summary>

```python
sc.pp.pca(adata, n_comps=50)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata)
```

</details>

---

### 8. Clustering & UMAP

<details>
<summary>Show code</summary>

```python
sc.tl.leiden(adata, resolution=0.8, key_added="leiden_bin")
adata.obs["cluster"] = adata.obs["leiden_bin"].astype("category")
sc.pl.umap(adata, color=["cluster","total_counts","pct_counts_mt"])
```

</details>

---

### 9. Spatial Neighbors & Enrichment

<details>
<summary>Show code</summary>

```python
sq.gr.spatial_neighbors(adata)
sq.gr.nhood_enrichment(adata, cluster_key="cluster")
sq.pl.nhood_enrichment(adata, cluster_key="cluster")
```

</details>

---

### 10. Marker Discovery

<details>
<summary>Show code</summary>

```python
sc.tl.rank_genes_groups(adata, "cluster", method="t-test")
sc.pl.rank_genes_groups(adata, n_genes=10, sharey=False)
```

</details>

---

## Saving Results

```python
OUT = os.path.join(RUN_DIR, f"visium_hd_{BIN}um_processed.h5ad")
adata.write(OUT, compression="gzip")
```

---

## Notes & Tips

* For **human data**, change mito gene prefix to `MT-`.
* Consider **grid neighbors** for HD (`coord_type="grid"`).
* Scaling before PCA is optional; include it if you want unit variance.

---

## References

* Wolf et al., *Genome Biology*, 2018 (Scanpy)
* Palla et al., *Nature Methods*, 2022 (Squidpy)
* 10x Genomics [Visium HD documentation](https://www.10xgenomics.com/)

````

---

# **Option 2: Split into `analysis.py` + `config.yaml`**

This makes the workflow cleaner:

### `config.yaml`
```yaml
# Config for Visium HD analysis
run_dir: "Visium_HD_data"
bin: "008"
min_counts: 200
max_counts: 15000
mito_cutoff: 25
ribo_cutoff: 30
n_hvgs: 3000
n_pcs: 50
n_neighbors: 15
umap_pcs: 30
resolution: 0.8
````

### `analysis.py`

```python
import os, json
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
from PIL import Image
import yaml

# Load config
with open("config.yaml") as f:
    cfg = yaml.safe_load(f)

RUN_DIR = cfg["run_dir"]
BIN = cfg["bin"]
BIN_DIR = os.path.join(RUN_DIR, f"binned_outputs/square_{BIN}um")
SPAT_DIR = os.path.join(BIN_DIR, "spatial")
MATRIX_H5 = os.path.join(BIN_DIR, "filtered_feature_bc_matrix.h5")
lib_id = f"square_{BIN}um"

# Load counts
adata = sc.read_10x_h5(MATRIX_H5)
adata.var_names_make_unique()

# ... [insert spatial, QC, filtering, normalization, PCA/UMAP, clustering as before, using cfg params]

OUT = os.path.join(RUN_DIR, f"visium_hd_{BIN}um_processed.h5ad")
adata.write(OUT, compression="gzip")
```
