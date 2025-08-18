# 10x Visium HD (Square 8 µm) — End-to-End Python Analysis

This is a **step-by-step walkthrough** of a complete analysis pipeline for 10x **Visium HD** binned data in Python.  
The guide mirrors the full analysis script, but breaks it into logical blocks with beginner-friendly explanations so you can understand **what’s happening at each step** and re-run it confidently.

> **What you’ll get**
>
> * Quality control (QC): metrics, filtering, normalization, and highly variable gene (HVG) selection  
> * Dimensionality reduction (PCA/UMAP) and clustering (Leiden)  
> * Spatial neighbors + neighborhood enrichment  
> * Marker discovery and summary plots  
> * A saved `.h5ad` file ready for downstream work  

---

## 📂 Data File Structure

After running the 10x Genomics pipeline for Visium HD, your data folder should look something like this:

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

````

**What each file is:**
- `filtered_feature_bc_matrix.h5`: the main **gene expression matrix** (spots × genes).  
- `tissue_positions.parquet`: coordinates telling you where each spot/barcode lies on the tissue image.  
- `scalefactors_json.json`: scaling factors that map spot coordinates to image resolution.  
- `tissue_hires_image.png` and `tissue_lowres_image.png`: the tissue histology images at high and low resolution.  

We’ll load all of these into Scanpy/Squidpy so that we can analyze expression data *and* visualize it in tissue space.

---

## 0) Environment & Install

> Run this once to install the required Python packages. It’s best to use a fresh conda environment with Python ≥3.9.

```python
# install (first time):
# !pip install scanpy squidpy anndata h5py pandas numpy matplotlib pillow pyarrow leidenalg
````

---

## 1) Imports & Basic Setup

> Import the libraries you’ll need:
>
> * **Scanpy** for core single-cell/spatial analysis
> * **Squidpy** for spatial graphs and enrichment
> * **Pandas/Numpy** for data wrangling
> * **Matplotlib/PIL** for visualization and image handling

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

## 2) Paths for Your Run

> Tell Python where your data lives.
>
> * `RUN_DIR` = your experiment folder
> * `BIN` = the bin size (e.g. `"008"` means 8 µm spots)
> * Other paths are built automatically.

```python
# paths
RUN_DIR = "Visium_HD_data"   # <-- change to your run folder
BIN = "008"                  # '008' for 8 µm
BIN_DIR = os.path.join(RUN_DIR, f"binned_outputs/square_{BIN}um")
SPAT_DIR = os.path.join(BIN_DIR, "spatial")
MATRIX_H5 = os.path.join(BIN_DIR, "filtered_feature_bc_matrix.h5")
lib_id = f"square_{BIN}um"
```

---

## 3) Load the Count Matrix

> Load the gene expression counts into an AnnData object.
> Make gene names unique (important if multiple genes share the same name).

```python
adata = sc.read_10x_h5(MATRIX_H5)
adata.var_names_make_unique()
print(adata)
```

---

## 4) Attach Spatial Coordinates & Images

> Add the pixel coordinates for each spot/barcode, and attach the histology images + scaling factors.
> This allows plotting clusters and gene expression **directly on the tissue image**.

```python
# read positions
pos = pd.read_parquet(os.path.join(SPAT_DIR, "tissue_positions.parquet"))

if "barcode" not in pos.columns:
    pos = pos.rename(columns={pos.columns[0]: "barcode"})
pos = pos.set_index("barcode")

# align to adata.obs_names
_base = lambda s: s.to_series().str.replace(r"-\d+$", "", regex=True)
if not pos.index.isin(adata.obs_names).all():
    pos.index = _base(pos.index)
    adata.obs["__base"] = _base(adata.obs_names)
    pos = pos.reindex(adata.obs["__base"])
else:
    pos = pos.reindex(adata.obs_names)

# pixel coords
pxr, pxc = "pxl_row_in_fullres", "pxl_col_in_fullres"
adata.obsm["spatial"] = np.c_[pos[pxc].to_numpy(), pos[pxr].to_numpy()]  # [x=col, y=row]

# scalefactors + images
with open(os.path.join(SPAT_DIR, "scalefactors_json.json"), "r") as fh:
    scales = json.load(fh)
img_hires = np.array(Image.open(os.path.join(SPAT_DIR, "tissue_hires_image.png")))
img_low   = np.array(Image.open(os.path.join(SPAT_DIR, "tissue_lowres_image.png")))

adata.uns["spatial"] = {
    lib_id: {"images": {"hires": img_hires, "lowres": img_low}, "scalefactors": scales, "metadata": {}},
    "library_id": [lib_id],
}
```

---

## 5) (Optional) Parse Grid Coordinates from Barcodes

> Some HD barcodes encode row/column directly. If so, we can extract that to build grid-based neighbors later.

```python
coords = adata.obs_names.to_series().str.extract(r"s_(\d{3})um_(\d+)_(\d+)(?:-\d+)?$")
coords.columns = ["um","array_row","array_col"]
if coords.notna().all().all():
    adata.obs[["array_row","array_col"]] = coords[["array_row","array_col"]].astype(float).values
```

---

## 6) QC Features (Mitochondrial & Ribosomal)

> Add simple QC annotations:
>
> * **Mitochondrial genes** (high % can mean stressed or dying cells)
> * **Ribosomal genes** (very high % can be problematic)

```python
adata.var["mt"]   = adata.var_names.str.startswith(("mt-","mt."))
adata.var["ribo"] = adata.var_names.str.match(r"^(Rps|Rpl)\d", na=False)

sc.pp.calculate_qc_metrics(adata, qc_vars=["mt","ribo"], inplace=True)
print("spots:", adata.n_obs, "| genes:", adata.n_vars)
print(adata.obs[["total_counts","n_genes_by_counts","pct_counts_mt","pct_counts_ribo"]].describe().T)
```

---

## 7) Initial QC Visualization

> Quick violin plots of QC metrics help you see if there are obvious outliers.

```python
qc_keys = [k for k in ["n_genes_by_counts","total_counts","pct_counts_mt","pct_counts_ribo"] if k in adata.obs.columns]
if qc_keys:
    sc.pl.violin(adata, qc_keys, jitter=0.3, multi_panel=True)
```

---

## 8) Data-Adaptive Filtering

> Filter out low-quality spots:
>
> * Use percentiles for total counts (to remove very empty or overloaded spots)
> * Cap mitochondrial/ribosomal percentages
> * Remove genes expressed in very few spots

```python
low, high = adata.obs['total_counts'].quantile([0.02, 0.99])
min_counts = max(200, int(low))
max_counts = max(int(high), 15000)
print(f"Using min_counts={min_counts}, max_counts={max_counts}")

sc.pp.filter_cells(adata, min_counts=min_counts)
sc.pp.filter_cells(adata, max_counts=max_counts)

adata = adata[adata.obs['pct_counts_mt'] < 25].copy()
adata = adata[adata.obs['pct_counts_ribo'] < 30].copy()

sc.pp.filter_genes(adata, min_counts=10)
```

---

## 9) Recompute QC After Subsetting

> Always re-check QC after filtering to see what remains.

```python
adata.var["mt"]   = adata.var_names.str.startswith(("mt-","mt."))
adata.var["ribo"] = adata.var_names.str.match(r"^(Rps|Rpl)\d", na=False)
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt","ribo"], inplace=True)

print("Remaining spots:", adata.n_obs, "| genes:", adata.n_vars)
qc_keys = [k for k in ["n_genes_by_counts","total_counts","pct_counts_mt","pct_counts_ribo"] if k in adata.obs.columns]
if adata.n_obs > 0 and qc_keys:
    sc.pl.violin(adata, qc_keys, jitter=0.3, multi_panel=True)
```

---

## 10) Normalize & Log-Transform

> Normalize each spot to 10,000 counts, then log-transform to stabilize variance.

```python
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
```

---

## 11) HVG Selection, PCA, Neighbors, UMAP

> * Select highly variable genes (HVGs)
> * Run PCA (dimensionality reduction)
> * Build a neighbor graph
> * Compute UMAP embedding for visualization

```python
sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat_v3")
adata = adata[:, adata.var.highly_variable].copy()

sc.pp.pca(adata, n_comps=50)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata)
```

⚠️ **Note:** This pipeline does **not scale** the data before PCA.
That means PCA emphasizes highly variable genes. If you want each gene to contribute equally (like in Seurat), add:

```python
sc.pp.scale(adata, max_value=10)
sc.pp.pca(adata, n_comps=50)
```

---

## 12) Leiden Clustering & UMAP

> Cluster spots using the Leiden algorithm, then visualize them on UMAP.

```python
import leidenalg
sc.tl.leiden(adata, resolution=0.8, key_added="leiden_bin")
adata.obs["cluster"] = adata.obs["leiden_bin"].astype("category")

sc.pl.umap(adata, color=["cluster","total_counts","pct_counts_mt"], wspace=0.4)
```

---

## 13) Spatial Neighbors & Neighborhood Enrichment

> Build a spatial neighbor graph and test whether clusters are enriched next to each other.

```python
sq.gr.spatial_neighbors(adata)
sq.gr.nhood_enrichment(adata, cluster_key="cluster")
sq.pl.nhood_enrichment(adata, cluster_key="cluster")
```

---

## 14) Spatial Gene & Cluster Maps

> Plot gene expression and clusters directly on the tissue image.

```python
for g in ["Olfm1", "Plp1", "Mbp"]:
    print(g, g in adata.var_names)

sq.pl.spatial_scatter(
    adata,
    color=["Olfm1", "Plp1", "Mbp", "cluster"],
    library_id="square_008um",   # match your bin size
    size=2
)
```

---

## 15) Marker Discovery & Plots

> Identify cluster marker genes and plot them in summary figures.

```python
sc.tl.rank_genes_groups(adata, "cluster", method="t-test")

sc.pl.rank_genes_groups(adata, n_genes=10, sharey=False)
sc.pl.rank_genes_groups_heatmap(adata, n_genes=10, groupby="cluster", show_gene_labels=True)

rg_df = sc.get.rank_genes_groups_df(adata, group=None) \
          .sort_values(["group","pvals_adj","scores"], ascending=[True, True, False])
topN = 6
top_df = rg_df.groupby("group", as_index=False, sort=False).head(topN)
ordered_genes = top_df.groupby("group")["names"].apply(list).explode().drop_duplicates().tolist()

sc.pl.dotplot(adata, var_names=ordered_genes, groupby="cluster", standard_scale="var", dendrogram=False)
sc.pl.heatmap(adata, var_names=ordered_genes, groupby="cluster",
              swap_axes=True, vmin=-2, vmax=2, cmap="viridis", show_gene_labels=True)
```

---

## 16) Save Processed Object

> Save your processed AnnData object so you don’t have to repeat everything next time.

```python
OUT = os.path.join(RUN_DIR, f"visium_hd_{BIN}um_processed.h5ad")
adata.write(OUT, compression="gzip")
print("Saved:", OUT)
```

---

## Notes & Tips

* **Species gene prefixes**: Mouse = `mt-` or `mt.`, human = `MT-`.
* **Scaling before PCA**: Not included by default (PCA emphasizes high-variance genes). Add `sc.pp.scale` if you want all genes weighted equally.
* **QC thresholds**: The script uses percentiles; adjust to your dataset.
* **Grid neighbors**: For HD, `coord_type="grid"` may better capture adjacency.
* **Plotting**: If running on a cluster or headless server, use `plt.show()` to render plots.

---

### Citations

* **Scanpy**: Wolf et al., *Genome Biology*, 2018
* **Squidpy**: Palla et al., *Nature Methods*, 2022
* **10x Genomics Visium HD**: [official documentation](https://www.10xgenomics.com/)
* **Example Data can be downloaded from**: [Visium_HD_3'_mouse_brain_data](https://www.10xgenomics.com/datasets/visium-hd-three-prime-mouse-brain-fresh-frozen)

