# %%
import anndata2ri
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

# R interface
from rpy2.robjects import pandas2ri

pandas2ri.activate()
anndata2ri.activate()
# %%
adata_raw = sc.read_h5ad("/data/trm-sc-analysis/data/external/openproblems_bmmc_multiome_genes_filtered.h5ad")
adata_raw.layers["logcounts"] = adata_raw.X

label_key = "cell_type"
batch_key = "batch"
# %%
# They just do this to reduce computational time
keep_batches = ["s1d3", "s2d1", "s3d7"]
adata = adata_raw[adata_raw.obs[batch_key].isin(keep_batches)].copy()
# %%
# This is just getting the gene expression only
adata = adata[:, adata.var["feature_types"] == "GEX"].copy()
sc.pp.filter_genes(adata, min_cells=1)
# %%
# It's unclear if adata.layers['counts'] is un-normalized or if it just needs to be normalized
adata.X = adata.layers["counts"].copy()
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
adata.layers["logcounts"] = adata.X.copy()
# %%
sc.pp.highly_variable_genes(adata)
sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
# %%
adata.uns[batch_key + "_colors"] = [
    "#1b9e77",
    "#d95f02",
    "#7570b3",
]  # Set custom colours for batches
sc.pl.umap(adata, color=[label_key, batch_key], wspace=1)
# %%
sc.pp.highly_variable_genes(
    adata, n_top_genes=2000, flavor="cell_ranger", batch_key=batch_key
)
# %%
n_batches = adata.var["highly_variable_nbatches"].value_counts()
ax = n_batches.plot(kind="bar")
# %%
adata_hvg = adata[:, adata.var["highly_variable"]].copy()
# %% [markdown]
"""
# Seurat Method
"""
# %%
adata_seurat = adata_hvg.copy()
# Convert categorical columns to strings
adata_seurat.obs[batch_key] = adata_seurat.obs[batch_key].astype(str)
adata_seurat.obs[label_key] = adata_seurat.obs[label_key].astype(str)
# Delete uns as this can contain arbitrary objects which are difficult to convert
del adata_seurat.uns
adata_seurat
# %%
import anndata2ri
import rpy2.robjects as ro
from rpy2.robjects.packages import importr

importr('Seurat')
integrateSeurat = ro.functions.wrap_r_function(
    ro.r(
        """
        function(adata_seurat, batch_key){
        seurat <- as.Seurat(adata_seurat, counts = "counts", data = "logcounts")
        batch_list <- SplitObject(seurat, split.by = batch_key)
        anchors <- FindIntegrationAnchors(batch_list, anchor.features = rownames(seurat))
        integrated <- IntegrateData(anchors)
        # Extract the integrated expression matrix
        integrated_expr <- GetAssayData(integrated)
        # Make sure the rows and columns are in the same order as the original object
        integrated_expr <- integrated_expr[rownames(seurat), colnames(seurat)]
        # Transpose the matrix to AnnData format
        integrated_expr <- t(integrated_expr)
        print(integrated_expr[1:10, 1:10])
        }
        """
    ),
    "integrateSeurat",
)

with (
    ro.default_converter + ro.pandas2ri.converter + anndata2ri.converter
).context():
    outs = integrateSeurat(adata_seurat, batch_key)

# seurat <- as.Seurat(adata_seurat, counts = "counts", data = "logcounts")