# %%
import scanpy as sc
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
adata.var["feature_types"].value_counts()