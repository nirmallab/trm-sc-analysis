# %%
import anndata
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import scanpy as sc
from scipy.stats import median_abs_deviation
import seaborn as sns
# %%
scDir = Path('../../data/external/GSE277165')
adatas = []
usedPrefix = set()
c =0
for scPath in scDir.iterdir():
    prefix = '_'.join(scPath.name.split("_")[0:2])+"_"
    if '_ATAC_' in scPath.name or 'Nae' in scPath.name:
        continue
    if prefix in usedPrefix:
        continue
    print(prefix.split('_')[1].strip('_'))
    adata = sc.read_10x_mtx(path=scDir, prefix=prefix)
    adata.obs['batch'] = prefix.split('_')[1].strip('_')
    adatas.append(adata)
    usedPrefix.add(prefix)
    c += 1
adata = anndata.concat(adatas)
# %%
# mitochondrial genes
adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
# ribosomal genes
adata.var["ribo"] = adata.var_names.str.upper().str.startswith(("RPS", "RPL"))
# hemoglobin genes.
adata.var["hb"] = adata.var_names.str.upper().str.contains("^HB[^(P)]")
# %%
sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True
)
# %%
p1 = sns.displot(adata.obs["total_counts"], bins=100, kde=False)
# sc.pl.violin(adata, 'total_counts')
p2 = sc.pl.violin(adata, "pct_counts_mt")
p3 = sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")
# %%
def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier
# %%
adata.obs["outlier"] = (
    is_outlier(adata, "log1p_total_counts", 5)
    | is_outlier(adata, "log1p_n_genes_by_counts", 5)
    | is_outlier(adata, "pct_counts_in_top_20_genes", 5)
)
adata.obs.outlier.value_counts()
# %%
adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 3) | (
    adata.obs["pct_counts_mt"] > 8
)
adata.obs.mt_outlier.value_counts()
# %%
print(f"Total number of cells: {adata.n_obs}")
adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()

print(f"Number of cells after filtering of low quality cells: {adata.n_obs}")

adata.layers["counts"] = adata.X
# %%
import anndata2ri
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
# %%
data_mat = adata.X.T
# samples = adata.obs["sample"]
importr("SingleCellExperiment")
importr("scDblFinder")
importr("BiocParallel")

detect_doublets = ro.functions.wrap_r_function(
    ro.r(
        """
function(data_mat){
print("data loaded, proceeding with doublet detection")
# bp <- MulticoreParam(cores, RNGseed=9)
sce = scDblFinder(
    SingleCellExperiment(
        list(counts=data_mat),
    )
)
return(
    list(
        score=sce$scDblFinder.score,
        class=sce$scDblFinder.class
    )
)
}
"""
    ),
    "detect_doublets",
)
with (
    ro.default_converter + ro.pandas2ri.converter + anndata2ri.converter
).context():
    outs = detect_doublets(data_mat)

adata.obs["scDblFinder_score"] = outs["score"]
adata.obs["scDblFinder_class"] = outs["class"]

# Remove doublets based on downstream analysis
# adata = adata[adata.obs['scDblFinder_class'] == 'singlet']
# %%
scales_counts = sc.pp.normalize_total(adata, target_sum=None, inplace=False)
adata.layers["log1p_norm"] = sc.pp.log1p(scales_counts["X"], copy=True)
# %%
fig, axes = plt.subplots(1, 2, figsize=(10, 5))
p1 = sns.histplot(adata.obs["total_counts"], bins=100, kde=False, ax=axes[0])
axes[0].set_title("Total counts")
p2 = sns.histplot(adata.layers["log1p_norm"].sum(1), bins=100, kde=False, ax=axes[1])
axes[1].set_title("Shifted logarithm")
plt.show()
# %%
importr('scry')

featureSelection = ro.functions.wrap_r_function(
    ro.r(
        """
function(adata){
print("anndata loaded")
print(adata)
sce = devianceFeatureSelection(adata, assay="X")
return(rowData(sce)$binomial_deviance)
}
"""
    ),
    "featureSelection",
)
with (
    ro.default_converter + ro.pandas2ri.converter + anndata2ri.converter
).context():
    binomial_deviance = featureSelection(adata).T
# %%
idx = binomial_deviance.argsort()[-4000:]
mask = np.zeros(adata.var_names.shape, dtype=bool)
mask[idx] = True

adata.var["highly_deviant"] = mask
adata.var["binomial_deviance"] = binomial_deviance
# %%
sc.pp.highly_variable_genes(adata, layer="log1p_norm")
# %%
ax = sns.scatterplot(
    data=adata.var, x="means", y="dispersions", hue="highly_deviant", s=5
)
ax.set_xlim(None, 1.5)
ax.set_ylim(None, 3)
plt.show()
# %%
# %%
adata.X = adata.layers["log1p_norm"]
# %%
adata.var["highly_variable"] = adata.var["highly_deviant"]
sc.pp.pca(adata, svd_solver="arpack", use_highly_variable=True)
sc.pl.pca_scatter(adata, color="total_counts")
# %%
sc.pp.neighbors(adata)
sc.tl.umap(adata)

adata.write_h5ad('../../data/external/stubenvoll_processed.h5ad')
# %%
adata = sc.read_h5ad('../../data/external/stubenvoll_processed.h5ad')
sc.pl.umap(adata, color="batch")
# %%
sc.pl.umap(
    adata,
    color=["total_counts", "pct_counts_mt", "scDblFinder_score", "scDblFinder_class"],
)
# %%
sc.pl.umap(
    adata,
    color=["scDblFinder_class"],
)
# %%
adata.X = adata.layers["counts"].copy()
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
adata.layers["logcounts"] = adata.X.copy()

batch_key = 'batch'
sc.pp.highly_variable_genes(
    adata, n_top_genes=2000, flavor="cell_ranger", batch_key=batch_key
)
adata_hvg = adata[:, adata.var["highly_variable"]].copy()
adata_seurat = adata_hvg.copy()
# Convert categorical columns to strings
adata_seurat.obs[batch_key] = adata_seurat.obs[batch_key].astype(str)
adata_seurat.obs_names_make_unique()
adata_seurat.var_names_make_unique()
del adata_seurat.uns
# %%
importr('Seurat')
saveSeurat = ro.functions.wrap_r_function(
    ro.r(
        """
        function(adata, filePath){
        seurat <- as.Seurat(adata, counts = "counts", data = "logcounts")
        print('converted')
        SaveSeuratRds(seurat, filePath)
        print('Saved')
        }
        """
    ),
    "saveSeurat",
)

with (
    ro.default_converter + ro.pandas2ri.converter + anndata2ri.converter
).context():
    integrated_expr = saveSeurat(adata_seurat, '../../data/external/stubenvoll_processed_seurat.rds')
# %%
importr('Seurat')
loadSeurat = ro.functions.wrap_r_function(
    ro.r(
        """
        function(location){
        expr = readRDS(location)
        return(expr)
        }
        """
    ),
    "loadSeurat",
)

with (
    ro.default_converter + ro.pandas2ri.converter + anndata2ri.converter
).context():
    integrated_expr = loadSeurat('../../data/external/stubenvoll_processed_seurat_integration.rds')
# %%
adata_seurat.X = integrated_expr
adata_seurat.layers["seurat"] = integrated_expr
# %%
# Reset the batch colours because we deleted them earlier
sc.tl.pca(adata_seurat)
sc.pp.neighbors(adata_seurat)
sc.tl.umap(adata_seurat)
sc.pl.umap(adata_seurat, color='batch', wspace=1)
# %%
adata_seurat.write_h5ad('../../data/external/stubenvoll_processed_integrated.h5ad')
adata = adata_seurat
# %%
import celltypist
from celltypist import models
immuneHigh = models.Model.load(model = 'Immune_All_High.pkl')
immuneLow = models.Model.load(model = 'Immune_All_Low.pkl')

adata_celltypist = adata.copy()  
adata_celltypist.X = adata.layers["log1p_norm"]  
# normalize to 10,000 counts per cell (celltypist wants this)
sc.pp.normalize_total(
    adata_celltypist, target_sum=10**4
)  
sc.pp.log1p(adata_celltypist)
predictions_high = celltypist.annotate(
    adata_celltypist, model=immuneHigh, majority_voting=True
)
predictions_low = celltypist.annotate(
    adata_celltypist, model=immuneLow, majority_voting=True
)
predictions_high_adata = predictions_high.to_adata().copy()
predictions_low_adata = predictions_low.to_adata().copy()

sc.pl.umap(
    predictions_high_adata,
    color=["majority_voting", "conf_score"],
    frameon=False,
    sort_order=False,
    wspace=0.3,
)
sc.pl.umap(
    predictions_low_adata,
    color=["majority_voting", "conf_score"],
    frameon=False,
    sort_order=False,
    wspace=0.3,
)

adata.write_h5ad('/home/tyj566/mnt/h/lsp-analysis/tyler/jasonSingleCell/processed/stubenvoll_celltyped.h5ad')
# %%
