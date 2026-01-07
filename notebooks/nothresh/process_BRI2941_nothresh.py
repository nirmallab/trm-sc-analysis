# %%
import anndata
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import scanpy as sc
from scipy.stats import median_abs_deviation
import seaborn as sns
# %%
adatas = []
perSampleDir = Path('../../data/raw/cellRangerOuts/BRI-2941/per_sample_outs/')
for experimentDir in perSampleDir.iterdir():
    adataPath = experimentDir / 'count' / 'sample_filtered_feature_bc_matrix.h5'
    adataSample = sc.read_10x_h5(adataPath)
    adataSample.var_names_make_unique()
    adataSample.obs['sample'] = experimentDir.stem
    adatas.append(adataSample)
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
# %%
# Shifted logarithm normalization
scales_counts = sc.pp.normalize_total(adata, target_sum=None, inplace=False)
adata.layers["log1p_norm"] = sc.pp.log1p(scales_counts["X"], copy=True)
# %%
# The biomdal distribution after this is really quite strange
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
adata.X = adata.layers["log1p_norm"]
# %%
adata.var["highly_variable"] = adata.var["highly_deviant"]
sc.pp.pca(adata, svd_solver="arpack", use_highly_variable=True)
sc.pl.pca_scatter(adata, color="total_counts")
# %%
sc.pp.neighbors(adata)
sc.tl.umap(adata)
# %%
sc.pl.umap(adata, color="total_counts")
# %%
sc.pl.umap(adata, color="sample")
# %%
sc.pl.umap(
    adata,
    color=["total_counts", "pct_counts_mt", "scDblFinder_score", "scDblFinder_class"],
)
# %%
sc.pp.neighbors(adata, n_pcs=30)
sc.tl.umap(adata)
# %%
sc.tl.leiden(adata, key_added="leiden_res0_25", resolution=0.25)
sc.tl.leiden(adata, key_added="leiden_res0_5", resolution=0.5)
sc.tl.leiden(adata, key_added="leiden_res1", resolution=1.0)
# %%
sc.pl.umap(
    adata,
    color=["leiden_res0_25", "leiden_res0_5", "leiden_res1"],
    legend_loc="on data",
)
# %%
# %%
import celltypist
from celltypist import models
immuneHigh = models.Model.load(model="../../models/Immune_All_High_Mouse.pkl")
immuneLow = models.Model.load(model="../../models/Immune_All_Low_Mouse.pkl")

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
    wspace=0.8,
)
sc.pl.dendrogram(predictions_low_adata, groupby="majority_voting")
sc.pl.dendrogram(predictions_high_adata, groupby="majority_voting")

labelsLow = predictions_low_adata.obs[["majority_voting", "conf_score"]]
labelsLow.columns = [f'{colName}_low' for colName in labelsLow.columns]

labelsHigh = predictions_high_adata.obs[["majority_voting", "conf_score"]]
labelsHigh.columns = [f'{colName}_high' for colName in labelsHigh.columns]

adata.obs = adata.obs.join(labelsLow)
adata.obs = adata.obs.join(labelsHigh)
# %%
geneNames = adata.var.index.str.upper()

markers = ['CCR3', 'IL5R', 'CD66', 'SiglecF']

for marker in markers:
    marker = marker.upper()
    isGene = geneNames.str.startswith(marker)
    print(isGene.sum())
    print(adata.var.index[isGene])

markers = ['Ccr3', 'Il5ra', 'Siglecf']

sc.pl.umap(adata, color=markers)
# %%
thresh = 2500
overThresh = adata.obs['total_counts']>thresh

adata.obs['overThresh'] = overThresh

sc.pl.umap(adata, color=['Ccr3', 'overThresh', 'pct_counts_mt'])