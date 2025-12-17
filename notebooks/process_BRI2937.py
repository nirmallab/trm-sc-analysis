# %%
import anndata
import numpy as np
from pathlib import Path
import scanpy as sc
from scipy.stats import median_abs_deviation
import seaborn as sns
# %%
adatas = []
perSampleDir = Path('../data/raw/cellRangerOuts/BRI-2937/per_sample_outs/')
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
# %% SoupX - Skipping for now
# import logging
# import anndata2ri
# import rpy2.rinterface_lib.callbacks as rcb
# import rpy2.robjects as ro

# rcb.logger.setLevel(logging.ERROR)
# ro.pandas2ri.activate()
# anndata2ri.activate()
# %load_ext rpy2.ipython
# #%%
# %%R
# library(SoupX)
# # %%
# adata_pp = adata.copy()
# sc.pp.normalize_per_cell(adata_pp)
# sc.pp.log1p(adata_pp)
# # %%
# sc.pp.pca(adata_pp)
# sc.pp.neighbors(adata_pp)
# sc.tl.leiden(adata_pp, key_added="soupx_groups")

# # Preprocess variables for SoupX
# soupx_groups = adata_pp.obs["soupx_groups"]
# # %%
# del adata_pp
# # %%
# cells = adata.obs_names
# genes = adata.var_names
# data = adata.X.T
# # %%
# adata_raw = sc.read_10x_h5('/data/trm-sc-analysis/data/raw/cellRangerOuts/BRI-2937/count/raw_feature_bc_matrix.h5')
# # %%
# adata_raw.var_names_make_unique()
# data_tod = adata_raw.X.T
# # %%
# del adata_raw
# # %%
# %%R -i data -i data_tod -i genes -i cells -i soupx_groups -o out 

# # specify row and column names of data
# rownames(data) = genes
# colnames(data) = cells
# # ensure correct sparse format for table of counts and table of droplets
# data <- as(data, "sparseMatrix")
# data_tod <- as(data_tod, "sparseMatrix")

# # Generate SoupChannel Object for SoupX 
# sc = SoupChannel(data_tod, data, calcSoupProfile = FALSE)

# # Add extra meta data to the SoupChannel object
# soupProf = data.frame(row.names = rownames(data), est = rowSums(data)/sum(data), counts = rowSums(data))
# sc = setSoupProfile(sc, soupProf)
# # Set cluster information in SoupChannel
# sc = setClusters(sc, soupx_groups)

# # Estimate contamination fraction
# sc  = autoEstCont(sc, doPlot=FALSE)
# # Infer corrected table of counts and rount to integer
# out = adjustCounts(sc, roundToInt = TRUE)
# %%
import logging

import anndata as ad
import anndata2ri
import rpy2.rinterface_lib.callbacks as rcb
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
adata.write_h5ad('../data/interim/BRI-2937_qc.h5ad')
# %%
adata = sc.read_h5ad('../data/interim/BRI-2937_qc.h5ad')
# %%
p1 = sns.histplot(adata.obs["total_counts"], bins=100, kde=False)