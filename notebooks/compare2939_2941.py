# %%
import anndata
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import scanpy as sc
from scipy.stats import median_abs_deviation
import seaborn as sns
# %%
def initMetrics(adata):
    # mitochondrial genes
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    # ribosomal genes
    adata.var["ribo"] = adata.var_names.str.upper().str.startswith(("RPS", "RPL"))
    # hemoglobin genes.
    adata.var["hb"] = adata.var_names.str.upper().str.contains("^HB[^(P)]")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True
    )
    return adata
# %%
adatas = []
perSampleDir = Path('../data/raw/cellRangerOuts/BRI-2941/per_sample_outs/')
for experimentDir in perSampleDir.iterdir():
    adataPath = experimentDir / 'count' / 'sample_filtered_feature_bc_matrix.h5'
    adataSample = sc.read_10x_h5(adataPath)
    adataSample.var_names_make_unique()
    adataSample.obs['sample'] = experimentDir.stem
    adatas.append(adataSample)
adata2941 = anndata.concat(adatas)
adata2941 = initMetrics(adata2941)
# %%
adatas = []
perSampleDir = Path('../data/raw/cellRangerOuts/BRI-2939/per_sample_outs/')
for experimentDir in perSampleDir.iterdir():
    adataPath = experimentDir / 'count' / 'sample_filtered_feature_bc_matrix.h5'
    adataSample = sc.read_10x_h5(adataPath)
    adataSample.var_names_make_unique()
    adataSample.obs['sample'] = experimentDir.stem
    adatas.append(adataSample)
adata2939 = anndata.concat(adatas)
adata2939 = initMetrics(adata2939)
# %%
plt.figure()
plt.hist(adata2941.obs['total_counts'], bins=100)
plt.title('2941')

plt.figure()
plt.hist(adata2939.obs['total_counts'], bins=100)
plt.title('2939')

# %%
perSampleDir = Path('../data/raw/cellRangerOuts/BRI-2939/per_sample_outs/')
for experimentDir in perSampleDir.iterdir():
    adataPath = experimentDir / 'count' / 'sample_filtered_feature_bc_matrix.h5'
    adataSample = sc.read_10x_h5(adataPath)
    adataSample.var_names_make_unique()
    adataSample.obs['sample'] = experimentDir.stem
    adataSample = initMetrics(adataSample)
    plt.figure()
    plt.hist(adataSample.obs['total_counts'], bins=100)