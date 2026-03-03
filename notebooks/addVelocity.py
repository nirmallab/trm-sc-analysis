# %%
import anndata
import os
from pathlib import Path

import scanpy as sc
import scvelo as scv
# %%
experiment = 'BRI-2937'
experimentCountDir = Path(f'/home/tyj566/mnt/h/lsp-analysis/tyler/jasonSingleCell/{experiment}')
# sampleCountDir = f'/home/tyj566/mnt/h/lsp-analysis/tyler/jasonSingleCell/{sample}/per_sample_outs/4/count/velocyto'

adata = sc.read_h5ad(f'../data/processed/{experiment}.h5ad')

adatas = []
for sample in adata.obs['sample'].unique():
    adataSample = adata[adata.obs['sample'] == sample]
    veloSubFolder = f'per_sample_outs/{sample}/count/velocyto'
    veloPath = experimentCountDir / veloSubFolder
    veloOut = veloPath / os.listdir(veloPath)[0]
    loomDat = sc.read(veloOut)

    adataSample.var_names_make_unique()
    adataSample = scv.utils.merge(adataSample, loomDat)
    adatas.append(adataSample)

adata = anndata.concat(adatas, join='outer')
# %%
