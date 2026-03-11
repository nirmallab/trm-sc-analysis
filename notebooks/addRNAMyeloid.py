# %%
import anndata
import muon as mu
import os
import pandas as pd
from pathlib import Path
import scanpy as sc
import scirpy as ir
import scvelo as scv
from tqdm import tqdm
# %%
def geneSearch(gene):
    genes = adata.var.index
    isGene = genes.str.startswith(gene)
    print(genes[isGene])
# %%
adata37 = sc.read_h5ad('../data/processed/BRI-2937.h5ad')
adata39 = sc.read_h5ad('../data/processed/BRI-2939.h5ad')
adata41 = sc.read_h5ad('../data/processed/BRI-2941.h5ad')

adata37.obs['condition'] = '2937'
adata39.obs['condition'] = '2939'
adata41.obs['condition'] = '2941'

adata = anndata.concat([adata37, adata39, adata41])
# %%
adataMyeloid = sc.read_h5ad('../data/external/BRI-2937_myeloid_only.h5ad')
# %%
obs1 = adata.obs.copy()
obs2 = adataMyeloid.obs.copy()
obs2['myeloid_anno'] = obs2['myeloid_anno'].astype("string")
obsM = obs1.join(obs2['myeloid_anno'])

obsM['myeloid_anno'].value_counts()
adata.obs = obsM
# %%
sample2Experiment = {k:f'BRI-{v}' for k,v in zip(adata.obs['sample'], adata.obs['condition'])}


adatas = []
for sample in tqdm(adata.obs['sample'].unique()):
    experiment = sample2Experiment[sample]
    experimentCountDir = Path(f'/home/tyj566/mnt/h/lsp-analysis/tyler/jasonSingleCell/{experiment}')

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
scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata)
# %%
scv.tl.velocity(adata, mode='stochastic')
# %%
scv.tl.velocity_graph(adata)
# %%
adataSub = adata[~adata.obs['myeloid_anno'].isna()]
sc.pp.normalize_per_cell(adataSub)
sc.pp.log1p(adataSub)
# sc.pp.highly_variable_genes(adataCd8, flavor="cell_ranger", n_top_genes=5000)
sc.tl.pca(adataSub)
sc.pp.neighbors(adataSub)
sc.tl.umap(adataSub)
# %%
sc.pl.umap(adataSub, color='myeloid_anno')
# %%
scv.tl.velocity_graph(adataSub)
# %%
scv.pl.velocity_embedding(
    adataSub,
    basis='umap',
    arrow_length=3,
    arrow_size=2,
    color='myeloid_anno'
)
# %%
isPheno = ~adata.obs['myeloid_anno'].isin(['cDC1', 'pDCs', 'fibrobalst'])
adataSub = adata[~adata.obs['myeloid_anno'].isna() & isPheno]
# %%
sc.pp.normalize_per_cell(adataSub)
sc.pp.log1p(adataSub)
# sc.pp.highly_variable_genes(adataCd8, flavor="cell_ranger", n_top_genes=5000)
sc.tl.pca(adataSub)
sc.pp.neighbors(adataSub)
sc.tl.umap(adataSub)
# %%
sc.pl.umap(adataSub, color='myeloid_anno')
# %%
scv.tl.velocity_graph(adataSub)
# %%
scv.pl.velocity_embedding(
    adataSub,
    basis='umap',
    arrow_length=3,
    arrow_size=2,
    color='myeloid_anno'
)
# %%