# %%
import anndata
import os
from pathlib import Path
from tqdm import tqdm

import scanpy as sc
import scvelo as scv

# %%
adata37 = sc.read_h5ad('../data/processed/BRI-2937.h5ad')
adata39 = sc.read_h5ad('../data/processed/BRI-2939.h5ad')
adata41 = sc.read_h5ad('../data/processed/BRI-2941.h5ad')

adata37.obs['condition'] = '2937'
adata39.obs['condition'] = '2939'
adata41.obs['condition'] = '2941'

adata = anndata.concat([adata37, adata39, adata41])

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
cellTypesT = ['NK cells',
'CD16- NK cells',
'Cycling NK cells',
'gamma-delta T cells',
'ILC3',
'Double-positive thymocytes',
'Tem/Trm cytotoxic T cells',
'Regulatory T cells']
adataT = adata[adata.obs['majority_voting_low'].isin(cellTypesT)].copy()
adataT.obs['Cd8Pos'] =   (adataT[:,'Cd8a'].X>0.2).todense()
adataT = adata[adata.obs['majority_voting_low'].isin(cellTypesT)]

adataT.obs['Cd8Pos'] =   (adataT[:,'Cd8a'].X>0.2).todense()
adataCd8 = adataT[adataT.obs['Cd8Pos']]
sc.pp.normalize_per_cell(adataCd8)
sc.pp.log1p(adataCd8)
# sc.pp.highly_variable_genes(adataCd8, flavor="cell_ranger", n_top_genes=5000)
sc.tl.pca(adataCd8)
sc.pp.neighbors(adataCd8)
sc.tl.umap(adataCd8)
# %%
sc.tl.leiden(adataCd8, resolution=.5, random_state=1234)
sc.pl.umap(adataCd8, color='leiden',add_outline=False,legend_loc="on data", legend_fontoutline=4)
sc.pl.umap(adataCd8, color='condition')
sc.pl.umap(adataCd8, color=['Pdcd1', 'Itgae'], wspace=0.15)
# %%
scv.tl.velocity_graph(adataCd8)
# %%
scv.pl.velocity_embedding(
    adataCd8,
    basis='umap',
    arrow_length=3,
    arrow_size=2,
    color='leiden'
)