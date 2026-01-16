# %%
import anndata
from pathlib import Path
import muon as mu
import numpy as np
import scanpy as sc
import scirpy as ir
import seaborn as sns
from matplotlib import cm as mpl_cm
from matplotlib import pyplot as plt
from plotnine import (
    ggplot,
    aes,
    geom_bar,
    labs
)
# %%
def subsetCellType(adata, cellType, voting='low'):
    return adata[adata.obs[f'majority_voting_{voting}'] == cellType]

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
sc.pp.filter_genes(adata, min_cells=10)
sc.pp.filter_cells(adata, min_genes=100)
sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor="cell_ranger", n_top_genes=5000)
sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
# %%
sc.pl.umap(adata, color=['condition', 'majority_voting_low'])
# %%
sc.pl.umap(adata, color=['majority_voting_low'])
# %%
plt.figure(figsize=(10,10))

cellTypesT = ['NK cells',
'CD16- NK cells',
'Cycling NK cells',
'gamma-delta T cells',
'ILC3',
'Double-positive thymocytes',
'Tem/Trm cytotoxic T cells',
'Regulatory T cells']
adataT = adata[adata.obs['majority_voting_low'].isin(cellTypesT)]
# %%
(
    ggplot(adataT.obs, aes("factor(condition)", fill="factor(majority_voting_low)"))
    + geom_bar(position="fill")
    + labs(title='T-Cell Associated Phenotypes',
           x='Experiment',
           y='Proportion',
           fill='Phenotype')
)
# %%
(
    ggplot(adata.obs, aes("factor(condition)", fill="factor(majority_voting_high)"))
    + geom_bar(position="fill")
    + labs(title='Broad Cell Phenotypes',
           x='Experiment',
           y='Proportion',
           fill='Phenotype')
)
# %%
sc.pp.normalize_per_cell(adataT)
sc.pp.log1p(adataT)
sc.pp.pca(adataT, svd_solver="arpack", use_highly_variable=True)
sc.tl.umap(adataT)
# %%
sc.pl.umap(adataT, color = ['majority_voting_low', 'condition'])
# %%
sc.pl.violin(adataT, keys='Cxcl9', groupby='majority_voting_low', rotation=90)
sc.pl.violin(adataT, keys='Cxcl10', groupby='majority_voting_low', rotation=90)
# sc.pl.violin(adataT, keys='Cxcl11', groupby='majority_voting_low', rotation=90)
# %%
adataReg = subsetCellType(adataT, cellType = 'Tem/Trm cytotoxic T cells')
sc.pl.violin(adataReg, keys='Cxcl10', groupby='condition', rotation=90)
# %%
sc.pl.umap(adataT, color = ['Cxcl10', 'condition'])
# %%
sc.pl.umap(adataT, color=['Cd8a', 'Cd4', 'Foxp3'])
# %%
geneExpr = np.array(adataT[:,'Foxp3'].X.todense()[:,0]).ravel()
plt.hist(geneExpr, bins=50)
# %%
adataT.obs['CD8Pos'] =   (adataT[:,'Cd8a'].X>0.2).todense()
adataT.obs['CD4PosFoxp3Neg'] =   (adataT[:,'Cd4'].X>0.2).todense() & (adataT[:,'Foxp3'].X<0.2).todense()

sc.pl.umap(adataT, color=['CD8Pos', 'CD4PosFoxp3Neg'])
# %%
adataTSub = adataT[adataT.obs['CD8Pos'] | adataT.obs['CD4PosFoxp3Neg']]
adataTSub.obs['cat'] = 'CD8Pos'
adataTSub.obs['cat'] = adataTSub.obs['cat'].astype(str)
adataTSub.obs.loc[adataTSub.obs['CD4PosFoxp3Neg'], 'cat'] = 'CD4PosFoxp3Neg'
sc.tl.rank_genes_groups(adataTSub, 'cat', method='wilcoxon')
# %%
sc.pl.rank_genes_groups_dotplot(adataTSub, n_genes=10)
# %%
rankingDfCD8Pos = sc.get.rank_genes_groups_df(adataTSub, group='CD8Pos')
rankingDfCD4 = sc.get.rank_genes_groups_df(adataTSub, group='CD4PosFoxp3Neg')
rankingDfCD8Pos = rankingDfCD8Pos.set_index('names')
# %%
geneSearch('Ccr5')
# %%
genes = ['Ifng', 'Gzmb', 'Fasl']
print(rankingDfCD8Pos.loc[genes])
conditions = adata.obs['condition'].unique()
for condition in conditions:
    sc.pl.dotplot(adataTSub[adataTSub.obs['condition'] == condition], genes, groupby='cat', ax=ax)
# %%
sc.pl.umap(adataT, color=['majority_voting_low', 'Cd101', 'Itgae'])
sc.pl.umap(adataT, color=['Cd69', 'Cd8a', 'Itgae', 'Rgs1', 'Cd101'])
# %%
sc.pl.umap(adataT, color=['Gzmk', 'Ccr5'])
