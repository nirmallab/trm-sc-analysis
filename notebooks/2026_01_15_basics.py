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
    labs,
    scale_fill_manual
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

adata37 = sc.read_h5ad('../data/processed/BRI-2937_tcr.h5ad')
adata39 = sc.read_h5ad('../data/processed/BRI-2939_tcr.h5ad')
adata41 = sc.read_h5ad('../data/processed/BRI-2941_tcr.h5ad')

adata37.obs['condition'] = '2937'
adata39.obs['condition'] = '2939'
adata41.obs['condition'] = '2941'

adata_tcr = anndata.concat([adata37, adata39, adata41])
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
sc.pl.umap(adata, color=['condition', 'majority_voting_high'])
# %%
# Add in eosinophils
sc.tl.leiden(adata, resolution=0.05, random_state=1234)
sc.pl.umap(adata, color='leiden')

adata.obs['celltype_fine'] = adata.obs['majority_voting_low'].astype(str)
adata.obs['celltype_broad'] = adata.obs['majority_voting_high'].astype(str)

isEo = adata.obs['leiden'].astype(str) == '2'
adata.obs.loc[isEo, 'celltype_fine'] = 'Eosinophil'
adata.obs.loc[isEo, 'celltype_broad'] = 'Eosinophil'

sc.pl.umap(adata, color=['celltype_fine', 'celltype_broad'])

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
    ggplot(adata.obs, aes("factor(condition)", fill="factor(celltype_broad)"))
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
fig, ax = plt.subplots()
sc.pl.violin(adataReg, keys='Cxcl10', groupby='condition', rotation=90, ax=ax, show=False)
ax.set_title('Tem/Trm cytotoxic T cells')
plt.show()
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
    fig, ax = plt.subplots()
    sc.pl.dotplot(adataTSub[adataTSub.obs['condition'] == condition], genes, groupby='cat', ax=ax, show=False)
    ax.set_title(condition)
# %%
sc.pl.umap(adataT, color=['majority_voting_low'])
# %%
# TRM genes
# sc.pl.umap(adataT, color=['Cd69', 'Cd8a', 'Itgae', 'Rgs1', 'Cd101', 'Cd69', 'Havcr2', 'Pdcd1'])
sc.pl.umap(adataT, color=['Itgae', 'Cd101', 'Pdcd1', 'Havcr2' ,'Cd69', 'Rgs1'])
# %%
# %%
# celltypist markers
sc.pl.umap(adataT, color=['Gzmk', 'Cd8a', 'Ccl5'])
# %%
# TEM genes
sc.pl.umap(adataT, color=['Cx3cr1', 'Ccr5', 'Gzma', 'Il7r','Il4', 'Il5', 'Ifng'])
"""
Sell: low
Il4: Not right cells
"""
# %%
# for res in np.arange(0.1, 1, 0.2):
#     sc.tl.leiden(adataT, resolution=res ,key_added=f'leiden_{res}')
#     sc.pl.umap(adataT, color=f'leiden_{res}')
sc.tl.leiden(adataT, resolution=1.5, random_state=1234)
sc.pl.umap(adataT, color='leiden',add_outline=True,legend_loc="on data", legend_fontoutline=2)
# %%
obs = adataT.obs
obsT = obs[obs['leiden'].astype(str).isin(['11', '7', '16', '3', '2'])]
obsT['leiden'] = obsT['leiden'].astype(str)
# %%
palette = {
    '2937': '#1f77b4',  # blue
    '2939': '#ff7f0e',  # orange
    '2941': '#2ca02c',  # green
}

(
    ggplot(obsT, aes("factor(leiden)", fill="factor(condition)"))
    + geom_bar(position="fill")
    + scale_fill_manual(values=palette)
    + labs(
        title='Broad Cell Phenotypes',
        x='Leiden',
        y='Proportion',
        fill='Phenotype'
    )
)
# %%
sc.tl.dendrogram(adataT, groupby='leiden')
sc.pl.dendrogram(adataT, groupby='leiden')

print(adataT.obs['leiden'].value_counts())
# %%
sc.pl.umap(adataT, color='condition')
# %%
# Treg investigation
sc.pl.umap(adataT, color=['Foxp3', 'Ctla4', 'Il2ra', 'Cd4'])

# %%
adata_tcr.obs_names_make_unique()
adataT.obs_names_make_unique()
mdata = mu.MuData({"gex": adataT, "airr": adata_tcr})
# %%
ir.pp.index_chains(mdata)
ir.tl.chain_qc(mdata)
ir.pp.ir_dist(mdata)
ir.tl.define_clonotypes(mdata, receptor_arms="all", dual_ir="primary_only")
ir.pp.ir_dist(
    mdata,
    metric="tcrdist",
    sequence="aa",
    cutoff=15,
)
ir.tl.define_clonotype_clusters(mdata, sequence="aa", metric="tcrdist", receptor_arms="all", dual_ir="any")
ir.tl.clonotype_network(mdata, min_cells=3, sequence="aa", metric="tcrdist")
# %%
_ = ir.pl.clonotype_network(mdata, color="gex:condition", label_fontsize=9, panel_size=(12, 12), base_size=20)
# %%
ir.tl.clonal_expansion(mdata)
# %%
mu.pl.embedding(mdata, basis="gex:umap", color=["airr:clonal_expansion", "airr:clone_id_size"])
# %%
ir.tl.clonotype_convergence(mdata, key_coarse="cc_aa_tcrdist", key_fine="clone_id")
# %%
mu.pl.embedding(mdata, "gex:umap", color="airr:is_convergent")
# %%
isT = mdata.mod["gex"].obs["majority_voting_low"].isin(cellTypesT)
mdataT = mdata[:, :].copy()
mdataT.mod["gex"] = mdata.mod["gex"][isT, :]
# %%
ir.pl.clonal_expansion(
    mdata,
    groupby="gex:majority_voting_low"
)
