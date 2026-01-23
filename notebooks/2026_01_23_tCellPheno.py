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
# %%
def geneSearch(gene):
    genes = adata.var.index
    isGene = genes.str.startswith(gene)
    print(genes[isGene])

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
cellTypesT = ['NK cells',
'CD16- NK cells',
'Cycling NK cells',
'gamma-delta T cells',
'ILC3',
'Double-positive thymocytes',
'Tem/Trm cytotoxic T cells',
'Regulatory T cells']
adataT = adata[adata.obs['majority_voting_low'].isin(cellTypesT)]
sc.pp.normalize_per_cell(adataT)
sc.pp.log1p(adataT)
sc.pp.pca(adataT, svd_solver="arpack", use_highly_variable=True)
sc.tl.umap(adataT)
# %%
sc.pl.umap(adataT, color='Cd8a')
# %%
sc.tl.leiden(adataT, resolution=1.5, random_state=1234)

sc.pl.umap(adataT, color='leiden',add_outline=True,legend_loc="on data", legend_fontoutline=2)
# %%
sc.tl.dendrogram(adataT, groupby='leiden')
sc.pl.dendrogram(adataT, groupby='leiden')
# %%
sc.pl.violin(adataT, 'Cd8a', groupby='leiden')
# %%
adataT.obs['Cd8Pos'] =   (adataT[:,'Cd8a'].X>0.2).todense()
adataCd8 = adataT[adataT.obs['Cd8Pos']]
sc.pp.normalize_per_cell(adataCd8)
sc.pp.log1p(adataCd8)
sc.pp.highly_variable_genes(adataCd8, flavor="cell_ranger", n_top_genes=5000)
sc.tl.pca(adataCd8)
sc.pp.neighbors(adataCd8)
sc.tl.umap(adataCd8)
# %%
sc.tl.leiden(adataCd8, resolution=.25)
sc.pl.umap(adataCd8, color='leiden',add_outline=False,legend_loc="on data", legend_fontoutline=4)
sc.pl.umap(adataCd8, color=['Pdcd1', 'Itgae', 'majority_voting_low', 'leiden'], wspace=0.5)
# %%
sc.tl.rank_genes_groups(adataCd8, 'leiden', method='wilcoxon')
# %%
sc.pl.rank_genes_groups(adataCd8)

# %%
sc.pl.umap(adataCd8, color=['majority_voting_low', 'leiden'], wspace=0.5)
# %%
adata_tcr.obs_names_make_unique()
adataCd8.obs_names_make_unique()
mdata = mu.MuData({"gex": adataCd8, "airr": adata_tcr})

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
ir.tl.clonal_expansion(mdata)
# %%
ax = ir.pl.group_abundance(mdata, groupby="airr:chain_pairing", target_col="gex:leiden")
# %%
