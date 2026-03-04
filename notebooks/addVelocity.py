# %%
import anndata
import muon as mu
import os
from pathlib import Path
import scanpy as sc
import scirpy as ir
import scvelo as scv
from tqdm import tqdm
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
# %%
adata37 = sc.read_h5ad('../data/processed/BRI-2937_tcr.h5ad')
adata39 = sc.read_h5ad('../data/processed/BRI-2939_tcr.h5ad')
adata41 = sc.read_h5ad('../data/processed/BRI-2941_tcr.h5ad')

adata37.obs['condition'] = '2937'
adata39.obs['condition'] = '2939'
adata41.obs['condition'] = '2941'

adata_tcr = anndata.concat([adata37, adata39, adata41])
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
ir.tl.clonotype_network(mdata, min_cells=3, sequence="aa", metric="tcrdist")

# %%
_ = ir.pl.clonotype_network(mdata, color="gex:sample", label_fontsize=9, panel_size=(7, 7), base_size=20)
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