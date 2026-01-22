# %%
import anndata
from pathlib import Path
import muon as mu
import numpy as np
import scanpy as sc
import scirpy as ir
from cycler import cycler
from matplotlib import cm as mpl_cm
from matplotlib import pyplot as plt
# %%
adatas = []
perSampleDir = Path('../data/raw/cellRangerOuts/BRI-2939/per_sample_outs/')
for experimentDir in perSampleDir.iterdir():
    adataPath = experimentDir / 'vdj_t' / 'filtered_contig_annotations.csv'
    adatas.append(ir.io.read_10x_vdj(adataPath))
adata_tcr = anndata.concat(adatas)
adata_tcr.write_h5ad('../data/processed/BRI-2939_tcr.h5ad')
adata = sc.read_h5ad('../data/processed/BRI-2939.h5ad')
# %%
mdata = mu.MuData({"gex": adata, "airr": adata_tcr})
mdata3k = ir.datasets.wu2020_3k()

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
_ = ir.pl.clonotype_network(mdata, color="gex:sample", label_fontsize=9, panel_size=(7, 7), base_size=20)
# %%
ir.tl.clonal_expansion(mdata)
# %%
mu.pl.embedding(mdata, basis="gex:umap", color=["airr:clonal_expansion", "airr:clone_id_size"])

# %%
_ = ir.pl.clonal_expansion(mdata, target_col="clone_id", groupby="gex:majority_voting_low", breakpoints=(1, 2, 5))
# %%
sc.pl.umap(adata, color='majority_voting_low')
mu.pl.embedding(mdata, basis="gex:umap", color=["Cd3e"])
# %%