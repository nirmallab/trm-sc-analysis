# %%
import muon as mu
import numpy as np
import scanpy as sc
import scirpy as ir
from cycler import cycler
from matplotlib import cm as mpl_cm
from matplotlib import pyplot as plt

# temporary fix for deprecated matplotlib functionality
import IPython.display
from matplotlib_inline.backend_inline import set_matplotlib_formats

IPython.display.set_matplotlib_formats = set_matplotlib_formats

sc.set_figure_params(figsize=(4, 4))
sc.settings.verbosity = 2  # verbosity: errors (0), warnings (1), info (2), hints (3)
# %%
mdata = ir.datasets.wu2020_3k()
# %%
# Preprocess gex
sc.pp.filter_genes(mdata["gex"], min_cells=10)
sc.pp.filter_cells(mdata["gex"], min_genes=100)
sc.pp.normalize_per_cell(mdata["gex"])
sc.pp.log1p(mdata["gex"])
sc.pp.highly_variable_genes(mdata["gex"], flavor="cell_ranger", n_top_genes=5000)
sc.tl.pca(mdata["gex"])
sc.pp.neighbors(mdata["gex"])
mdata["gex"].obsm["X_umap"] = mdata["gex"].obsm["X_umap_orig"]

mapping = {
    "nan": "other",
    "3.1-MT": "other",
    "4.1-Trm": "CD4_Trm",
    "4.2-RPL32": "CD4_RPL32",
    "4.3-TCF7": "CD4_TCF7",
    "4.4-FOS": "CD4_FOSS",
    "4.5-IL6ST": "CD4_IL6ST",
    "4.6a-Treg": "CD4_Treg",
    "4.6b-Treg": "CD4_Treg",
    "8.1-Teff": "CD8_Teff",
    "8.2-Tem": "CD8_Tem",
    "8.3a-Trm": "CD8_Trm",
    "8.3b-Trm": "CD8_Trm",
    "8.3c-Trm": "CD8_Trm",
    "8.4-Chrom": "other",
    "8.5-Mitosis": "other",
    "8.6-KLRB1": "other",
}
mdata["gex"].obs["cluster"] = mdata["gex"].obs["cluster_orig"].map(mapping)

mdata.update()
# %%
mu.pl.embedding(
    mdata,
    basis="gex:umap",
    color=["gex:sample", "gex:patient", "gex:cluster"],
    ncols=3,
    wspace=0.7,
)
mu.pl.embedding(
    mdata,
    basis="gex:umap",
    color=["CD8A", "CD4", "FOXP3"],
    ncols=3,
    wspace=0.7,
)
# %%
ir.pp.index_chains(mdata)
# %%
ir.tl.chain_qc(mdata)
# %%
_ = ir.pl.group_abundance(mdata, groupby="airr:receptor_subtype", target_col="gex:source")
# %%
_ = ir.pl.group_abundance(mdata, groupby="airr:chain_pairing", target_col="gex:source")
# %%
print(
    "Fraction of cells with more than one pair of TCRs: {:.2f}".format(
        np.sum(mdata.obs["airr:chain_pairing"].isin(["extra VJ", "extra VDJ", "two full chains", "multichain"]))
        / mdata["airr"].n_obs
    )
)
# %%
mu.pl.embedding(mdata, basis="gex:umap", color="airr:chain_pairing", groups="multichain")
# %%
mu.pp.filter_obs(mdata, "airr:chain_pairing", lambda x: x != "multichain")
# %%
mu.pp.filter_obs(mdata, "airr:chain_pairing", lambda x: x != "multichain")
# %%
mu.pp.filter_obs(mdata, "airr:chain_pairing", lambda x: ~np.isin(x, ["orphan VDJ", "orphan VJ"]))
# %%
ax = ir.pl.group_abundance(mdata, groupby="airr:chain_pairing", target_col="gex:source")
# %%
# using default parameters, `ir_dist` will compute nucleotide sequence identity
ir.pp.ir_dist(mdata)
ir.tl.define_clonotypes(mdata, receptor_arms="all", dual_ir="primary_only")
# %%
ir.tl.clonotype_network(mdata, min_cells=2)
# %%
mdata.obs.groupby("gex:source", dropna=False).size()
# %%
_ = ir.pl.clonotype_network(mdata, color="gex:source", base_size=20, label_fontsize=9, panel_size=(7, 7))
# %%
ir.pp.ir_dist(
    mdata,
    metric="tcrdist",
    sequence="aa",
    cutoff=15,
)
# %%
ir.tl.define_clonotype_clusters(mdata, sequence="aa", metric="tcrdist", receptor_arms="all", dual_ir="any")
# %%
ir.tl.clonotype_network(mdata, min_cells=3, sequence="aa", metric="tcrdist")
# %%
_ = ir.pl.clonotype_network(mdata, color="gex:patient", label_fontsize=9, panel_size=(7, 7), base_size=20)
# %%
with ir.get.airr_context(mdata, "junction_aa", ["VJ_1", "VDJ_1", "VJ_2", "VDJ_2"]):
    cdr3_ct_159 = (
        # TODO astype(str) is required due to a bug in pandas ignoring `dropna=False`. It seems fixed in pandas 2.x
        mdata.obs.loc[lambda x: x["airr:cc_aa_tcrdist"] == "159"]
        .astype(str)
        .groupby(
            [
                "VJ_1_junction_aa",
                "VDJ_1_junction_aa",
                "VJ_2_junction_aa",
                "VDJ_2_junction_aa",
                "airr:receptor_subtype",
            ],
            observed=True,
            dropna=False,
        )
        .size()
        .reset_index(name="n_cells")
    )
cdr3_ct_159
# %%
# Force to have the save variable(V) gene
ir.tl.define_clonotype_clusters(
    mdata,
    sequence="aa",
    metric="tcrdist",
    receptor_arms="all",
    dual_ir="any",
    same_v_gene=True,
    key_added="cc_aa_tcrdist_same_v",
)
# find clonotypes with more than one `clonotype_same_v`
ct_different_v = mdata.obs.groupby("airr:cc_aa_tcrdist").apply(lambda x: x["airr:cc_aa_tcrdist_same_v"].nunique() > 1)
ct_different_v = ct_different_v[ct_different_v].index.values.tolist()
ct_different_v
# %%
with ir.get.airr_context(mdata, "v_call", ["VJ_1", "VDJ_1"]):
    ct_different_v_df = (
        mdata.obs.loc[
            lambda x: x["airr:cc_aa_tcrdist"].isin(ct_different_v),
            [
                "airr:cc_aa_tcrdist",
                "airr:cc_aa_tcrdist_same_v",
                "VJ_1_v_call",
                "VDJ_1_v_call",
            ],
        ]
        .sort_values("airr:cc_aa_tcrdist")
        .drop_duplicates()
        .reset_index(drop=True)
    )
ct_different_v_df
# %%
ir.tl.clonal_expansion(mdata)
# %%
mu.pl.embedding(mdata, basis="gex:umap", color=["airr:clonal_expansion", "airr:clone_id_size"], wspace=0.3)
# %%
_ = ir.pl.clonal_expansion(mdata, target_col="clone_id", groupby="gex:cluster", breakpoints=(1, 2, 5), normalize=False)
# %%
ir.pl.clonal_expansion(mdata, target_col="clone_id", groupby="gex:cluster")