# %%
import celltypist
from celltypist import models
from pathlib import Path
import scanpy as sc
# %%
immuneHigh = models.Model.load(model="../models/Immune_All_High_Mouse.pkl")
immuneLow = models.Model.load(model="../models/Immune_All_Low_Mouse.pkl")
# %%
def annotateCells(adata):
    # Prep data for cell typist
    adata_celltypist = adata.copy()  
    adata_celltypist.X = adata.layers["log1p_norm"]  
    # normalize to 10,000 counts per cell (celltypist wants this)
    sc.pp.normalize_total(
        adata_celltypist, target_sum=10**4
    )  
    sc.pp.log1p(adata_celltypist)
    predictions_high = celltypist.annotate(
        adata_celltypist, model=immuneHigh, majority_voting=True
    )
    predictions_low = celltypist.annotate(
        adata_celltypist, model=immuneLow, majority_voting=True
    )
    predictions_high_adata = predictions_high.to_adata().copy()
    predictions_low_adata = predictions_low.to_adata().copy()

    sc.pl.umap(
        predictions_high_adata,
        color=["majority_voting", "conf_score"],
        frameon=False,
        sort_order=False,
        wspace=0.3,
    )
    sc.pl.umap(
        predictions_low_adata,
        color=["majority_voting", "conf_score"],
        frameon=False,
        sort_order=False,
        wspace=0.3,
    )
    sc.pl.dendrogram(predictions_low_adata, groupby="majority_voting")
    sc.pl.dendrogram(predictions_high_adata, groupby="majority_voting")
    
    labelsLow = predictions_low_adata.obs[["majority_voting", "conf_score"]]
    labelsLow.columns = [f'{colName}_low' for colName in labelsLow.columns]

    labelsHigh = predictions_high_adata.obs[["majority_voting", "conf_score"]]
    labelsHigh.columns = [f'{colName}_high' for colName in labelsHigh.columns]

    adata.obs = adata.obs.join(labelsLow)
    adata.obs = adata.obs.join(labelsHigh)
    return adata
# %%
samples = ['BRI-2937', 'BRI-2939', 'BRI-2941']
saveDir = Path('../data/processed')
for sample in samples:
    adataPath = Path(f'../data/interim/{sample}/{sample}_reduced.h5ad')
    adata = sc.read_h5ad(adataPath)

    adataLabeled = annotateCells(adata)

    newPath = saveDir / f'{sample}_annotated.h5ad'

    adata.write_h5ad(newPath)
# %%
