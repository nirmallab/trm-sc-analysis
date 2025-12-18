adata_pp = adata.copy()
sc.pp.normalize_per_cell(adata_pp)
sc.pp.log1p(adata_pp)
# %%
sc.pp.pca(adata_pp)
sc.pp.neighbors(adata_pp)
sc.tl.leiden(adata_pp, key_added="soupx_groups")

# Preprocess variables for SoupX
soupx_groups = adata_pp.obs["soupx_groups"]
# %%
del adata_pp
# %%
cells = adata.obs_names
genes = adata.var_names
data = adata.X.T
# %%
# Get raw data
adatasRaw = []
for experimentDir in perSampleDir.iterdir():
    adataPath = experimentDir / 'count' / 'sample_raw_feature_bc_matrix.h5'
    adataSample = sc.read_10x_h5(adataPath)
    adataSample.var_names_make_unique()
    adataSample.obs['sample'] = experimentDir.stem
    adatasRaw.append(adataSample)
adata_raw = anndata.concat(adatasRaw)
adata_raw.var_names_make_unique()
data_tod = adata_raw.X.T
# %%
del adata_raw
# %%
import anndata2ri
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
importr('SoupX')
# %%
soupxAmbient = ro.functions.wrap_r_function(
    ro.r(
        """
function(data, data_tod, genes, cells, soupx_groups){
rownames(data) = genes
colnames(data) = cells
# ensure correct sparse format for table of counts and table of droplets
data <- as(data, "sparseMatrix")
data_tod <- as(data_tod, "sparseMatrix")

# Generate SoupChannel Object for SoupX 
sc = SoupChannel(data_tod, data, calcSoupProfile = FALSE)

# Add extra meta data to the SoupChannel object
soupProf = data.frame(row.names = rownames(data), est = rowSums(data)/sum(data), counts = rowSums(data))
sc = setSoupProfile(sc, soupProf)
# Set cluster information in SoupChannel
sc = setClusters(sc, soupx_groups)

# Estimate contamination fraction
sc  = autoEstCont(sc, doPlot=FALSE)
# Infer corrected table of counts and rount to integer
out = adjustCounts(sc, roundToInt = TRUE)
return(out)
}
"""
    ),
    "soupxAmbient",
)
with (
    ro.default_converter + ro.pandas2ri.converter + anndata2ri.converter
).context():
    soupxOut = soupxAmbient(data, data_tod, genes, cells, soupx_groups)
# %% 
adata.layers["counts"] = adata.X
adata.layers["soupX_counts"] = soupxOut.T
adata.X = adata.layers["soupX_counts"]
# %%
print(f"Total number of genes: {adata.n_vars}")
# Min 20 cells - filters out 0 count genes
sc.pp.filter_genes(adata, min_cells=20)
print(f"Number of genes after cell filter: {adata.n_vars}")