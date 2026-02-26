library(Seurat)
batch_key = "batch"
seurat = LoadSeuratRds('./testSeurat.rds')
print('Seurat object')
batch_list <- SplitObject(seurat, split.by = batch_key)
print('Split Object')
anchors <- FindIntegrationAnchors(batch_list, anchor.features = rownames(seurat))
print('Finding anchors')
integrated <- IntegrateData(anchors)
# Extract the integrated expression matrix
integrated_expr <- GetAssayData(integrated)
# Make sure the rows and columns are in the same order as the original object
integrated_expr <- integrated_expr[rownames(seurat), colnames(seurat)]
# Transpose the matrix to AnnData format
integrated_expr <- t(integrated_expr)
print(integrated_expr[1:10, 1:10])
saveRDS(integrated_expr, file = "saveRDSIntegrated_expr.rds")
SaveSeuratRds(integrated_expr, './integratedExpression.rds')
