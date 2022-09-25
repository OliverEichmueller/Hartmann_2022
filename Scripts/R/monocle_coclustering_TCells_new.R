# Transfer ref to monocle, calculate gene modules


library(monocle3)
cds.ref <- readRDS('TCells/monocle/cds_ref_211023.Rds')

# re-do for T-cells --------------------------------------------------------

# extract genes
genes <- as.data.frame(rownames(sc.obj.TCells@assays$RNA), 
                       row.names = rownames(sc.obj.TCells@assays$RNA))
colnames(genes) <- "gene_short_name"

# extract cells
cells <- as.data.frame(
  sc.obj.TCells@assays[["RNA"]]@data@Dimnames[[2]], 
  row.names = sc.obj.TCells@assays[["RNA"]]@data@Dimnames[[2]])
colnames(cells) <- "barcode"

# extract expression matrix
expression_matrix <- sc.obj.TCells@assays[["RNA"]]@data
expression_matrix <- expression_matrix[rownames(sc.obj.TCells@assays$RNA), ]

# Assemble cell data set object
cds.integration <- new_cell_data_set(expression_matrix,
                                     cell_metadata = cells,
                                     gene_metadata = genes)

# transfer dataset information

cds.integration@colData$orig.ident <- 
  sc.obj.TCells@meta.data[cds.integration@colData$barcode, 'orig.ident']
cds.integration@colData$integrated_snn_res.0.3 <- 
  sc.obj.TCells@meta.data[cds.integration@colData$barcode, 'integrated_snn_res.0.3']
cds.integration@colData$integrated_snn_res.0.25 <- 
  sc.obj.TCells@meta.data[cds.integration@colData$barcode, 'integrated_snn_res.0.25']



# preprocess cds with variable features from seurat
cds.integration@int_colData@listData$reducedDims$PCA <- 
  sc.obj.TCells@reductions[["pca"]]@cell.embeddings[row.names(cds.integration@colData),]

# perform umap dimensionality reduction

cds.integration@int_colData@listData$reducedDims$UMAP <- 
  sc.obj.TCells@reductions[["umap"]]@cell.embeddings[row.names(cds.integration@colData),]


plot_cells(cds.integration, color_cells_by  = 'integrated_snn_res.0.25')

cds.integration@preprocess_aux$gene_loadings <- sc.obj.TCells@reductions$pca@feature.loadings

# Generate and Aggregate based on modules of reference -------------------------------------

modules_TC <- find_gene_modules(cds.integration, resolution = .05, random_seed = 1984)

modules_TC$module %>% table()
modules_TC$id %>% unique() 
write.csv(modules_TC, 'TCells/monocle/modules_TC.csv')
modules_TC <- read.csv('TCells_new/monocle/modules_TC.csv', row.names = 1)

aggr_mod_ref <- aggregate_gene_expression(cds.ref, gene_group_df = modules_TC, 
                          cell_group_df = cds.ref@colData[,c("barcode", "functional.cluster")])

aggr_mod_Tc <- aggregate_gene_expression(cds.integration, gene_group_df = modules_TC, 
                                          cell_group_df = cds.integration@colData[,c("barcode", "integrated_snn_res.0.25")])


col2 = colorRampPalette(c('#67001F', '#B2182B', '#D6604D', '#F4A582',
                          '#FDDBC7', '#FFFFFF', '#D1E5F0', '#92C5DE',
                          '#4393C3', '#2166AC', '#053061'))

pdf('TCells_new/monocle/corrplot_modules_res.0.25.pdf')
pl = corrplot::corrplot(cor(cbind(aggr_mod_ref, aggr_mod_Tc) %>% as.matrix())
                   , order = 'hclust', addrect = 8, col = rev(col2(200)))
dev.off()


pdf('TCells_new/monocle/modules_res025.pdf')
pl <- corrplot::corrplot(t(as.matrix(aggr_mod_Tc)), is.corr = F, col = rev(col2(200)))
dev.off()

write.xlsx(list(`module genes lymphoid cells` = modules_TC %>% arrange(module),
                          `aggregated modules lymphoid` = as.data.frame(aggr_mod_Tc),
                          `aggregated modules Andreatta` = as.data.frame(aggr_mod_ref),
                          `Correlation of modules` = cor(cbind(aggr_mod_ref, aggr_mod_Tc) %>% as.matrix())),
                     rowNames = TRUE, file = 'TCells_new/monocle/Supp_Table_corlymphoid.xlsx',
           overwrite = T)

saveRDS(cds.integration, 'TCells_new/monocle/cds_TCells_211028_filtered.Rds')
cds.integration <- readRDS('TCells_new/monocle/cds_TCells_211028_filtered.Rds')
saveRDS(cds.ref, 'TCells/monocle/cds_ref_211023.Rds')
cds.ref <- readRDS('TCells_new/monocle/cds_ref_211023.Rds')

