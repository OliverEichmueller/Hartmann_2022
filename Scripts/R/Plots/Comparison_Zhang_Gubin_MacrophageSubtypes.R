# Process Zhang data for comparison to macrophages

# load required libraries ------------------------------------------------------


add.flag <- function(pheatmap,
                     kept.labels,
                     repel.degree) {
  
  # repel.degree = number within [0, 1], which controls how much 
  #                space to allocate for repelling labels.
  ## repel.degree = 0: spread out labels over existing range of kept labels
  ## repel.degree = 1: spread out labels over the full y-axis
  
  heatmap <- pheatmap$gtable
  
  new.label <- heatmap$grobs[[which(heatmap$layout$name == "row_names")]] 
  
  # keep only labels in kept.labels, replace the rest with ""
  new.label$label <- ifelse(new.label$label %in% kept.labels, 
                            new.label$label, "")
  
  # calculate evenly spaced out y-axis positions
  repelled.y <- function(d, d.select, k = repel.degree){
    # d = vector of distances for labels
    # d.select = vector of T/F for which labels are significant
    
    # recursive function to get current label positions
    # (note the unit is "npc" for all components of each distance)
    strip.npc <- function(dd){
      if(!"unit.arithmetic" %in% class(dd)) {
        return(as.numeric(dd))
      }
      
      d1 <- strip.npc(dd$arg1)
      d2 <- strip.npc(dd$arg2)
      fn <- dd$fname
      return(lazyeval::lazy_eval(paste(d1, fn, d2)))
    }
    
    full.range <- sapply(seq_along(d), function(i) strip.npc(d[i]))
    selected.range <- sapply(seq_along(d[d.select]), function(i) strip.npc(d[d.select][i]))
    
    return(unit(seq(from = max(selected.range) + k*(max(full.range) - max(selected.range)),
                    to = min(selected.range) - k*(min(selected.range) - min(full.range)), 
                    length.out = sum(d.select)), 
                "npc"))
  }
  new.y.positions <- repelled.y(new.label$y,
                                d.select = new.label$label != "")
  new.flag <- segmentsGrob(x0 = new.label$x,
                           x1 = new.label$x + unit(0.15, "npc"),
                           y0 = new.label$y[new.label$label != ""],
                           y1 = new.y.positions)
  
  # shift position for selected labels
  new.label$x <- new.label$x + unit(0.2, "npc")
  new.label$y[new.label$label != ""] <- new.y.positions
  
  # add flag to heatmap
  heatmap <- gtable::gtable_add_grob(x = heatmap,
                                     grobs = new.flag,
                                     t = 4, 
                                     l = 4
  )
  
  # replace label positions in heatmap
  heatmap$grobs[[which(heatmap$layout$name == "row_names")]] <- new.label
  
  # plot result
  grid.newpage()
  grid.draw(heatmap)
  
  # return a copy of the heatmap invisibly
  invisible(heatmap)
}



color_paper <- list(
  `rt-1` = "#ffc080",
  `rt-2` = "#c06000",
  `rtc4-1` = "#90bff9",
  `rtc4-2` = "#0000c0",
  `untreated-1` = "#a0a0a4"
)


library(EnhancedVolcano)
library(pheatmap)
library(monocle3)
library(grid)
library(zoo)
library(EnhancedVolcano)
library(progeny)
library(Seurat)
library(tidyr)
library(SeuratDisk)
library(dplyr)
require(MarkdownReports)


OutDir <- 'Zhang_2020_processed/'
dir.create(OutDir)

# CD40 dataset -----------------------------------------------------------------

load('Zhang_2020_data/CSF1R_expression.rda')
load('Zhang_2020_data/CSF1R_metadata.rda')

Zhang_Csf1r <- CreateSeuratObject(counts = CSF1R_expression_data$raw_counts, project = "Zhang_CSF1R",
                                 meta.data = CSF1R_metadata)



Zhang_Csf1r[['percent.mt']] <- PercentageFeatureSet(Zhang_Csf1r, pattern = "^mt\\-")
Zhang_Csf1r[['percent.ribo']] <- PercentageFeatureSet(Zhang_Csf1r, pattern = "^Rp")

Zhang_Csf1r <- NormalizeData(Zhang_Csf1r, normalization.method = "LogNormalize", 
                                        scale.factor = 10000)
Zhang_Csf1r <- ScaleData(Zhang_Csf1r)
Zhang_Csf1r <- FindVariableFeatures(Zhang_Csf1r)
Zhang_Csf1r <- RunPCA(Zhang_Csf1r)

Zhang_Csf1r <- RunUMAP(Zhang_Csf1r, dims = 1:50, 
                                  n.components = 2)
Zhang_Csf1r@misc$UMAP_recalc <- Zhang_Csf1r@reductions$umap
umaps <- as.matrix(
  Zhang_Csf1r@meta.data[,c("Global_UMAP_1", "Global_UMAP_2")])

colnames(umaps) <- c("UMAP_1", "UMAP_2")

Zhang_Csf1r@misc$UMAP_orig <- CreateDimReducObject(umaps)
Zhang_Csf1r@misc$UMAP_subcl <-
  CreateDimReducObject(as.matrix(Zhang_Csf1r@meta.data[,c("Sub_UMAP_1", "Sub_UMAP_2")]))
Zhang_Csf1r@reductions$umap <- Zhang_Csf1r@misc$UMAP_orig
Zhang_Csf1r@reductions$umap <- Zhang_Csf1r@misc$UMAP_subcl

DimPlot(Zhang_Csf1r, group.by = "Sub_Cluster")
FeaturePlot(Zhang_Csf1r, features = c(c("Cd8b1", "Il18r1", "Il18rap", 
                                                   "percent.mt", "percent.ribo")))

saveRDS(Zhang_Csf1r, 'Zhang_2020_processed/Zhang_Csf1r_processed_211117.Rds')


# Transfer to monocle for correlations -------

# extract genes
genes <- as.data.frame(rownames(Zhang_Csf1r@assays$RNA), 
                       row.names = rownames(Zhang_Csf1r@assays$RNA))
colnames(genes) <- "gene_short_name"

# extract cells
cells <- as.data.frame(
  Zhang_Csf1r@assays[["RNA"]]@data@Dimnames[[2]], 
  row.names = Zhang_Csf1r@assays[["RNA"]]@data@Dimnames[[2]])
colnames(cells) <- "barcode"

# extract expression matrix
expression_matrix <- Zhang_Csf1r@assays[["RNA"]]@data
expression_matrix <- expression_matrix[rownames(Zhang_Csf1r@assays$RNA), ]

# Assemble cell data set object
cds.Zhang_Csf1r <- new_cell_data_set(expression_matrix,
                                  cell_metadata = cells,
                                  gene_metadata = genes)

# transfer dataset information
cds.Zhang_Csf1r@colData <- cds.Zhang_Csf1r@colData %>%
  as.data.frame() %>%
  left_join(CSF1R_metadata, by = c("barcode" = "CellName")) %>%
  DataFrame(row.names = .$barcode)

# preprocess cds with variable features from seurat
cds.Zhang_Csf1r@int_colData@listData$reducedDims$PCA <- 
  Zhang_Csf1r@reductions[["pca"]]@cell.embeddings[row.names(cds.Zhang_Csf1r@colData),]

# perform umap dimensionality reduction

cds.Zhang_Csf1r@int_colData@listData$reducedDims$UMAP <- 
  Zhang_Csf1r@misc$UMAP_subcl@cell.embeddings[row.names(cds.Zhang_Csf1r@colData),]

plot_cells(cds.Zhang_Csf1r, color_cells_by = "Global_Cluster")

saveRDS(cds.Zhang_Csf1r, 'Zhang_2020_processed/Zhang_Csf1r_monocle_211117.Rds')
myelo_cells <- CSF1R_metadata %>% 
  filter(Global_Cluster == "Myeloid cell") %>% magrittr::use_series(CellName)

cds.Zhang_Csf1r_m <- cds.Zhang_Csf1r[,myelo_cells]
plot_cells(cds.Zhang_Csf1r_m, color_cells_by = "Global_Cluster")
saveRDS(cds.Zhang_Csf1r_m, 'Zhang_2020_processed/Zhang_Csf1r_monocle_Myelo_211117.Rds')
# calculate modules and perform correlation ------------------------------------

cds.MCells <- readRDS('Macro_DC_2/Comparison_Cell_Paper/cds_MCells.Rds')
# saveRDS(cds.MCells, 'Macro_DC_2/Comparison_Cell_Paper/cds_MCells.Rds')

labs <- c("TAM-1", "TAM-2", "TAM-3", "TAM-4", "TAM-5", "NC-Monocytes", "TAM-6", "cDC2", "TAM-7",
          "Monocytes", "cDC1-1", "TRM", "cDC1-2")
cds.MCells@colData$integrated_snn_res.0.5_named <- 
  factor(cds.MCells@colData$integrated_snn_res.0.5, labels = labs)

sc.obj.MCells$integrated_snn_res.0.5_named <- factor(
  sc.obj.MCells$integrated_snn_res.0.5, labels = labs)

plot <- DimPlot(sc.obj.MCells, group.by = "integrated_snn_res.0.5_named", 
        cols = alpha(c(pals::tol(12), "grey"), .5), pt.size = .75, combine = F)
plot <- plot[[1]]+ coord_equal(ratio = .8) +
  ggtitle("Myeloid cells", subtitle = "(27865 cells)") + theme( axis.line = element_blank(),
                                    title = element_text(hjust = 0.5),
                                    plot.subtitle = element_text(hjust = 0.5),
                                    axis.text = element_blank(),
                                    axis.ticks = element_blank(),
                                    axis.title = element_blank(),
                                    panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
                                    panel.grid.minor = element_blank(), 
                                    panel.grid.major = element_blank(),
                                    plot.background = element_rect(fill = "transparent",colour = NA))

renamed_vec <- c("untreated-1", "C12-iF", "C12-oF", "C12/aC4-iF", "C12/aC4-oF")
names(renamed_vec) <- c("untreated-1", "rt-1", "rt-2", "rtc4-1", "rtc4-2")

png('Macro_DC_2/UMAP_clusters_repelled_myeloidcells2.png', bg = "transparent", width = 600, height = 450)
LabelClusters(plot = plot, 'integrated_snn_res.0.5_named', 
              arrow = arrow(length = unit(1,"cm"), type = "closed", angle = 0),
              box.padding = 5, size =5,
              color = "black", max.overlaps = Inf)
dev.off()
plot <- DimPlot(sc.obj.MCells, group.by = "orig.ident",
        shuffle = TRUE, 
        cols = alpha(color_paper, .5), pt.size = .5, combine = FALSE)

png('Macro_DC_2/UMAP_myeloid_per_ident.png', bg = "transparent", width = 500, height = 400)
plot[[1]] + coord_equal(ratio = .8) +
  scale_color_manual(labels = as_labeller(renamed_vec), values=alpha(as.vector(color_paper), .5)) +
  ggtitle("Myeloid cells") + theme( axis.line = element_blank(),
                                    title = element_text(hjust = 0.5),
                                    plot.subtitle = element_text(hjust = 0.5),
                                    axis.text = element_blank(),
                                    axis.ticks = element_blank(),
                                    axis.title = element_blank(),
                                    panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
                                    panel.grid.minor = element_blank(), 
                                    panel.grid.major = element_blank(),
                                    plot.background = element_rect(fill = "transparent",colour = NA))
dev.off()
plot <- DimPlot(sc.obj.MCells, group.by = "integrated_snn_res.0.5_named",shuffle = T, 
                cols = alpha(c(pals::tol(12), "grey"), .5), pt.size = .75, combine = F)

png('Macro_DC_2/UMAP_myeloid_nolab.png', bg = "transparent", width = 600, height = 450)
plot[[1]] + coord_equal(ratio = .8) +
  ggtitle("Myeloid cells", subtitle = "(27865 cells)") + theme( axis.line = element_blank(),
                                                                title = element_text(hjust = 0.5),
                                                                plot.subtitle = element_text(hjust = 0.5),
                                                                axis.text = element_blank(),
                                                                axis.ticks = element_blank(),
                                                                axis.title = element_blank(),
                                                                panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
                                                                panel.grid.minor = element_blank(), 
                                                                panel.grid.major = element_blank(),
                                                                plot.background = element_rect(fill = "transparent",colour = NA))
dev.off()

cluster_cols <- c(pals::tol(12), "grey")
cluster_cols[-c(1,2,4,5,10)] = "grey"
plot <- DimPlot(sc.obj.integration, group.by = "integrated_snn_res.0.2",
        shuffle = F, label = F, pt.size = .1, combine = FALSE)
png('Macro_DC_2/Subclustering_selection_Macro2.png', bg = "transparent", width = 1000, height = 800)
plot[[1]] + coord_equal(ratio = .8) +
  scale_color_manual(values = alpha(cluster_cols, .5),guide = "none",) +
  ggtitle("") +
  theme(axis.line = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),
         axis.title = element_blank(), title = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA)) + NoLegend()
dev.off()


FeaturePlot(sc.obj.MCells, features = c("Cd209a", "Itgax", "Clec9a", "Ccl22"))


plot_cells(cds.MCells, color_cells_by = "integrated_snn_res.0.5")

cds.Zhang_Csf1r_m <- preprocess_cds(cds.Zhang_Csf1r_m, use_genes = VariableFeatures(Zhang_Csf1r))

modules_MC <- find_gene_modules(cds.Zhang_Csf1r_m, resolution = .05, 
                                reduction_method = "UMAP")


modules_MC$module %>% table()
modules_MC2$module %>% table()
modules_MC$id %>% unique() 
write.csv(modules_MC, 'Zhang_2020_processed/modules_Zhang_Myeloid.csv')

aggr_mod_ref <- aggregate_gene_expression(cds.Zhang_Csf1r_m, gene_group_df = modules_MC, 
                                          cell_group_df = cds.Zhang_Csf1r_m@colData[,c("barcode", "Sub_Cluster")])


aggr_mod_Mc <- aggregate_gene_expression(cds.MCells, gene_group_df = modules_MC, 
                                         cell_group_df = cds.MCells@colData[,c("barcode", "integrated_snn_res.0.5_named")])

cor_bound <- cor(cbind(aggr_mod_Mc, aggr_mod_ref) %>% as.matrix())
pl = corrplot::corrplot(cor_bound, order = 'hclust')
dev.off()
corrplot::corrplot(cor_bound[c(1:7, 9:10, 12, 21:28), 
                             c(1:7, 9:10, 12, 21:28)], order = 'hclust', addrect = 10)

pdf('Zhang_2020_processed/Correlation_Zhang_Myeloid.pdf', width = 10, height = 10)
corrplot::corrplot(cor_bound[c("Monocytes", "NC-Monocytes", "TAM-3", "TAM-4", "TAM-5",
                               "TAM-6", "TAM-7","TAM-1", "TAM-2", "TRM"), 
                             c(22,23,21,28,24,27,26,25)])
dev.off()
corrplot::corrplot(cor_bound[c("0", "1", "2", "3", "4", "5", "6", "8" ,"9", "11"), c(21:28)])

pdf('Zhang_2020_processed/Correlation_Zhang_DCs.pdf', width = 5, height = 5)
corrplot::corrplot(cor_bound[c("cDC2", "cDC1-1", "cDC1-2"), c(17:20)])
dev.off()


openxlsx::write.xlsx(list(`Modules Zhang et al.` = aggr_mod_ref, `Modules Hartmann et al.` = aggr_mod_Mc,
                          `Correlation to Zhang` = cor_bound), row.names = TRUE, overwrite = T,
                     file = 'Zhang_2020_processed/Modules_Cor_Zhang.xlsx')




# Similar correlation to Gubin et al. paper ------------------------------------

cds.ref.combined_MCells <- readRDS('Macro_DC_2/Comparison_Cell_Paper/cds.ref.combined_MCells.Rds')

modules_MC_gubin <- read.csv('Macro_DC_2/Comparison_Cell_Paper/modules_ref_Macro.csv')
modules_MC_gubin <- modules_MC_gubin[,-1]
labs <- c("1:Macro", "2:Macro", "3:Macro", "4:Macro", "8:NK", "7", "5:Macro",
          "9:CD8", "10:Mki-L", "11:Treg", "12:DC", "13:Neutro", "14", "16", "17:pDC",
          "18", "19")

cds.ref.combined_MCells@colData$GraphCluster_fix <-
  factor(cds.ref.combined_MCells@colData$GraphCluster, 
         labels=labs)
cells_keep <- cds.ref.combined_MCells@colData %>%
  as.data.frame() %>%
  filter(GraphCluster_fix %in% c("1:Macro", "2:Macro", "3:Macro", "4:Macro", "5:Macro")) %>%
  mutate(bc = row.names(.)) %>%
  magrittr::use_series(bc)

cds.ref.combined_MCells <- cds.ref.combined_MCells[,cells_keep]

cellgroupdf <- as.data.frame(cds.ref.combined_MCells@colData[,c("barcode", "GraphCluster_fix")])
cellgroupdf <- cellgroupdf %>% mutate(barcode = row.names(.))
aggr_mod_ref <- aggregate_gene_expression(cds.ref.combined_MCells, 
                                          gene_group_df = modules_MC_gubin, 
                                          cell_group_df = cellgroupdf)


aggr_mod_Mc <- aggregate_gene_expression(cds.MCells, gene_group_df = modules_MC_gubin, 
                                         cell_group_df = cds.MCells@colData[,c("barcode", 
                                                                               "integrated_snn_res.0.5_named")])

plot_cells(cds.ref.combined_MCells, color_cells_by  = 'GraphCluster_fix', reduction_method = 'tSNE',
                  cell_size = 1, group_label_size = 10)

cor_bound <- cor(cbind(aggr_mod_Mc, aggr_mod_ref) %>% as.matrix())
corrplot::corrplot(cor_bound, order = 'hclust')
pdf('Macro_DC_2/Comparison_Cell_Paper/correlation_Gubin_name_new.pdf', width = 10, height = 10)
corrplot::corrplot(cor_bound[c("Monocytes", "NC-Monocytes", "TAM-3", "TAM-4", "TAM-5",
                               "TAM-6", "TAM-7","TAM-1", "TAM-2", "TRM"),
                             c("1:Macro", "4:Macro","5:Macro","3:Macro",   "2:Macro" )])
dev.off()


openxlsx::write.xlsx(list(`Modules Gubin et al.` = aggr_mod_ref, `Modules Hartmann et al.` = aggr_mod_Mc,
                          `Correlation to Gubin` = cor_bound), row.names = TRUE, overwrite = T,
                     file = 'Macro_DC_2/Comparison_Cell_Paper/Modules_Cor_Gubin.xlsx')


# Calculate markers fuer res05 -------------------------------------------------
sc.obj.integration <- readRDS('processed.DS_2/processed_RTC4_integration_211115_scores.Rds')
all_sum <- sc.obj.integration@meta.data %>% dplyr::count(orig.ident)

sc.obj.MCells@active.ident <- sc.obj.MCells$integrated_snn_res.0.5_named
markers_res05 <- FindAllMarkers(sc.obj.MCells, assay = "RNA")

# plot heatmap

markers_res05 %>% filter(avg_log2FC>.6) %>%
  dplyr::group_by(cluster)  %>% 
  filter(p_val_adj<1e-50) %>% 
  ungroup()


goi_unique <- markers_res05 %>% filter(avg_log2FC>0.6) %>%
  dplyr::group_by(cluster)  %>% 
  filter(p_val_adj<1e-50) %>% 
  ungroup() %>%
  magrittr::use_series(gene)


goi_new <- markers_res05 %>% 
  dplyr::group_by(cluster)  %>% 
  filter(avg_log2FC>0.6, p_val_adj<1e-50) %>% 
  ungroup() %>%
  magrittr::use_series(gene) %>%
  unique()

avg_goi_M <- AverageExpression(sc.obj.MCells, features = goi_new, 
                               group.by = 'integrated_snn_res.0.5_named')


data.frame(real = 0:9, new = c(0,1,3,4,2,5,6,7,9, 8))

MCell_markers$order <- NA
cluster_order <- c("Monocytes", "NC-Monocytes","TAM-3", "TAM-4", "TAM-5", "TAM-6",
                   "TAM-2", "TAM-1", "TAM-7", "TRM", "cDC2", "cDC1-1", "cDC1-2")
order_vec <- 1:13

for (i in 1:13) {
  markers_res05$order[markers_res05$cluster==cluster_order[i]] <- order_vec[i]
}


gene_order <- markers_res05 %>% 
  dplyr::group_by(cluster)  %>% 
  filter(avg_log2FC>0.6, p_val_adj<1e-50) %>% 
  dplyr::arrange(p_val_adj) %>%
  ungroup() %>%
  select(gene,order) %>%
  filter(!duplicated(gene)==TRUE) %>%
  arrange(order) %>% magrittr::use_series(gene) %>% unique()


avg_goi_M_df <- avg_goi_M$RNA[gene_order,cluster_order]

vec <- (alpha(c(as.vector(pals::tol(12)), "grey"), .8))
names(vec) <- labs
annot_col <- list()
annot_col[['cluster']] = vec
ph1 <- pheatmap(avg_goi_M_df, 
                cluster_rows = F,
                # annotation_row = data.frame(row.names = goi_df$gene, cluster = as.factor(goi_df$cluster)),
                annotation_col = data.frame(row.names = labs, cluster = as.factor(labs)),
                color = rev(pals::brewer.rdbu(100)),
                annotation_colors = annot_col, main = "DEGs Macrophages and DCs",
                show_rownames = T, scale = 'row', cluster_cols = F,angle_col = 45, 
                width = 10, height = 10)


# annotate only selected genes
pl1 <- add.flag(ph1, c('Cd209a', 'Clec9a', 'Ccl22', 'Ly6c2', 'Nr4a1', 'Itgal',
                       'Mafb', 'Maf','Ccl12', 'Mgl2',
                       'Vegfa'), 
                repel.degree = .1)

pl2 <- DimPlot(sc.obj.MCells, group.by = 'integrated_snn_res.0.5_named',
               cols = (alpha(c(as.vector(pals::tol(12)),"grey"), .5)), label = TRUE, label.size = 10)
png('Macro_DC_2/UMAP_Heatmap_Temp_res05_MCells.png', width = 1500, height = 750)
cowplot::plot_grid(pl2, pl1,rel_widths = c(1,1))
dev.off()


openxlsx::write.xlsx(markers_res05 %>% 
                        dplyr::group_by(cluster)  %>% 
                        filter(avg_log2FC>0.6,  p_val_adj<1e-50) %>% 
                        ungroup(), 
                     file = 'Macro_DC_2/Markers_4Heatmap_Mcells_res05_211117.xlsx'
                     , asTable = T, overwrite = T, row.names = TRUE)

pl_monocytes <- markers_res05 %>% 
  dplyr::group_by(cluster)  %>% 
  filter(cluster == "Monocytes") %>% 
  ungroup()

pdf(paste0('Macro_DC_2/VolcanoPlots/MonocytesvsRest.pdf'), width = 20, height = 20)
EnhancedVolcano(pl_monocytes, 
                lab = pl_monocytes$gene,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75, title = "Monocytes_vs_rest", 
                subtitle = "All Conditions combined")

dev.off()

# Re plot scores per group -----------------------------

sc.obj.MCells_gsm <-
  readRDS('processed.DS_2/processed_RTC4_MCells_Gubin_integration_211115_scores.Rds')
# saveRDS(sc.obj.MCells_gsm, 'processed.DS_2/processed_RTC4_MCells_Gubin_integration_211115_scores.Rds')
labs <- c("TAM-1", "TAM-2", "TAM-3", "TAM-4", "TAM-5", "NC-Monocytes", "TAM-6", "cDC2", "TAM-7",
          "Monocytes", "cDC1-1", "TRM", "cDC1-2")

sc.obj.MCells_gsm$integrated_snn_res.0.5_named <- factor(
  sc.obj.MCells_gsm$integrated_snn_res.0.5, labels = labs)
sc.obj.MCells_gsm@meta.data$fixed_cluster <- 
  paste0(sc.obj.MCells_gsm@meta.data$integrated_snn_res.0.5_named,
         sc.obj.MCells_gsm@meta.data$GraphCluster_fix, recycle0 = T) %>%
  stringr::str_replace("NA", "")

scores_df_M <- sc.obj.MCells_gsm@meta.data %>%
  filter(fixed_cluster %in% c("TRM", "Monocytes", "TAM-1", "TAM-2",
                              "TAM-3", "TAM-4", "TAM-5", "TAM-6", "TAM-7",
                              "1:Macro", "2:Macro", "3:Macro",
                              "4:Macro", "5:Macro")) %>%
  dplyr::group_by(orig.ident) %>%
  summarise(
    id = unique(paste0(orig.ident, "_Macrophages")),
    GO.0070555_1 = median(GO.0070555_1),
    GO.0002548_1 = median(GO.0002548_1),
    GO.0071356_1 = median(GO.0071356_1),
    GO.0048247_1 = median(GO.0048247_1),
    GO.0071621_1 = median(GO.0071621_1),
    GO.0019221_1 = median(GO.0019221_1),
    GO.0070371_1 = median(GO.0070371_1),
    GO.1990869_1 = median(GO.1990869_1)
  )


renamed_vec <- c(
  GO.0070555_1 = "Response to IL-1",
  GO.0002548_1 = "Monocyte Chemotaxis",
  GO.0071356_1 = "cellular response to tumor necrosis factor",
  GO.0048247_1 = "lymphocyte chemotaxis",
  GO.0071621_1 = "granulocyte chemotaxis",
  GO.0019221_1 = "cytokine-mediated signaling pathway",
  GO.0070371_1 = "ERK1 and ERK2 cascade",
  GO.1990869_1 = "cellular response to chemokine"
)

scores_df_M2 <- scores_df_M %>% 
  rename_with(.cols = 3:10, .fn = function(x) x = renamed_vec) %>% t() %>%
  as.data.frame() %>% janitor::row_to_names(2, remove_rows_above = T)

annot <- data.frame(row.names = colnames(scores_df_M2),
                    sample = scores_df_M[,'orig.ident'],
                    source = c(rep("Gubin et al.", 4), rep("Hartmann et al.", 5)))
colnames(annot) <- c("sample", "source")


sample_col <- as.vector(c(pals::tol(4), as.named.vector(color_paper)))
names(sample_col) <- scores_df_M$orig.ident %>% unique()
source_col <- as.vector(pals::tol(2))
names(source_col) <- annot$source %>% unique()

annot_col <- list(sample = sample_col, source = source_col)

pdf('Macro_DC_2/DiffExp/GO_heatmaps/GO_heatmap_comparison_Gubin.pdf', width = 10, height = 4)
print(scores_df_M2 %>% 
        mutate(across(.fns = function(x) x = as.numeric(x))) %>%
        ungroup() %>% 
        t() %>% scale() %>% t() %>% clip.outliers() %>%
        pheatmap::pheatmap(color = pals::coolwarm(100), cluster_cols = F,
                           gaps_col = 4, 
                           annotation_col = annot, cluster_rows = F,
                           show_colnames = F, annotation_colors = annot_col))
dev.off()

df_gubin_comparison <- 
  scores_df_M2 %>% 
  mutate(across(.fns = function(x) x = as.numeric(x))) %>%
  ungroup() %>% 
  t() %>% scale() %>% t() %>% clip.outliers() %>%
  as.data.frame()



scores_df_M <- sc.obj.MCells@meta.data %>%
  dplyr::group_by(integrated_snn_res.0.5_named, orig.ident) %>%
  filter(integrated_snn_res.0.5_named %in% c("TRM", "Monocytes", "TAM-1", "TAM-2",
                                             "TAM-3", "TAM-4", "TAM-5", "TAM-6", "TAM-7")) %>%
  summarise(
    id = unique(paste0(orig.ident, "_", integrated_snn_res.0.5_named)),
    GO.0070555_1 = median(GO.0070555_1),
    GO.0002548_1 = median(GO.0002548_1),
    GO.0071356_1 = median(GO.0071356_1),
    GO.0048247_1 = median(GO.0048247_1),
    GO.0071621_1 = median(GO.0071621_1),
    GO.0019221_1 = median(GO.0019221_1),
    GO.0070371_1 = median(GO.0070371_1),
    GO.1990869_1 = median(GO.1990869_1)
  )


scores_df_M2 <- scores_df_M %>% 
  rename_with(.cols = 4:11, .fn = function(x) x = renamed_vec) %>% t() %>%
  as.data.frame() %>% janitor::row_to_names(3, remove_rows_above = T)


annot <- data.frame(row.names = colnames(scores_df_M2),
                    sample = scores_df_M$orig.ident,
                    celltype = scores_df_M$integrated_snn_res.0.5_named)
colnames(annot) <- c("sample", "cluster")


celltype_col <-c(as.vector(pals::tol(12)), "grey")
names(celltype_col) <- scores_df_M$integrated_snn_res.0.5_named %>% unique()
sample_col <- as.named.vector(color_paper)
annot_col <- list(cluster = celltype_col, sample = sample_col)

pdf('Macro_DC_2/DiffExp/GO_heatmaps/GO_heatmap_selection.pdf', width = 15, height = 4)
print(scores_df_M2 %>%
        mutate(across(.fns = function(x) x = as.numeric(x))) %>%
        ungroup() %>%
        t() %>% scale() %>% t() %>% clip.outliers() %>%
        pheatmap::pheatmap(color = pals::coolwarm(100), cluster_cols = F,
                           gaps_col = seq(0,45, by = 5), 
                           annotation_col = annot,
                           show_colnames = F, annotation_colors = annot_col))
dev.off()

df_selection <- scores_df_M2 %>%
  mutate(across(.fns = function(x) x = as.numeric(x))) %>%
  ungroup() %>%
  t() %>% scale() %>% t() %>% clip.outliers() %>%
  as.data.frame()

df_all <- scores_df_M2 %>%
  mutate(across(.fns = function(x) x = as.numeric(x))) %>%
  ungroup() %>%
  t() %>% scale() %>% t() %>% clip.outliers() %>%
  as.data.frame()


# Do on whole integrated dataset for comparison

scores_df_M <- sc.obj.integration@meta.data %>%
  dplyr::group_by(integrated_snn_res.0.2_named, orig.ident) %>%
  filter(integrated_snn_res.0.2_named %in% c("Cd8 mito-low", "Naive TC", "NK", "B", "Cd8 mito-high",
                                             "Treg", "Neutro", "pDCs")) %>%
  summarise(
    id = unique(paste0(orig.ident, "_", integrated_snn_res.0.2_named)),
    GO.0070555_1 = median(GO.0070555_1),
    GO.0002548_1 = median(GO.0002548_1),
    GO.0071356_1 = median(GO.0071356_1),
    GO.0048247_1 = median(GO.0048247_1),
    GO.0071621_1 = median(GO.0071621_1),
    GO.0019221_1 = median(GO.0019221_1),
    GO.0070371_1 = median(GO.0070371_1),
    GO.1990869_1 = median(GO.1990869_1)
  )


scores_df_M2 <- scores_df_M %>% 
  rename_with(.cols = 4:11, .fn = function(x) x = renamed_vec) %>% t() %>%
  as.data.frame() %>% janitor::row_to_names(3, remove_rows_above = T)


annot <- data.frame(row.names = colnames(scores_df_M2),
                    sample = scores_df_M$orig.ident,
                    celltype = scores_df_M$integrated_snn_res.0.2_named)
colnames(annot) <- c("sample", "cluster")


celltype_col <-c(as.vector(pals::tol(12)), "grey")[c(3,6:9,11:12)]
names(celltype_col) <- scores_df_M$integrated_snn_res.0.2_named %>% unique()
sample_col <- as.named.vector(color_paper)
annot_col <- list(cluster = celltype_col, sample = sample_col)

pdf('Macro_DC_2/DiffExp/GO_heatmaps/GO_heatmap_OtherCellTypes.pdf', width = 15, height = 4)
print(scores_df_M2 %>%
        mutate(across(.fns = function(x) x = as.numeric(x))) %>%
        ungroup() %>%
        t() %>% scale() %>% t() %>% clip.outliers() %>%
        as.data.frame() %>%
        pheatmap::pheatmap(color = pals::coolwarm(100), cluster_cols = F,
                           gaps_col = seq(0,30, by = 5), 
                           annotation_col = annot,
                           show_colnames = F, annotation_colors = annot_col))
dev.off()

df_othercelltypes <- scores_df_M2 %>%
  mutate(across(.fns = function(x) x = as.numeric(x))) %>%
  ungroup() %>%
  t() %>% scale() %>% t() %>% clip.outliers() 

openxlsx::write.xlsx(list(`All Clusters` = df_all, `Selection Clusters` = df_selection,
                          `Other Cell Types` = df_othercelltypes, 
                          `Comparison Gubin Macrophages` = df_gubin_comparison),
                     row.names = TRUE, overwrite = T,
                     file = 'Macro_DC_2/DiffExp/GO_heatmaps/Heatmap_SourceData_211118.xlsx')




# Stuff ========
color_paper2 = color_paper
names(color_paper2) <- c("C12-iF", "C12-oF", "C12/aC4-iF", "C12/aC4-oF",
                         "untreated-1")
p <- sc.obj.MCells@meta.data %>%
  group_by(orig.ident) %>%
  dplyr::count(integrated_snn_res.0.5_named) %>%
  left_join(all_sum, by = "orig.ident") %>%
  ungroup() %>%
  mutate(pct = round(n.x/n.y*100,2),
         orig.ident = factor(orig.ident),
         sample = factor(orig.ident, 
                         levels = c("untreated-1", "rt-1", "rt-2", "rtc4-1", "rtc4-2"),
                         labels = c("untreated-1", "C12-iF", "C12-oF", "C12/aC4-iF", 
                                    "C12/aC4-oF"))) %>%
  ggplot(aes(x = sample,  y = pct, fill = sample)) +
  geom_col() + 
  geom_text(aes(label = round(pct,1)), hjust = 0.5, position=position_dodge(width=0.9), 
            vjust=-0.25, size = 3) +
  facet_wrap(~integrated_snn_res.0.5_named, nrow = 1) +
  ylim(0,25) +
  scale_fill_manual(values = color_paper2, guide = "none") +
  theme_pubr(x.text.angle = 45)




p <- sc.obj.MCells@meta.data %>%
  group_by(orig.ident) %>%
  dplyr::count(integrated_snn_res.0.5_named) %>%
  left_join(all_sum, by = "orig.ident") %>%
  ungroup() %>%
  mutate(pct = round(n.x/n.y*100,2),
         orig.ident = factor(orig.ident),
         sample = factor(orig.ident, 
                         levels = c("untreated-1", "rt-1", "rt-2", "rtc4-1", "rtc4-2"),
                         labels = c("untreated-1", "C12-iF", "C12-oF", "C12/aC4-iF", 
                                    "C12/aC4-oF"))) %>%
  filter(integrated_snn_res.0.5_named %in% labs[samples_take]) %>%
  mutate(integrated_snn_res.0.5_named = factor(
    integrated_snn_res.0.5_named, levels = labs[samples_take][reorder_takes])) %>%
  ggplot(aes(x = sample,  y = pct, fill = sample)) +
  geom_col() + 
  geom_text(aes(label = round(pct,1)), hjust = 0.5, position=position_dodge(width=1), 
            vjust=-.75, size = 3) +
  facet_wrap(~integrated_snn_res.0.5_named, nrow = 2) +
  ylim(0,25) + ylab("Percentage of Dataset") + xlab(label = "") +
  scale_fill_manual(values = color_paper2, guide = "none") +
  theme_pubr(x.text.angle = 45) + theme(axis.text.x = element_text(size = 10))


samples_take <- c(1:7,9,10,12)
reorder_takes <- c(1:5,7,8,6,9,10)
labs[samples_take][reorder_takes]

g <- ggplot_gtable(ggplot_build(p))
strip_both <- which(grepl('strip-t', g$layout$name))
reorder_takes2 <- reorder_takes[c(6:10, 1:5)]
fills <- alpha(c(as.vector(pals::tol(12)), "grey"),.5)[samples_take][reorder_takes2]
k <- 1
for (i in strip_both) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g)
dev.off()

p <- sc.obj.MCells@meta.data %>%
  group_by(orig.ident) %>%
  dplyr::count(integrated_snn_res.0.5_named) %>%
  left_join(all_sum, by = "orig.ident") %>%
  ungroup() %>%
  mutate(pct = round(n.x/n.y*100,2),
         orig.ident = factor(orig.ident),
         sample = factor(orig.ident, 
                         levels = c("untreated-1", "rt-1", "rt-2", "rtc4-1", "rtc4-2"),
                         labels = c("untreated-1", "C12-iF", "C12-oF", "C12/aC4-iF", 
                                    "C12/aC4-oF"))) %>%
  filter(integrated_snn_res.0.5_named %in% labs[(1:13)[-samples_take]]) %>%
  ggplot(aes(x = sample,  y = pct, fill = sample)) +
  geom_col() + 
  geom_text(aes(label = round(pct,1)), hjust = 0.5, position=position_dodge(width=1), 
            vjust=-.75, size = 3) +
  facet_wrap(~integrated_snn_res.0.5_named, nrow = 1) +
  ylim(0,25) + ylab("Percentage of Dataset") + xlab(label = "") +
  scale_fill_manual(values = color_paper2, guide = "none") +
  theme_pubr(x.text.angle = 45) + theme(axis.text.x = element_text(size = 10))


g <- ggplot_gtable(ggplot_build(p))
strip_both <- which(grepl('strip-t', g$layout$name))

fills <- alpha(c(as.vector(pals::tol(12)), "grey"),.5)[(1:13)[-samples_take]]
k <- 1
for (i in strip_both) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g)
dev.off()







sc.obj.MCells@meta.data %>%
  filter(integrated_snn_res.0.5_named %in% c("TAM-2", "TAM-4", "TAM-5", "TAM-6")) %>%
  dplyr::count(orig.ident) %>%
  left_join(all_sum, by = "orig.ident") %>%
  ungroup() %>%
  mutate(pct = round(n.x/n.y*100,2),
         orig.ident = factor(orig.ident),
         sample = factor(orig.ident, 
                         levels = c("untreated-1", "rt-1", "rt-2", "rtc4-1", "rtc4-2"),
                         labels = c("untreated-1", "C12-iF", "C12-oF", "C12/aC4-iF", 
                                    "C12/aC4-oF"))) %>%
  ggplot(aes(x = sample,  y = pct, fill = sample)) +
  geom_col() + 
  geom_text(aes(label = round(pct,1)), hjust = 0.5, position=position_dodge(width=1), 
            vjust=-.75, size = 3) +
  ylim(0,50) + ylab("Percentage of Dataset") + xlab(label = "") +
  scale_fill_manual(values = color_paper2, guide = "none") +
  ggtitle("Lineage towards TAM-5", "Combination of TAM-2, TAM-4, \nTAM-6 and TAM-5") +
  theme_pubr(x.text.angle = 45) + theme(axis.text.x = element_text(size = 10))
cells <- row.names(sc.obj.MCells@meta.data %>%
                     filter(integrated_snn_res.0.5_named %in% c("TAM-2", "TAM-4", "TAM-5", "TAM-6")))

sc.obj.MCells@meta.data[cells,"TAM-5_lineage"] = sc.obj.MCells@meta.data[cells,"integrated_snn_res.0.5_named"]
sc.obj.MCells@meta.data$`TAM-5_lineage` <-
  factor(sc.obj.MCells@meta.data$`TAM-5_lineage`, labels = levels(sc.obj.MCells@meta.data$integrated_snn_res.0.5_named)[c(2,3,5,7)])



plot <- DimPlot(sc.obj.MCells, group.by = "TAM-5_lineage",split.by = "orig.ident", cols = fills[c(2,3,5,7)],
        na.value = alpha("grey", 0.2), combine = FALSE)

plot <- plot[[1]] + facet_wrap(~orig.ident, labeller = as_labeller(renamed_vec), ncol = 2) +
  ggtitle("Lineage towards TAM-5")+
  theme( axis.line = element_blank(), title = element_text(hjust = 0.5), 
         plot.subtitle = element_text(hjust = 0.5),
         axis.text = element_blank(),
         axis.ticks = element_blank(),
         axis.title = element_blank(),
         panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
         panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(),
         plot.background = element_rect(fill = "transparent",colour = NA))


png(filename = 'Macro_DC_2/UMAP_lineageTAM5.png', width = 750, height = 750, bg = "transparent")
plot
dev.off()

sc.obj.MCells <- readRDS('processed.DS_2/processed_RTC4_MCells_211115_scores.Rds')

genes <- c("Itgam", "Csf3r", "Siglech", "Cd209a", "Clec9a", "Ccl22", "Ly6c2", "Ly6g", "Chil3", 
           "Sell", "Nr4a1", "Itgal", "Ccl2", "Ccl9", "Apoe", "Fcgr4", "Pparg", "Fcer1g", "Irf8",
           "Xcr1", 
           "Ccl12", "Mgl2", "Vegfa", "H2-Ab1", "Mrc1", "Arg1", "Spp1", "Adgre1", 
           "Cxcl10", "Chil3", "Il1b", "Nos2")

plot_list <- FeaturePlot(sc.obj.MCells, features = genes,order = FALSE, pt.size = 1,
                         cols = c(alpha("grey", .5), "red"), max.cutoff = 'q99', combine = F) 

plot_list <- lapply(plot_list, function(x) x = x + coord_equal(ratio = 0.7) +
                      theme( axis.line = element_blank(), title = element_text(hjust = 0.5,
                                                                               size = 30), 
                             plot.subtitle = element_text(hjust = 0.5),
                             axis.text = element_blank(),
                             axis.ticks = element_blank(),
                             axis.title = element_blank(),
                             panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
                             panel.grid.minor = element_blank(), 
                             panel.grid.major = element_blank(),
                             plot.background = element_rect(fill = "transparent",colour = NA)))
for (i in 1:length(genes)) {
  png(paste0('Macro_DC_2/UMAP_genes/', genes[i], "_umap.png"), width = 600, height = 450,
      bg = "transparent")
  print(plot_list[[i]])
  dev.off()
}




