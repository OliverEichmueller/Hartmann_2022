# Transfer Ctrl to monocle
library(Seurat)
library(ggplot2)
library(ggpubr)

# Thu Nov 18 08:43:59 2021 ------------------------------


# Transfer to monocle for correlations -------
levs <- levels(sc.obj.integration$integrated_snn_res.0.2_named)

sc.obj.integration$integrated_snn_res.0.2_named <-
  factor(sc.obj.integration$integrated_snn_res.0.2_named,
         labels = c("TAM-1", "TAM-2", levs[3:13]))

# extract genes
genes <- as.data.frame(rownames(sc.obj.integration@assays$RNA), 
                       row.names = rownames(sc.obj.integration@assays$RNA))
colnames(genes) <- "gene_short_name"

# extract cells
cells <- as.data.frame(
  sc.obj.integration@assays[["RNA"]]@data@Dimnames[[2]], 
  row.names = sc.obj.integration@assays[["RNA"]]@data@Dimnames[[2]])
colnames(cells) <- "barcode"

# extract expression matrix
expression_matrix <- sc.obj.integration@assays[["RNA"]]@data
expression_matrix <- expression_matrix[rownames(sc.obj.integration@assays$RNA), ]

# Assemble cell data set object
cds.integration <- new_cell_data_set(expression_matrix,
                                     cell_metadata = cells,
                                     gene_metadata = genes)

# transfer dataset information
cds.integration@colData <- cds.integration@colData %>%
  as.data.frame() %>%
  left_join(sc.obj.integration@meta.data %>% mutate(barcode = row.names(.)), by = c("barcode")) %>%
  DataFrame(row.names = .$barcode)

# preprocess cds with variable features from seurat
cds.integration@int_colData@listData$reducedDims$PCA <- 
  sc.obj.integration@reductions[["pca"]]@cell.embeddings[row.names(cds.integration@colData),]

# perform umap dimensionality reduction

cds.integration@int_colData@listData$reducedDims$UMAP <- 
  sc.obj.integration@misc$reductions.backup$umap_2D@cell.embeddings[row.names(cds.integration@colData),]

cds.integration@preprocess_aux$gene_loadings <- sc.obj.integration@reductions$pca@feature.loadings

plot_cells(cds.integration, color_cells_by = "integrated_snn_res.0.2_named")

saveRDS(cds.integration, 'processed.DS_2/Integration_all_monocle_211117.Rds')

DC_cells1 <- cds.MCells@colData %>% as.data.frame() %>%
  filter(integrated_snn_res.0.5_named %in% c("cDC2", "cDC1-1", "cDC1-2")) %>%
  magrittr::use_series(barcode)
DC_cells2 <- cds.integration@colData %>% as.data.frame() %>%
  filter(integrated_snn_res.0.2_named %in% c("pDC")) %>%
  magrittr::use_series(barcode)
DC_cellsref <- cds.Zhang_Csf1r_m@colData %>% as.data.frame() %>%
  filter(Sub_Cluster %in% c("mM03_pDC-Siglech", "mM04_cDC2-Cd209a", "mM05_cDC2-Itgax",
                            "mM06_cDC1-Clec9a", "mM07_cDc1-Ccl22")) %>%
  magrittr::use_series(barcode)

cds.DCs <- combine_cds(list(cds.MCells[,DC_cells1], cds.integration[,DC_cells2]))
cds.DCs@colData$DC_names <- stringr::str_replace(paste0(cds.DCs@colData$integrated_snn_res.0.2_named, cds.DCs@colData$integrated_snn_res.0.5_named)
                     , "NA", "")

modules_MC <- read.csv('Zhang_2020_processed/modules_Zhang_Myeloid.csv')
modules_MC$module %>% table
modules_MC <- modules_MC[,-1]

aggr_mod_ref <- aggregate_gene_expression(cds.Zhang_Csf1r_m, gene_group_df = modules_MC, 
                                          cell_group_df = cds.Zhang_Csf1r_m@colData[DC_cellsref,c("barcode", "Sub_Cluster")])

cds.DCs@colData$barcode <- row.names(cds.DCs@colData)

aggr_mod_Mc <- aggregate_gene_expression(cds.DCs, gene_group_df = modules_MC, 
                                         cell_group_df = cds.DCs@colData[,c("barcode", "DC_names")])

cor_bound <- cor(cbind(aggr_mod_Mc, aggr_mod_ref) %>% as.matrix())
pl = corrplot::corrplot(cor_bound, order = 'hclust')
dev.off()


pdf('Zhang_2020_processed/Correlation_Zhang_DC.pdf', width = 5, height = 5)
corrplot::corrplot(cor_bound[c(1:4), 
                             c(5:8)])
dev.off()
saveRDS(cds.DCs, 'processed.DS_2/DCs_all_monocle_211117.Rds')


# plot integration all for main figure -----------------------------------------


labs <- c("TAM-1", "TAM-2", "CD8 T cells 1", "NC-Mono", "cDC1-1", "naive T cells",
          "NK cells", "B cells", "CD8 T cells 2", "cDC1-2", "T reg", "Neutrophils", "pDC")

sc.obj.integration$integrated_snn_res.0.2_named <-
  factor(sc.obj.integration$integrated_snn_res.0.2,
         labels = labs)

plot <- DimPlot(sc.obj.integration, group.by = "integrated_snn_res.0.2_named", 
                cols = alpha(c(pals::tol(12), "grey"), .5), pt.size = .75, combine = F)
plot <- plot[[1]]+ coord_equal(ratio = .8) +
  ggtitle("Tumor infiltrating cells", subtitle = "(37441 cells)") + 
  theme( axis.line = element_blank(), title = element_text(hjust = 0.5), 
         plot.subtitle = element_text(hjust = 0.5),
         axis.text = element_blank(),
         axis.ticks = element_blank(),
         axis.title = element_blank(),
         panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
         panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(),
         plot.background = element_rect(fill = "transparent",colour = NA))


png('processed.plots_2/UMAP_color_fix.png', bg = "transparent", width = 600, height = 450)
LabelClusters(plot = plot, 'integrated_snn_res.0.2_named', 
              arrow = arrow(length = unit(1,"cm"), type = "closed", angle = 0),
              box.padding = 5, size =5,
              color = "black", max.overlaps = Inf)
dev.off()
png('processed.plots_2/UMAP_color_fix_nolabels.png', bg = "transparent", width = 800, height = 600)
plot
dev.off()




p <- sc.obj.integration@meta.data %>%
  group_by(orig.ident) %>%
  dplyr::count(integrated_snn_res.0.2_named) %>%
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
  facet_wrap(~integrated_snn_res.0.2_named, nrow = 1) +
  ylim(0,50) + ylab("Percentage of Dataset") + xlab(label = "") +
  scale_fill_manual(values = color_paper2, guide = "none") +
  theme_pubr(x.text.angle = 45) + theme(axis.text.x = element_text(size = 10))


g <- ggplot_gtable(ggplot_build(p))
strip_both <- which(grepl('strip-t', g$layout$name))

fills <- alpha(c(as.vector(pals::tol(12)), "grey"),.5)
k <- 1
for (i in strip_both) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g)
dev.off()

sc.obj.integration <- readRDS('processed.DS_2/processed_RTC4_integration_211115_scores.Rds')

genes <- unique(c("Itgam", "Csf3r", "Siglech", "Cd209a", "Clec9a", "Ccl22", "Ly6c2", "Nr4a1",
           "Ccl12", "Mgl2", "Vegfa", "H2-Ab1", "Mrc1", "Klrb1c", "Cd3e", "Cd19", "Foxp3",
           "Cd8b1", "Cd4", "Cxcl10","Csf1r","Itgam","Arg1",
                             "Apoe","Mrc1",
                             "Havcr2","Lcn2","Mmp9",
                             "Ccl22","Il12b","Xcr1","Cd8b1",
                             "Ifng","Cd8a", "Il18r1",
                             "Il2ra","Foxp3","S1pr1", "Tcf7","Gzma","Klrb1c",
                             "Pax5","Cd22","Cd19","Siglech",
                             "Cd209d","Tex2"))

plot_list <- FeaturePlot(sc.obj.integration, features = unique(genes),order = FALSE,
            cols = c(alpha("grey", .5), "red"), max.cutoff = 'q99', combine = F) 
  
plot_list <- lapply(plot_list, function(x) x = x + coord_equal(ratio = 0.8) +
  theme( axis.line = element_blank(), title = element_text(size = 30, hjust = 0.5), 
         plot.subtitle = element_text(hjust = 0.5),
         axis.text = element_blank(),
         axis.ticks = element_blank(),
         axis.title = element_blank(),
         panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
         panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(),
         plot.background = element_rect(fill = "transparent",colour = NA)))
for (i in 1:length(genes)) {
  png(paste0('processed.plots_2/UMAP_genes/', genes[i], "_umap.png"), width = 600, height = 450,
      bg = "transparent")
  print(plot_list[[i]])
  dev.off()
}




labs <- c("Myeloid", "Myeloid", "Lymphoid", "Myeloid", "Myeloid", "Lymphoid",
          "Lymphoid", "B cells", "Lymphoid", "Myeloid", "Lymphoid", "Neutrophils", "pDC")


sc.obj.integration$integrated_snn_res.0.2_grouped <- factor(
  sc.obj.integration$integrated_snn_res.0.2, labels = labs)

DimPlot(sc.obj.integration, group.by = "integrated_snn_res.0.2_grouped")


pl1 <- sc.obj.integration@meta.data %>%
  dplyr::group_by(orig.ident) %>%
  dplyr::count(integrated_snn_res.0.2_grouped) %>%
  dplyr::group_by(orig.ident, integrated_snn_res.0.2_grouped) %>%
  summarise(sum_group = sum(n)) %>%
  mutate(rel = round(sum_group/sum(sum_group)*100,2)) %>%
  arrange(desc(integrated_snn_res.0.2_grouped)) %>%
  mutate(prop = rel / sum(rel) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop) %>%
  ggplot(aes(x="", y = rel, fill = integrated_snn_res.0.2_grouped)) +
  geom_bar(stat="identity", width = 1, color = "white")+
  coord_polar("y", start = 0, direction = 1)+
  scale_fill_manual(values = alpha(c(pals::brewer.set1(3)[-3], pals::tol(12)[c(8,12)], "grey"), 0.75)
                    ,guide = "none", aesthetics = c("color", "fill"))+
  geom_text_repel(aes(label = rel, y=ypos, color = integrated_snn_res.0.2_grouped), min.segment.length = .2,
                  nudge_x = .8, box.padding = 0, alpha = 1)+
  facet_grid(cols = vars(orig.ident), labeller = as_labeller(renamed_vec)) + 
  theme_void()  + theme(strip.text = element_text(size = 15))

pl2 <- DimPlot(sc.obj.integration, group.by = "integrated_snn_res.0.2_grouped", shuffle = TRUE,
               cols = alpha(c(pals::brewer.set1(3)[-3], pals::tol(12)[c(8,12)], "grey"), 0.25),
               label = F, repel = TRUE) + ggtitle("Broad cell groups") + coord_equal(ratio = 0.8) +
  NoAxes()

plot_grid(pl2, pl1, ncol = 1, rel_heights = c(1,0.3))

png('processed.plots_2/UMAP_broad_clusters.png', width = 600, height = 500)
pl2
dev.off()

png('processed.plots_2/UMAP_broad_clusters_nolab.png', width = 600, height = 500)
pl2
dev.off()


pdf('processed.plots_2/pie_broad_clusters.pdf', width = 10, height = 5)
pl1
dev.off()
