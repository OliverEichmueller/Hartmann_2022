# Plots for T-cell figure
## Sun Nov 21 13:14:16 2021 ------------------------------


color_paper <- list(
  `rt-1` = "#ffc080",
  `rt-2` = "#c06000",
  `rtc4-1` = "#90bff9",
  `rtc4-2` = "#0000c0",
  `untreated-1` = "#a0a0a4"
)
library(biomaRt)
library(ggpubr)
library(CodeAndRoll2)
library(monocle3)
library(ComplexHeatmap)
library(EnhancedVolcano)
library(Seurat)
library(UCell)
library(dplyr)
library(ggplot2)

sc.obj.integration <- readRDS('processed.DS_2/processed_RTC4_integration_211115_scores.Rds')
all_sum <- sc.obj.integration@meta.data %>% dplyr::count(orig.ident)

sc.obj.TCells <- readRDS('processed.DS_2/processed_RTC4_TCells_integration_scores_211111.Rds')


OutDir <- 'TCells_new/Plots_final/'
dir.create(OutDir)


cluster_cols <- c(pals::tol(12), "grey")
cluster_cols[-c(3,6,7,9,11)] = "grey"
plot <- DimPlot(sc.obj.integration, group.by = "integrated_snn_res.0.2",
                shuffle = F, label = F, pt.size = .1, combine = FALSE)
png('TCells_new/Plots_final/Subclustering_selection_TCells.png', bg = "transparent", width = 1000, height = 800)
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
labs <- c("Dividing Eff", "NK", "Cd8-EM", "naive Cd4, Cd8", "Treg", "Cd8 Act/Exh", "IL-18-R. Cd8", 
          "naive Cd8", "Th1", "Dividing IL-18-R. Cd8")
sc.obj.TCells$integrated_snn_res.0.25_Final <- factor(
  sc.obj.TCells$integrated_snn_res.0.25_named, labels = labs
)


plot <- DimPlot(sc.obj.TCells, group.by = "integrated_snn_res.0.25_Final", 
                cols = alpha(c(pals::tol(10), "grey"), .5), pt.size = .75, combine = F)
plot <- plot[[1]]+ coord_equal(ratio = .8) +
  ggtitle("Lymphoid cells", subtitle = "(8006 cells)") + theme( axis.line = element_blank(),
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

png('TCells_new/Plots_final/UMAP_clusters_repelled_TCells.png', bg = "transparent", width = 600, height = 450)
LabelClusters(plot = plot, 'integrated_snn_res.0.25_Final', 
              arrow = arrow(length = unit(1,"cm"), type = "closed", angle = 0),
              box.padding = 5, size =5,
              color = "black", max.overlaps = Inf)
dev.off()
plot <- DimPlot(sc.obj.TCells, group.by = "orig.ident",
                shuffle = TRUE, 
                cols = alpha(color_paper, .5), pt.size = .5, combine = FALSE)

png('TCells_new/Plots_final/UMAP_Lymphoid_per_ident.png', bg = "transparent", width = 500, height = 400)
plot[[1]] + coord_equal(ratio = .8) +
  scale_color_manual(labels = as_labeller(renamed_vec), values=alpha(as.vector(color_paper), .5)) +
  ggtitle("Lymphoid cells") + theme( axis.line = element_blank(),
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
plot <- DimPlot(sc.obj.TCells, group.by = "integrated_snn_res.0.25_Final",shuffle = T, 
                cols = alpha(c(pals::tol(10)), .5), pt.size = .75, combine = F)

png('TCells_new/Plots_final/UMAP_lymphoid_nolab.png', bg = "transparent", width = 600, height = 450)
plot[[1]] + coord_equal(ratio = .8) +
  ggtitle("Lymphoid cells", subtitle = "(8006 cells)") + theme( axis.line = element_blank(),
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



# Plot %s ---------------------

color_paper2 <- color_paper[c("untreated-1", "rt-1", "rt-2", "rtc4-1", "rtc4-2")]
names(color_paper2) <- c("untreated-1", "C12-iF", "C12-oF", "C12/aC4-iF", 
                         "C12/aC4-oF")

p <- sc.obj.TCells@meta.data %>%
  group_by(orig.ident) %>%
  dplyr::count(integrated_snn_res.0.25_Final) %>%
  left_join(all_sum, by = "orig.ident") %>%
  ungroup() %>%
  mutate(pct = round(n.x/n.y*100,2),
         orig.ident = factor(orig.ident),
         sample = factor(orig.ident, 
                         levels = c("untreated-1", "rt-1", "rt-2", "rtc4-1", "rtc4-2"),
                         labels = c("untreated-1", "C12-iF", "C12-oF", "C12/aC4-iF", 
                                    "C12/aC4-oF")),
         integrated_snn_res.0.25_Final = 
           factor(integrated_snn_res.0.25_Final,
                  levels = c( "NK", "Th1","Treg",  "naive Cd4, Cd8","naive Cd8", 
                              "Dividing Eff",  "Cd8 Act/Exh", "Cd8-EM",  "Dividing IL-18-R. Cd8", "IL-18-R. Cd8")
                                                  )) %>%
  ggplot(aes(x = sample,  y = pct, fill = sample)) +
  geom_col() + 
  geom_text(aes(label = round(pct,1)), hjust = 0.5, position=position_dodge(width=1), 
            vjust=-.75, size = 3) +
  facet_wrap(~integrated_snn_res.0.25_Final, nrow = 2) +
  ylim(0,12) + ylab("Percentage of Dataset") + xlab(label = "") +
  scale_fill_manual(values = color_paper2, guide = "none") +
  theme_pubr(x.text.angle = 45) + theme(axis.text.x = element_text(size = 10))


g <- ggplot_gtable(ggplot_build(p))
strip_both <- which(grepl('strip-t', g$layout$name))
data.frame(labs)
fills <- alpha(c(as.vector(pals::tol(10))),.5)[c(1, 6, 3, 10, 7, 2, 9, 5, 4, 8)]
k <- 1
for (i in strip_both) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g)
dev.off()



# Check GO-Terms again ---------

ensembl <- useEnsembl("ensembl",  dataset = "mmusculus_gene_ensembl")


genes <- getBM(attributes = c('external_gene_name'), #  'ensembl_transcript_id', 'go_id'
               filters = "go_parent_term",  uniqueRows = TRUE,
               values = "GO:0050863", mart = ensembl)[,1]

sc.obj.TCells <- AddModuleScore(sc.obj.TCells, 
                                features = list(GO.0071346 = intersect(row.names(sc.obj.TCells), genes)),
                                name = "GO.0050863_")

sc.obj.TCells@active.ident <- sc.obj.TCells$integrated_snn_res.0.25_Final

VlnPlot(sc.obj.TCells, features = "GO.0050863_1", split.by = "orig.ident", pt.size = 0)


sc.obj.TCells@meta.data %>%
  filter(orig.ident %in% c("rt-2", "rtc4-2")) %>%
  filter(integrated_snn_res.0.25_Final %in% c("naive Cd8", "naive Cd4, Cd8")) %>%
  dplyr::group_by(integrated_snn_res.0.25_Final, orig.ident) %>%
  summarise(
    GO.0002697_1 = median(GO.0002697_1),
    # GO.0002683_1 = median(GO.0002683_1),
    GO.0030098_1 = median(GO.0030098_1),
    GO.0001819_1 = median(GO.0001819_1),
    GO.0001906_1 = median(GO.0001906_1),
    GO.0050863_1 = median(GO.0050863_1),
    GO.0002286_1 = median(GO.0002286_1)
  ) %>%
  tidyr::gather("Score", "Expression", -integrated_snn_res.0.25_Final, -orig.ident) %>%
  group_by(Score) %>% 
  mutate(Expression = scale(Expression)) %>%
  ggplot(aes(x = orig.ident, y = integrated_snn_res.0.25_Final, fill = Expression)) +
  geom_tile() + scale_fill_viridis_c() + facet_wrap(~Score) +
  theme_pubr(x.text.angle = 45)

# GO.0002697_1 = regulation of immune effector process
# GO.0002683_1 = negative regulation of immune system process
# GO.0030098_1 = lymphocyte differentiation
# GO.0002286_1 = T cell activation involved in immune response
# GO.0001819_1 = positive regulation of cytokine production
# GO.0001906_1 = cell killing
# GO.0050863_1 = regulation of T cell activation



scores_df_T <- sc.obj.TCells@meta.data %>%
  filter(orig.ident %in% c("rt-2", "rtc4-2")) %>%
  filter(integrated_snn_res.0.25_Final %in% c("naive Cd8", "naive Cd4, Cd8")) %>%
  dplyr::group_by(integrated_snn_res.0.25_Final, orig.ident) %>%
  summarise(
    cluster = unique(paste(integrated_snn_res.0.25_Final, orig.ident, sep = "_")),
    GO.0002697_1 = median(GO.0002697_1),
    GO.0030098_1 = median(GO.0030098_1),
    GO.0001819_1 = median(GO.0001819_1),
    GO.0001906_1 = median(GO.0001906_1),
    GO.0050863_1 = median(GO.0050863_1),
    GO.0002286_1 = median(GO.0002286_1)
  )


renamed_vec <- c(
  GO.0002697_1 = "regulation of immune effector process",
  GO.0030098_1 = "lymphocyte differentiation",
  GO.0002286_1 = "T cell activation involved in immune response",
  GO.0001819_1 = "positive regulation of cytokine production",
  GO.0001906_1 = "cell killing",
  GO.0050863_1 = "regulation of T cell activation"
)

scores_df_T2 <- scores_df_T %>% 
  rename_with(.cols = 4:9, .fn = function(x) x = renamed_vec) %>% t() %>%
  as.data.frame() %>% janitor::row_to_names(3, remove_rows_above = T)

annot <- data.frame(row.names = colnames(scores_df_T2),
                    sample = rep(c("C12-oF", "C12/aC4-oF"), 2),
                    CellType = scores_df_T[,1])
colnames(annot) <- c("sample", "CellType")


sample_col <- as.vector(as.named.vector(color_paper)[c(2,4)])
names(sample_col) <- c("C12-oF", "C12/aC4-oF")
CellType_col <- rep(as.vector(pals::tol(10))[c(4,8)])
names(CellType_col) <- annot$CellType %>% unique()

annot_col <- list(sample = sample_col, CellType = CellType_col)

pdf('TCells_new/Plots_final/GSEA_naiveTCells.pdf', width = 10, height = 4)
print(scores_df_T2 %>% 
        mutate(across(.fns = function(x) x = as.numeric(x))) %>%
        ungroup() %>% 
        t() %>% scale() %>% t() %>% clip.outliers() %>%
        pheatmap::pheatmap(color = pals::coolwarm(100), cluster_cols = F,
                           gaps_col = 4, 
                           annotation_col = annot, cluster_rows = F,
                           show_colnames = F, annotation_colors = annot_col))
dev.off()

# plot genes --------


genes <- c("Cd4", "Cd8a", "Cd8b1", "Klrb1c", "Mki67", "Il18r1", "Il18rap",
           "Foxp3", "Lag3", "Pdcd1", "Gzma", "Havcr2", "S1pr1", "Bhlhe40", "Ctla4", 
           "Cd3g", "Nkg7", "Cxcr6", "Ifng", "Lef1", "Il6ra", "Sell", "Il7r", "Tnfrsf4",
           "Itgb8", "Rgs16", "Litaf", "Ccr7", "Tnfsf8", "Tnfsf11")

plot_list <- FeaturePlot(sc.obj.TCells, features = genes,order = FALSE, pt.size = 1,
                         cols = c(alpha("grey", .5), "red"), max.cutoff = 'q99', combine = F) 

plot_list <- lapply(plot_list, function(x) x = x + coord_equal(ratio = 0.8) +
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
  png(paste0('TCells_new/Plots_final/UMAP_genes/TCells_', genes[i], "_umap.png"), width = 600, height = 450,
      bg = "transparent")
  print(plot_list[[i]])
  dev.off()
}

genes2 <- c("Gzma", "Cd8a")
plot_list <- FeaturePlot(sc.obj.TCells, features = genes2,order = FALSE, 
                         split.by = "orig.ident",
            cols = c(alpha("grey", .5), "red"), max.cutoff = 'q99', combine = F, pt.size = 1)

plot_list[[1]] + NoLegend()

sc.obj.TCells$Sample <- factor(sc.obj.TCells$orig.ident, labels = 
                                 c("C12-iF", "C12-oF", "C12/aC4-iF", "C12/aC4-oF", "untreated-1"))

png('TCells_new/Plots_final/UMAP_genes/Gzma_Cd8a_perident.png', width = 1000, height = 400)
FeaturePlot(sc.obj.TCells, features = genes2,order = FALSE, 
            split.by = "Sample",
            cols = c(alpha("grey", .5), "red"), max.cutoff = 'q99', combine = TRUE, pt.size = 1) &
  NoAxes()
dev.off()
FeaturePlot(sc.obj.TCells, features = genes2,order = FALSE, 
            split.by = "Sample",
            cols = c(alpha("grey", .5), "red"), max.cutoff = 'q99', combine = TRUE, pt.size = 1) +
  NoAxes()


# Plot Scores

genes <- colnames(sc.obj.TCells@meta.data)[c(30,36:52)]

plot_list <- FeaturePlot(sc.obj.TCells, features = genes,order = FALSE, pt.size = 1,
                         cols = c(alpha("grey", .5), "red"), max.cutoff = 'q99', combine = F) 

plot_list <- lapply(plot_list, function(x) x = x + coord_equal(ratio = 0.8) +
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
  png(paste0('TCells_new/Plots_final/UMAP_genes/Scores_TCells_', genes[i], "_umap.png"), width = 600, height = 450,
      bg = "transparent")
  print(plot_list[[i]])
  dev.off()
}

sc.obj.TCells <- AddModuleScore_UCell(sc.obj.TCells, features = list(NK = "Klrb1c"),
                                      name = "_Score")

plot <- VlnPlot(sc.obj.TCells, features = "CD8_UCell", pt.size = 0, cols = alpha(pals::tol(10),.75), 
        combine = F)

pdf('TCells_new/Plots_final/Scores/CD8_Score_violin.pdf')
plot[[1]] + geom_hline(yintercept = 0.5) + ggtitle("CD8 Score", subtitle = "Expression of Cd8a and Cd8b1") +
  xlab("") + theme(plot.subtitle = element_text(hjust = 0.5))
dev.off()

plot <- VlnPlot(sc.obj.TCells, features = "NK_Score", pt.size = 0, cols = alpha(pals::tol(10),.75), 
        combine = F)

pdf('TCells_new/Plots_final/Scores/NK_Score_violin.pdf')
plot[[1]] + geom_hline(yintercept = 0.25) + ggtitle("NK Score", subtitle = "Expression of Klrb1c (Nk1.1)") +
  xlab("") + theme(plot.subtitle = element_text(hjust = 0.5))
dev.off()

plot <- VlnPlot(sc.obj.TCells, features = "CD8_LAG3_Score", pt.size = 0, cols = alpha(pals::tol(10),.75), 
        combine = F)

pdf('TCells_new/Plots_final/Scores/CD8_LAG3_Score_violin.pdf')
plot[[1]] + geom_hline(yintercept = 0.66) + ggtitle("CD8 + LAG3 Score", 
                                                    subtitle = "Expression of Cd8a, Cd8b1 and Lag3") +
  xlab("") + theme(plot.subtitle = element_text(hjust = 0.5))
dev.off()

saveRDS(sc.obj.TCells, 'processed.DS_2/processed_RTC4_TCells_integration_scores_211121.Rds')


# Pie Chart Tregs --------------------------------------------------------------

pl <- sc.obj.TCells@meta.data %>%
  dplyr::group_by(orig.ident) %>%
  dplyr::count(integrated_snn_res.0.25_Final) %>%
  mutate(CellType = ifelse(integrated_snn_res.0.25_Final == "Treg", "Treg", "Lymphoid"),
         CellType2 = ifelse(integrated_snn_res.0.25_Final == "Treg", "Treg", 
                            ifelse(integrated_snn_res.0.25_Final == "NK","NK", "T cells"))) %>%
  dplyr::group_by(orig.ident, CellType) %>%
  summarise(sum_group = sum(n)) %>%
  mutate(rel = round(sum_group/sum(sum_group)*100,2)) %>%
  ggplot(aes(x="", y = rel, fill = CellType)) +
  geom_bar(stat="identity", width = 1)+
  coord_polar("y", start = 0)+
  scale_fill_manual(values = c("grey", pals::tol(10)[5]))+
  geom_text_repel(aes(label = rel, y=rep(c(65,5),5)), min.segment.length = .2,
                  nudge_x = 0.7)+
  facet_grid(cols = vars(orig.ident), labeller = as_labeller(renamed_vec)) + 
  theme_void()  + theme(strip.text = element_text(size = 15))

pdf('TCells_new/Plots_final/perc_Tregs_Lymphoid.pdf', width = 7.5, height = 3)
pl
dev.off()

pl <- sc.obj.TCells@meta.data %>%
  dplyr::group_by(orig.ident) %>%
  dplyr::count(integrated_snn_res.0.25_Final) %>%
  mutate(CellType = ifelse(integrated_snn_res.0.25_Final == "Treg", "Treg", "Lymphoid"),
         CellType2 = ifelse(integrated_snn_res.0.25_Final == "Treg", "Treg", 
                            ifelse(integrated_snn_res.0.25_Final == "NK","NK", "T cells"))) %>%
  dplyr::group_by(orig.ident, CellType2) %>%
  summarise(sum_group = sum(n)) %>%
  mutate(rel = round(sum_group/sum(sum_group)*100,2)) %>%
  ggplot(aes(x="", y = rel, fill = CellType2)) +
  geom_bar(stat="identity", width = 1)+
  coord_polar("y", start = 0)+
  scale_fill_manual(values = c(pals::tol(10)[2], "grey", pals::tol(10)[5]))+
  geom_text_repel(aes(label = rel, y=rep(c(95,45,5),5)), min.segment.length = .2,
                  nudge_x = 0.7)+
  facet_grid(cols = vars(orig.ident), labeller = as_labeller(renamed_vec)) + 
  theme_void()  + theme(strip.text = element_text(size = 15))

pdf('TCells_new/Plots_final/perc_Tregs_TCells_NK.pdf', width = 7.5, height = 3)
pl
dev.off()

pl <- sc.obj.TCells@meta.data %>%
  dplyr::group_by(orig.ident) %>%
  dplyr::count(integrated_snn_res.0.25_Final) %>%
  mutate(CellType = ifelse(integrated_snn_res.0.25_Final == "Treg", "Treg", "Lymphoid"),
         CellType2 = ifelse(integrated_snn_res.0.25_Final == "Treg", "Treg", 
                            ifelse(integrated_snn_res.0.25_Final == "NK","NK", "T cells"))) %>%
  filter(CellType2 != "NK") %>%
  dplyr::group_by(orig.ident, CellType2) %>%
  summarise(sum_group = sum(n)) %>%
  mutate(rel = round(sum_group/sum(sum_group)*100,2)) %>%
  ggplot(aes(x="", y = rel, fill = CellType2)) +
  geom_bar(stat="identity", width = 1)+
  coord_polar("y", start = 0)+
  scale_fill_manual(values = c("grey", pals::tol(10)[5]))+
  geom_text_repel(aes(label = rel, y=rep(c(45,5),5)), min.segment.length = .2,
                  nudge_x = 0.7)+
  facet_grid(cols = vars(orig.ident), labeller = as_labeller(renamed_vec)) + 
  theme_void()  + theme(strip.text = element_text(size = 15))

pdf('TCells_new/Plots_final/perc_Tregs_TCells_noNK.pdf', width = 7.5, height = 3)
pl
dev.off()



# Re-Plot Enhanced Volcano for Tumor 2 differences -----------------------------

cluster5_degs <- openxlsx::read.xlsx('TCells_new/diff_exp/VolcanoPlots/Tumor2_rtvscombi/DE_per_cluster_t2.xlsx',
                    sheet = "cluster_5")

cluster7_degs <- openxlsx::read.xlsx('TCells_new/diff_exp/VolcanoPlots/Tumor2_rtvscombi/DE_per_cluster_t2.xlsx',
                    sheet = "cluster_7")
pdf('TCells_new/Plots_final/DEG_stuff/Cluster5_CD8ActExh_volcano.pdf', width = 10, height = 10)
EnhancedVolcano(cluster2_degs, 
                lab = cluster2_degs$gene,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75, title = paste0("CD8 Act./Exh"), 
                subtitle = "C12 vs. C12/aC4 in out-of-field tumors")
dev.off()

pdf('TCells_new/Plots_final/DEG_stuff/Cluster7_naiveCD8_volcano.pdf', width = 10, height = 10)
EnhancedVolcano(cluster7_degs, 
                lab = cluster7_degs$gene,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75, title = paste0("naive CD8"), 
                subtitle = "C12 vs. C12/aC4 in out-of-field tumors")
dev.off()


library(enrichplot)

cluster7 <- clusterProfiler::enrichGO(cluster7_degs %>% filter(avg_log2FC< -1, p_val_adj<0.05) %>%
                                       dplyr::select(gene, avg_log2FC, p_val_adj) %>% magrittr::use_series(gene),
                                     OrgDb = 'org.Mm.eg.db', keyType = "SYMBOL", ont = "BP")

cluster7_pw <- pairwise_termsim(cluster7)
dotplot(cluster7, showCategory = 60)
pdf('TCells_new/Plots_final/DEG_stuff/Cluster7_naiveCD8_tree.pdf', height = 20, width = 10)
treeplot(cluster7_pw, showCategory = 40)
dev.off()

# Plot GO heatmap
genes <- getBM(attributes = c('external_gene_name'), #  'ensembl_transcript_id', 'go_id'
               filters = "go_parent_term",  uniqueRows = TRUE,
               values = "GO:0001819", mart = ensembl)[,1]

sc.obj.TCells <- AddModuleScore(sc.obj.TCells, 
                                features = list(GO.0071346 = intersect(row.names(sc.obj.TCells), genes)),
                                name = "GO.0001819_")

VlnPlot(sc.obj.TCells, features = "GO.0001819_1", split.by = "orig.ident", pt.size = 0)
# GO.0002697_1 = regulation of immune effector process
# GO.0002683_1 = negative regulation of immune system process
# GO.0001819_1 = positive regulation of cytokine production

# GO.0001906_1 = cell killing
# GO.0050863_1 = regulation of T cell activation
# GO.0002286_1 = T cell activation involved in immune response
# GO.0030098_1 = lymphocyte differentiation
# GO.0070098_1 = chemokine-mediated signaling pathway
# GO.0060326_1 = cell chemotaxis
# GO.0070231_1 = T cell apoptotic process
# GO.0070233_1 = negative regulation of T cell apoptotic process
# GO.0071356_1 = cellular response to tumor necrosis factor
# GO.0001819_1 = positive regulation of cytokine production



scores_df_T <- sc.obj.TCells@meta.data %>%
  filter(orig.ident %in% c("rt-2", "rtc4-2")) %>%
  filter(integrated_snn_res.0.25_Final %in% c("naive Cd8", "naive Cd4, Cd8")) %>%
  dplyr::group_by(integrated_snn_res.0.25_Final, orig.ident) %>%
  summarise(
    cluster = unique(paste(integrated_snn_res.0.25_Final, orig.ident, sep = "_")),
    GO.0001906_1 = median(GO.0001906_1),
    GO.0050863_1 = median(GO.0050863_1),
    GO.0002286_1 = median(GO.0002286_1),
    GO.0030098_1 = median(GO.0030098_1),
    GO.0070098_1 = median(GO.0070098_1),
    GO.0060326_1 = median(GO.0060326_1),
    GO.0070231_1 = median(GO.0070231_1),
    GO.0070233_1 = median(GO.0070233_1),
    GO.0071356_1 = median(GO.0071356_1),
    GO.0001819_1 = median(GO.0001819_1)
  )

sc.obj.TCells@meta.data %>% 
  dplyr::select( GO.0001906_1, GO.0050863_1, GO.0002286_1, GO.0030098_1, 
                 GO.0070098_1, GO.0060326_1, 
                 GO.0070231_1, GO.0070233_1, GO.0071356_1, GO.0001819_1, orig.ident, integrated_snn_res.0.25_Final) %>%
  filter(orig.ident %in% c("rt-2", "rtc4-2")) %>%
  filter(integrated_snn_res.0.25_Final %in% c("naive Cd8", "naive Cd4, Cd8")) %>%
  tidyr::gather(key = "GO", value = "Expression", -orig.ident, -integrated_snn_res.0.25_Final) %>%
  ggplot(aes(x = orig.ident, y = Expression, fill = orig.ident)) +
  geom_boxplot() + scale_fill_manual(values = color_paper) +
  facet_grid(cols = vars(GO), rows = vars(integrated_snn_res.0.25_Final))




renamed_vec <- c(
  GO.0001906_1 = "cell killing",
  GO.0050863_1 = "regulation of T cell activation",
  GO.0002286_1 = "T cell activation involved in immune response",
  GO.0030098_1 = "lymphocyte differentiation",
  GO.0070098_1 = "chemokine-mediated signaling pathway",
  GO.0060326_1 = "cell chemotaxis",
  GO.0070231_1 = "T cell apoptotic process",
  GO.0070233_1 = "negative regulation of T cell apoptotic process",
  GO.0071356_1 = "cellular response to tumor necrosis factor",
  GO.0001819_1 = "positive regulation of cytokine production"
)

scores_df_T2 <- scores_df_T %>% 
  rename_with(.cols = 4:13, .fn = function(x) x = renamed_vec) %>% t() %>%
  as.data.frame() %>% janitor::row_to_names(3, remove_rows_above = T)

annot <- data.frame(row.names = colnames(scores_df_T2),
                    sample = rep(c("C12-oF", "C12/aC4-oF"), 2),
                    CellType = scores_df_T[,1])
colnames(annot) <- c("sample", "CellType")


sample_col <- as.vector(as.named.vector(color_paper)[c(2,4)])
names(sample_col) <- c("C12-oF", "C12/aC4-oF")
CellType_col <- rep(as.vector(pals::tol(10))[c(4,8)])
names(CellType_col) <- annot$CellType %>% unique()

annot_col <- list(sample = sample_col, CellType = CellType_col)

pdf('TCells_new/Plots_final/GSEA_naiveTCells_new.pdf', width = 10, height = 4)
print(scores_df_T2 %>% 
        mutate(across(.fns = function(x) x = as.numeric(x))) %>%
        ungroup() %>% 
        t() %>% scale() %>% t() %>% clip.outliers() %>%
        pheatmap::pheatmap(color = pals::coolwarm(100), cluster_cols = F,
                           gaps_col = 4, 
                           annotation_col = annot, cluster_rows = F,
                           show_colnames = F, annotation_colors = annot_col))
dev.off()



# Re-plot All-cells integration with name correction


# DC_1 = Ccl22
# DC_2 = Clec9a

# Summary Plot MCells ----------------------------------------------------------

labs <- c("TAM-1", "TAM-2", "TAM-3", "TAM-4", "TAM-5", "NC-Monocytes", "TAM-6", "cDC2", "TAM-7",
          "Monocytes", "cDC1-1", "TRM", "cDC1-2")

sc.obj.MCells$integrated_snn_res.0.5_named <- factor(
  sc.obj.MCells$integrated_snn_res.0.5, labels = labs)

labs <- c("TAM-Group-1", "TAM-Group-2", "TAM-Group-1", "TAM-Group-2", "TAM-Group-2", 
          "NC-Monocytes", "TAM-Group-2", "DC", "TAM-Group-3",
          "Monocytes", "DC", "TRM", "DC")

sc.obj.MCells$integrated_snn_res.0.5_grouped <- factor(
  sc.obj.MCells$integrated_snn_res.0.5, labels = labs)

pl1 <- sc.obj.MCells@meta.data %>%
  dplyr::group_by(orig.ident) %>%
  filter(integrated_snn_res.0.5_grouped != "DC") %>%
  dplyr::count(integrated_snn_res.0.5_grouped) %>%
  dplyr::group_by(orig.ident, integrated_snn_res.0.5_grouped) %>%
  summarise(sum_group = sum(n)) %>%
  mutate(rel = round(sum_group/sum(sum_group)*100,2)) %>%
  arrange(desc(integrated_snn_res.0.5_grouped)) %>%
  mutate(prop = rel / sum(rel) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop) %>%
  ggplot(aes(x="", y = rel, fill = integrated_snn_res.0.5_grouped)) +
  geom_bar(stat="identity", width = 1, color = "white")+
  coord_polar("y", start = 0, direction = 1)+
  scale_fill_manual(values = alpha(c(pals::brewer.set1(5)[-4], pals::tol(12)[c(10,12)]), 0.75),
                    guide = "none")+
  geom_text_repel(aes(label = rel, y=ypos), min.segment.length = .2,
                  nudge_x = .7, box.padding = 0)+
  facet_grid(cols = vars(orig.ident), labeller = as_labeller(renamed_vec)) + 
  theme_void()  + theme(strip.text = element_text(size = 15))

pl2 <- DimPlot(sc.obj.MCells, group.by = "integrated_snn_res.0.5_grouped",
        cols = alpha(c(pals::brewer.set1(5)[1:3], "grey",
                       pals::brewer.set1(5)[5],
                       pals::tol(12)[c(10,12)]), 0.5),
        label = F, repel = TRUE) + ggtitle("Myeloid cell groups") + coord_equal(ratio = 0.8) +
  NoAxes()

plot_grid(pl2, pl1, ncol = 1, rel_heights = c(1,0.3))

png('Macro_DC_2/UMAP_clusters_Myeloid_groups.png', width = 600, height = 500)
pl2
dev.off()

png('Macro_DC_2/UMAP_clusters_Myeloid_groups_nolab.png', width = 600, height = 500)
pl2
dev.off()

pdf('Macro_DC_2/Pie_clusters_Myeloid_groups_nolab.pdf', width = 10, height = 5)
pl1
dev.off()




# Re-plot All-cells integration with name correction ---------------------------


# DC_1 = Ccl22
# DC_2 = Clec9a

labs <- c("TAM-1", "TAM-2", "CD8 T cells 1", "NC-Mono", "cDC1-1", "naive T cells",
          "NK cells", "B cells", "CD8 T cells 2", "cDC1-2", "T reg", "Neutrophils", "pDC")

sc.obj.integration$integrated_snn_res.0.2_named <-
  factor(sc.obj.integration$integrated_snn_res.0.2,
         labels = labs)

DimPlot(sc.obj.integration, group.by = "integrated_snn_res.0.2_named", label = T,
        )


saveRDS(sc.obj.integration, 'processed.DS_2/processed_RTC4_integration_211115_scores.Rds')
saveRDS(sc.obj.MCells, 'processed.DS_2/processed_RTC4_MCells_211115_scores.Rds')
saveRDS(sc.obj.TCells, 'processed.DS_2/processed_RTC4_TCells_integration_scores_211121.Rds')

  