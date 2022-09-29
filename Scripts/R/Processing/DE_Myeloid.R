## Diff Exp of cell types M Cells
##------ Mon Nov  1 14:58:07 2021 ------##

## Oliver Eichmueller

# load required libraries ------------------------------------------------------
library(pheatmap)
library(EnhancedVolcano)
library(progeny)
library(Seurat)
library(tidyr)
library(SeuratDisk)
library(dplyr)
require(MarkdownReports)

source('~/R/Git_Projects/RTC4_2021/Scripts/Parameters.RTC4.processing.R')
source('~/Documents/GitHub/Seurat.utils/Functions/Plotting.dim.reduction.3D.R')
source('~/Documents/GitHub/Seurat.utils/00.Load.Seurat.Utils.LOCAL.R')
source("https://raw.githubusercontent.com/vertesy/CodeAndRoll/master/CodeAndRoll.R")

proc.Dir <- "~/R/Git_Projects/RTC4_2021/processed.DS_2/"



sc.obj.MCells <- 
  readRDS(file = paste0(proc.Dir, "processed_RTC4_MCells_integration_scores_211029.Rds"))



# DE Testing for tumor 1, 2 and untreated --------------------------------------

sc.obj.MCells$Radiation <- factor(sc.obj.MCells$orig.ident, labels = c("radiated.T1", "unrad.T2", "radiated.T1", "unrad.T2", "unrad.Ctrl"))

DimPlot(sc.obj.MCells, group.by = "Radiation")



clusters <- levels(sc.obj.MCells$integrated_snn_res.0.3)

M_cellsDE_per_cluster_rad <- list()

for (i in 1:length(clusters)) {
  sc.obj = subset(sc.obj.MCells, integrated_snn_res.0.3==clusters[i])
  sc.obj@active.ident = factor(sc.obj$Radiation)
  M_cellsDE_per_cluster_rad[[paste0("cluster_", clusters[i])]] <-
    FindAllMarkers(sc.obj, assay = "RNA", only.pos = TRUE)
}

M_cellsDE_per_cluster_rad_filtered <- 
  lapply(M_cellsDE_per_cluster_rad, 
         function(x) x = x %>% filter(p_val_adj <0.05) %>%
           dplyr::select(p_val_adj, cluster, gene) %>%
           pivot_wider(names_from = "cluster", values_from = "p_val_adj"))


M_cellsDE_per_cluster_rad_filtered <- plyr::ldply(M_cellsDE_per_cluster_rad_filtered, data.frame)
M_cellsDE_per_cluster_rad_filtered %>%
  filter(.id == "cluster_1")

sc.obj.MCells@active.ident <- sc.obj.MCells$integrated_snn_res.0.3
VlnPlot(sc.obj.MCells, features = "Gzmc", group.by = 'integrated_snn_res.0.3',
        split.by = "orig.ident", slot = 'data', assay = 'RNA')


FeaturePlot(sc.obj.MCells, features = 'Gzmc', split.by = "orig.ident")


# new DE just comparing Rad with no rad ----------------------------------------
clusters <- levels(sc.obj.MCells$integrated_snn_res.0.3)

M_cellsDE_per_cluster_radvsunrad <- list()

for (i in 1:length(clusters)) {
  sc.obj = subset(sc.obj.MCells, integrated_snn_res.0.3==clusters[i])
  sc.obj@active.ident = factor(sc.obj$Radiation)
  M_cellsDE_per_cluster_radvsunrad[[paste0("cluster_", clusters[i])]] <-
    FindMarkers(sc.obj,ident.1 = "radiated.T1", ident.2 = "unrad.T2", assay = "RNA")
}

M_cellsDE_per_cluster_radvsunrad_filtered <- 
  lapply(M_cellsDE_per_cluster_radvsunrad, 
         function(x) x = x %>%
           dplyr::select(p_val_adj, avg_log2FC))

clusters <- names(M_cellsDE_per_cluster_radvsunrad_filtered)
labs <- c("M0", "M2_1?", "M1_1", "MDSCs?", "M2_2?", "NC-Mono", "DC_1", 
          "DC_2", "DC_3", "DC_4")

i=10
png(paste0('Macro_DC_2/VolcanoPlots/T1vsT2/', labs[i], "_",
           clusters[i], "_volcano.png"), width = 800, height = 800)
EnhancedVolcano(M_cellsDE_per_cluster_radvsunrad_filtered[[clusters[i]]], 
                lab = rownames(M_cellsDE_per_cluster_radvsunrad_filtered[[clusters[i]]]),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75, title = paste0(clusters[i], " ", labs[i]), 
                subtitle = "T1 vs. T2 (irresp. of CTLA4)")

dev.off()

openxlsx::write.xlsx(lapply(M_cellsDE_per_cluster_radvsunrad_filtered,
                            function(x) x = x%>%
                              dplyr::mutate(gene =row.names(.))), file = 'Macro_DC_2/VolcanoPlots/T1vsT2/M_cellsDE_per_cluster_radvsunrad_filtered.xlsx'
                     , asTable = T, overwrite = T, row.names = TRUE)

M_cellsDE_ALL_radvsunrad_filtered <- 
  FindMarkers(sc.obj.MCells,ident.1 = "radiated.T1", ident.2 = "unrad.T2", assay = "RNA"
              , group.by = "Radiation")

png(paste0('Macro_DC_2/VolcanoPlots/T1vsT2/', "All_MacroDC_volcano.png"), width = 800, height = 800)
EnhancedVolcano(M_cellsDE_ALL_radvsunrad_filtered, 
                lab = rownames(M_cellsDE_ALL_radvsunrad_filtered),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75, title = "All Macro and DC", 
                subtitle = "T1 vs. T2 (irresp. of CTLA4)")

dev.off()

# Compare CTLA4 ----------------------------------------------------------------
sc.obj.MCells$CTLA4 <- factor(sc.obj.MCells$orig.ident, labels = c("no", "no", "yes", "yes", "Ctrl"))



clusters <- levels(sc.obj.MCells$integrated_snn_res.0.3)

MCells_DE_per_cluster_CTLA4vsno <- list()

for (i in 1:length(clusters)) {
  sc.obj = subset(sc.obj.MCells, integrated_snn_res.0.3==clusters[i])
  sc.obj@active.ident = sc.obj$CTLA4
  MCells_DE_per_cluster_CTLA4vsno[[paste0("cluster_", clusters[i])]] <-
    FindMarkers(sc.obj,ident.1 = "yes", ident.2 = "no", assay = "RNA")
}

MCells_DE_per_cluster_CTLA4vsno_filtered <- 
  lapply(MCells_DE_per_cluster_CTLA4vsno, 
         function(x) x = x %>%
           dplyr::select(p_val_adj, avg_log2FC))

clusters <- names(MCells_DE_per_cluster_CTLA4vsno_filtered)
labs <- c("M0", "M2_1?", "M1_1", "MDSCs?", "M2_2?", "NC-Mono", "DC_1", 
          "DC_2", "DC_3", "DC_4")

i=10
png(paste0('Macro_DC_2/VolcanoPlots/CTLA4vsno/', labs[i], "_",
           clusters[i], "_volcano.png"), width = 800, height = 800)
EnhancedVolcano(MCells_DE_per_cluster_CTLA4vsno_filtered[[clusters[i]]], 
                lab = rownames(MCells_DE_per_cluster_CTLA4vsno_filtered[[clusters[i]]]),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75, title = paste0(clusters[i], " ", labs[i]), 
                subtitle = "CTLA4 vs. no CTLA4 (irresp. of RT)")

dev.off()

openxlsx::write.xlsx(lapply(MCells_DE_per_cluster_CTLA4vsno_filtered,
                            function(x) x = x%>%
                              dplyr::mutate(gene =row.names(.))), file = 'Macro_DC_2/VolcanoPlots/CTLA4vsno/MCells_DE_per_cluster_CTLA4vsno_filtered.xlsx'
                     , asTable = T, overwrite = T, row.names = TRUE)

M_cellsDE_ALL_CTLA4vsno_filtered <- 
  FindMarkers(sc.obj.MCells,ident.1 = "yes", ident.2 = "no", assay = "RNA"
              , group.by = "CTLA4")

png(paste0('Macro_DC_2/VolcanoPlots/CTLA4vsno/', "All_MacroDC_volcano.png"), width = 800, height = 800)
EnhancedVolcano(M_cellsDE_ALL_CTLA4vsno_filtered, 
                lab = rownames(M_cellsDE_ALL_CTLA4vsno_filtered),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75, title = "All Macro and DC", 
                subtitle = "CTLA4 vs. no CTLA4 (irresp. of RT)")

dev.off()

openxlsx::write.xlsx(M_cellsDE_ALL_CTLA4vsno_filtered%>%
                       dplyr::mutate(gene =row.names(.)), file = 'Macro_DC_2/VolcanoPlots/CTLA4vsno/MCells_DE_All_ctla4.xlsx'
                     , asTable = T, overwrite = T, row.names = TRUE)


# compare just within tumor 1 --------------------------------------------------

clusters <- levels(sc.obj.MCells$integrated_snn_res.0.3)

MCells_DE_per_cluster_t1 <- list()

for (i in 1:length(clusters)) {
  sc.obj = subset(sc.obj.MCells, integrated_snn_res.0.3==clusters[i])
  sc.obj@active.ident = factor(sc.obj$orig.ident)
  MCells_DE_per_cluster_t1[[paste0("cluster_", clusters[i])]] <-
    FindMarkers(sc.obj,ident.1 = "rt-1", ident.2 = "rtc4-1", assay = "RNA")
}



clusters <- names(MCells_DE_per_cluster_t1)
labs <- c("M0", "M2_1?", "M1_1", "MDSCs?", "M2_2?", "NC-Mono", "DC_1", 
          "DC_2", "DC_3", "DC_4")

i=10
png(paste0('Macro_DC_2/VolcanoPlots/Tumor1_RTvsCombi/', labs[i], "_",
           clusters[i], "_volcano.png"), width = 800, height = 800)
EnhancedVolcano(MCells_DE_per_cluster_t1[[clusters[i]]], 
                lab = rownames(MCells_DE_per_cluster_t1[[clusters[i]]]),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75, title = paste0(clusters[i], " ", labs[i]), 
                subtitle = "RT vs. no RT+CTLA4 in Tumor 1")

dev.off()

openxlsx::write.xlsx(lapply(MCells_DE_per_cluster_t1,
                            function(x) x = x%>%
                              dplyr::mutate(gene =row.names(.))), file = 'Macro_DC_2/VolcanoPlots/Tumor1_RTvsCombi/MCells_DE_per_cluster_t1.xlsx'
                     , asTable = T, overwrite = T, row.names = TRUE)

MCells_DE_All_t1 <- 
  FindMarkers(sc.obj.MCells,ident.1 = "rt-1", ident.2 = "rtc4-1", assay = "RNA"
              , group.by = "orig.ident")

png(paste0('Macro_DC_2/VolcanoPlots/Tumor1_RTvsCombi/', "All_MacroDC_volcano.png"), width = 800, height = 800)
EnhancedVolcano(MCells_DE_All_t1, 
                lab = rownames(MCells_DE_All_t1),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75, title = "All Macro and DC", 
                subtitle = "RT vs. no RT+CTLA4 in Tumor 1")

dev.off()

openxlsx::write.xlsx(MCells_DE_All_t1%>%
                              dplyr::mutate(gene =row.names(.)), file = 'Macro_DC_2/VolcanoPlots/Tumor1_RTvsCombi/MCells_DE_All_t1.xlsx'
                     , asTable = T, overwrite = T, row.names = TRUE)


# compare tumor 1 to ctrl ------------------------------------------------------

clusters <- levels(sc.obj.MCells$integrated_snn_res.0.3)

MCells_DE_per_cluster_t1vsctrl <- list()

for (i in 1:length(clusters)) {
  sc.obj = subset(sc.obj.MCells, integrated_snn_res.0.3==clusters[i])
  sc.obj@active.ident = sc.obj$Radiation
  MCells_DE_per_cluster_t1vsctrl[[paste0("cluster_", clusters[i])]] <-
    FindMarkers(sc.obj,ident.1 = "radiated.T1", ident.2 = "unrad.Ctrl", assay = "RNA")
}



clusters <- names(MCells_DE_per_cluster_t1vsctrl)
labs <- c("M0", "M2_1?", "M1_1", "MDSCs?", "M2_2?", "NC-Mono", "DC_1", 
          "DC_2", "DC_3", "DC_4")

i=10
png(paste0('Macro_DC_2/VolcanoPlots/Tumor1_vs_Ctrl/', labs[i], "_",
           clusters[i], "_volcano.png"), width = 800, height = 800)
EnhancedVolcano(MCells_DE_per_cluster_t1vsctrl[[clusters[i]]], 
                lab = rownames(MCells_DE_per_cluster_t1vsctrl[[clusters[i]]]),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75, title = paste0(clusters[i], " ", labs[i]), 
                subtitle = "RT vs. Ctrl in Tumor 1")

dev.off()

openxlsx::write.xlsx(lapply(MCells_DE_per_cluster_t1vsctrl,
                            function(x) x = x%>%
                              dplyr::mutate(gene =row.names(.))), file = 'Macro_DC_2/VolcanoPlots/Tumor1_vs_Ctrl/MCells_DE_per_cluster_t1vsctrl.xlsx'
                     , asTable = T, overwrite = T, row.names = TRUE)

MCells_DE_All_t1vsctrl <- 
  FindMarkers(sc.obj.MCells,ident.1 = "radiated.T1", ident.2 = "unrad.Ctrl", assay = "RNA"
              , group.by = "Radiation")

png(paste0('Macro_DC_2/VolcanoPlots/Tumor1_vs_Ctrl/', "All_MacroDC_volcano.png"), width = 800, height = 800)
EnhancedVolcano(MCells_DE_All_t1vsctrl, 
                lab = rownames(MCells_DE_All_t1vsctrl),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75, title = "All Macro and DC", 
                subtitle = "RT vs. Ctrl in Tumor 1")

dev.off()

openxlsx::write.xlsx(MCells_DE_All_t1vsctrl%>%
                              dplyr::mutate(gene =row.names(.)), 
                     file = 'Macro_DC_2/VolcanoPlots/Tumor1_vs_Ctrl/MCells_DE_All_t1vsctrl.xlsx'
                     , asTable = T, overwrite = T, row.names = TRUE)



# compare tumor 2 to ctrl ------------------------------------------------------

clusters <- levels(sc.obj.MCells$integrated_snn_res.0.3)

MCells_DE_per_cluster_t2vsctrl <- list()

for (i in 1:length(clusters)) {
  sc.obj = subset(sc.obj.MCells, integrated_snn_res.0.3==clusters[i])
  sc.obj@active.ident = sc.obj$Radiation
  MCells_DE_per_cluster_t2vsctrl[[paste0("cluster_", clusters[i])]] <-
    FindMarkers(sc.obj,ident.1 = "unrad.T2", ident.2 = "unrad.Ctrl", assay = "RNA")
}



clusters <- names(MCells_DE_per_cluster_t2vsctrl)
labs <- c("M0", "M2_1?", "M1_1", "MDSCs?", "M2_2?", "NC-Mono", "DC_1", 
          "DC_2", "DC_3", "DC_4")

i=10
png(paste0('Macro_DC_2/VolcanoPlots/Tumor2_vs_Ctrl/', labs[i], "_",
           clusters[i], "_volcano.png"), width = 800, height = 800)
EnhancedVolcano(MCells_DE_per_cluster_t2vsctrl[[clusters[i]]], 
                lab = rownames(MCells_DE_per_cluster_t2vsctrl[[clusters[i]]]),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75, title = paste0(clusters[i], " ", labs[i]), 
                subtitle = "RT vs. Ctrl in Tumor 2")

dev.off()

openxlsx::write.xlsx(lapply(MCells_DE_per_cluster_t2vsctrl,
                            function(x) x = x%>%
                              dplyr::mutate(gene =row.names(.))), 
                     file = 'Macro_DC_2/VolcanoPlots/Tumor2_vs_Ctrl/MCells_DE_per_cluster_t2vsctrl.xlsx'
                     , asTable = T, overwrite = T, row.names = TRUE)

MCells_DE_All_t2vsctrl <- 
  FindMarkers(sc.obj.MCells,ident.1 = "unrad.T2", ident.2 = "unrad.Ctrl", assay = "RNA"
              , group.by = "Radiation")

png(paste0('Macro_DC_2/VolcanoPlots/Tumor2_vs_Ctrl/', "All_MacroDC_volcano.png"), width = 800, height = 800)
EnhancedVolcano(MCells_DE_All_t2vsctrl, 
                lab = rownames(MCells_DE_All_t2vsctrl),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75, title = "All Macro and DC", 
                subtitle = "RT vs. Ctrl in Tumor 2")

dev.off()

openxlsx::write.xlsx(MCells_DE_All_t2vsctrl%>%
                              dplyr::mutate(gene =row.names(.)), 
                     file = 'Macro_DC_2/VolcanoPlots/Tumor2_vs_Ctrl/MCells_DE_All_t2vsctrl.xlsx'
                     , asTable = T, overwrite = T, row.names = TRUE)


# compare just within tumor 2 --------------------------------------------------

clusters <- levels(sc.obj.MCells$integrated_snn_res.0.3)

MCells_DE_per_cluster_t2 <- list()

for (i in 1:length(clusters)) {
  sc.obj = subset(sc.obj.MCells, integrated_snn_res.0.3==clusters[i])
  sc.obj@active.ident = factor(sc.obj$orig.ident)
  MCells_DE_per_cluster_t2[[paste0("cluster_", clusters[i])]] <-
    FindMarkers(sc.obj,ident.1 = "rt-2", ident.2 = "rtc4-2", assay = "RNA")
}


clusters <- names(MCells_DE_per_cluster_t2)
labs <- c("M0", "M2_1?", "M1_1", "MDSCs?", "M2_2?", "NC-Mono", "DC_1", 
          "DC_2", "DC_3", "DC_4")

i=10
png(paste0('Macro_DC_2/VolcanoPlots/Tumor2_RTvsCombi/', labs[i], "_",
           clusters[i], "_volcano.png"), width = 800, height = 800)
EnhancedVolcano(MCells_DE_per_cluster_t2[[clusters[i]]], 
                lab = rownames(MCells_DE_per_cluster_t2[[clusters[i]]]),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75, title = paste0(clusters[i], " ", labs[i]), 
                subtitle = "RT vs. no RT+CTLA4 in Tumor 2")

dev.off()

openxlsx::write.xlsx(lapply(MCells_DE_per_cluster_t2,
                            function(x) x = x%>%
                              dplyr::mutate(gene =row.names(.))), file = 'Macro_DC_2/VolcanoPlots/Tumor2_RTvsCombi/MCells_DE_per_cluster_t2.xlsx'
                     , asTable = T, overwrite = T, row.names = TRUE)

MCells_DE_All_t2 <- 
  FindMarkers(sc.obj.MCells,ident.1 = "rt-2", ident.2 = "rtc4-2", assay = "RNA"
              , group.by = "orig.ident")

png(paste0('Macro_DC_2/VolcanoPlots/Tumor2_RTvsCombi/', "All_MacroDC_volcano.png"), width = 800, height = 800)
EnhancedVolcano(MCells_DE_All_t2, 
                lab = rownames(MCells_DE_All_t2),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75, title = "All Macro and DC", 
                subtitle = "RT vs. no RT+CTLA4 in Tumor 2")

dev.off()


openxlsx::write.xlsx(MCells_DE_All_t2%>%
                       dplyr::mutate(gene =row.names(.)), file = 'Macro_DC_2/VolcanoPlots/Tumor2_RTvsCombi/MCells_DE_All_t2.xlsx'
                     , asTable = T, overwrite = T, row.names = TRUE)
FeaturePlot(sc.obj.TCells, features = c('Il18r1'), split.by = 'orig.ident')

VlnPlot(sc.obj.TCells, features = "Tbx21", split.by = "orig.ident")

# Luisa genes
#Irf4, Klf4, Hdac1, 2, 3 and 11, Stat3. Stat6 and Trem2
sc.obj.MCells@active.assay <- 'RNA'
FeaturePlot(sc.obj.MCells, features = c("Irf4", "Klf4", "Hdac1", "Hdac2", "Hdac3",
                                        "Hdac11", "Stat3", "Stat6", "Trem2")
            , cols = c(alpha("grey", .5), "red"))
VlnPlot(sc.obj.MCells, features = c("Irf4", "Klf4", "Hdac1", "Hdac2", "Hdac3",
                                        "Hdac11", "Stat3", "Stat6", "Trem2"), 
        group.by = 'integrated_snn_res.0.3', pt.size = 0)

DimPlot(sc.obj.MCells, group.by = 'integrated_snn_res.0.3', label = T)


# Heatmap ----------------------------------------------------------------------

MCell_markers <- read.csv('Macro_DC_2/DiffExp_res03_filtered.csv')

MCell_markers %>% filter(avg_log2FC>.6) %>%
  dplyr::group_by(cluster)  %>% 
  filter(p_val_adj<1e-50) %>% 
  ungroup() %>% filter(cluster ==7)


goi_unique <- MCell_markers %>% filter(avg_log2FC>0.6) %>%
  dplyr::group_by(cluster)  %>% 
  filter(p_val_adj<1e-50) %>% 
  ungroup() %>%
  dplyr::count(gene) %>%
  # filter(n<=1) %>% 
  magrittr::use_series(gene)
# goi_unique <- goi_unique[!(stringr::str_detect(goi_unique, "^mt\\-|Rps|Rpl"))]

goi_new <- MCell_markers %>% 
  dplyr::group_by(cluster)  %>% 
  filter(avg_log2FC>0.6, gene %in% goi_unique, p_val_adj<1e-50) %>% 
  ungroup() %>%
  magrittr::use_series(gene) %>%
  unique()
labs <- c("M0", "M2_1?", "M1", "M2/MDSCs", "M2_2", "NC-Mono", "DC_1", "DC_2", "DC_3", "DC_4")
sc.obj.MCells$integrated_snn_res.0.3_names <-
  factor(sc.obj.MCells$integrated_snn_res.0.3, labels = labs)

avg_goi_M <- AverageExpression(sc.obj.MCells, features = goi_new, 
                               group.by = 'integrated_snn_res.0.3_names')


data.frame(real = 0:9, new = c(0,1,3,4,2,5,6,7,9, 8))

MCell_markers$order <- NA
vec <- 0:9
order_vec <- as.numeric( c(0,1,3,4,2,5,6,7,9, 8))
for (i in 1:10) {
  MCell_markers$order[MCell_markers$cluster==vec[i]] <- order_vec[i]
}


gene_order <- MCell_markers %>% 
  dplyr::group_by(cluster)  %>% 
  filter(avg_log2FC>0.6, gene %in% goi_unique, p_val_adj<1e-50) %>% 
  dplyr::arrange(p_val_adj) %>%
  ungroup() %>%
  select(gene,order) %>%
  filter(!duplicated(gene)==TRUE) %>%
  arrange(order) %>% magrittr::use_series(gene) %>% unique()

names(labs) <- 0:9
labs
avg_goi_M_df <- avg_goi_M$RNA[gene_order,c("M0","M2_1?","M2_2","M1","M2/MDSCs","NC-Mono",'DC_1',"DC_2","DC_4","DC_3")]

vec <- (alpha(c(as.vector(pals::tol(10))), .8))
names(vec) <- labs
annot_col <- list()
annot_col[['cluster']] = vec
ph1 <- pheatmap(avg_goi_M_df[gene_order,], 
                cluster_rows = F,
                # annotation_row = data.frame(row.names = goi_df$gene, cluster = as.factor(goi_df$cluster)),
                annotation_col = data.frame(row.names = labs, cluster = as.factor(labs)),
                color = rev(pals::brewer.rdbu(100)),
                annotation_colors = annot_col, main = "Clustering Heatmap Macrophages and DCs",
                show_rownames = T, scale = 'row', cluster_cols = F)


# annotate only selected genes
pl1 <- add.flag(ph1, c('Mrc1', 'Foxp3', 'Gzma', 'Ccl2', 'Cd19', 'Tex2',
                       'Siglech', 'Cd209d',"Havcr2", 'Cd4', 'Cd8a',
                       'Chil3', 'Camp', 'Lcn2', 'Ltf', 'Mmp9', 'Csf3r', 'Il1b', 'Ccl6'), 
                repel.degree = .1)

pl2 <- DimPlot(sc.obj.MCells, group.by = 'integrated_snn_res.0.3_names',
               cols = (alpha(c(as.vector(pals::tol(10))), .4)), label = TRUE, label.size = 10)
png('Macro_DC_2/UMAP_Heatmap_Temp_res03_MCells.png', width = 1500, height = 750)
cowplot::plot_grid(pl2, pl1,rel_widths = c(1,1))
dev.off()

MCell_markers$cluster_named <- factor(MCell_markers$cluster, labels = labs)



openxlsx::write.xlsx(MCell_markers %>% 
                       dplyr::group_by(cluster)  %>% 
                       filter(avg_log2FC>0.6, gene %in% goi_unique, p_val_adj<1e-50) %>%
                       select(-order), 
                     file = 'Macro_DC_2/Markers_4Heatmap_Mcells_211106.xlsx'
                     , asTable = T, overwrite = T, row.names = TRUE)



openxlsx::write.xlsx(data.frame(avg_goi_M_df) %>% mutate(gene = row.names(.)),
                     'Macro_DC_2/heatmap_df_MCells_211106.xlsx',
                     overwrite = T)
write.csv(avg_goi_M_df,
          'Macro_DC_2/heatmap_df_MCells_211106.csv')

