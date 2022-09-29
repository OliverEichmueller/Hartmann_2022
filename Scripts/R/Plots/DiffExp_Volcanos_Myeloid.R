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


# Compare CTLA4 within res05 ------------------------------------------------

clusters <- levels(sc.obj.MCells$integrated_snn_res.0.5_named)

MCells_DE_per_cluster_CTLA4vsno <- list()

for (i in 1:length(clusters)) {
  sc.obj = subset(sc.obj.MCells, integrated_snn_res.0.5_named==clusters[i])
  sc.obj@active.ident = sc.obj$CTLA4
  MCells_DE_per_cluster_CTLA4vsno[[paste0("cluster_", clusters[i])]] <-
    FindMarkers(sc.obj,ident.1 = "yes", ident.2 = "no", assay = "RNA")
}

MCells_DE_per_cluster_CTLA4vsno_filtered <- 
  lapply(MCells_DE_per_cluster_CTLA4vsno, 
         function(x) x = x %>%
           dplyr::select(p_val_adj, avg_log2FC))


openxlsx::write.xlsx(lapply(MCells_DE_per_cluster_CTLA4vsno_filtered,
                            function(x) x = x%>%
                              dplyr::mutate(gene =row.names(.))), 
                     file = 'Macro_DC_2/VolcanoPlots/CTLA4vsno/MCells_DE_per_cluster_CTLA4vsno_filtered_res05.xlsx'
                     , asTable = T, overwrite = T, row.names = TRUE)

clusters <- names(MCells_DE_per_cluster_CTLA4vsno_filtered)


labs <- levels(sc.obj.MCells$integrated_snn_res.0.5_named)

for (i in 1:length(labs)) {
  
  png(paste0('Macro_DC_2/VolcanoPlots/CTLA4vsno_res05/',
             clusters[i], "_volcano.png"), width = 800, height = 800)
  print(EnhancedVolcano(MCells_DE_per_cluster_CTLA4vsno_filtered[[clusters[i]]], 
                  lab = rownames(MCells_DE_per_cluster_CTLA4vsno_filtered[[clusters[i]]]),
                  x = 'avg_log2FC',
                  y = 'p_val_adj',
                  legendIconSize = 4.0,
                  drawConnectors = TRUE,
                  widthConnectors = 0.75, title = paste0(clusters[i], " ", labs[i]), 
                  subtitle = "CTLA4 vs. no CTLA4 (irresp. of RT)"))
  
  dev.off()
}

DimPlot(sc.obj.MCells, group.by = "integrated_snn_res.0.5_grouped", label = T)


# Compare grouped per cluster res05 ------------------------------------------------
M_cellsDE_ALL_CTLA4vsno_grouped <- 
  FindMarkers(sc.obj.MCells,ident.1 = "yes", ident.2 = "no", assay = "RNA"
              , group.by = "CTLA4")


clusters <- levels(sc.obj.MCells$integrated_snn_res.0.5_grouped)

MCells_DE_per_cluster_CTLA4vsno_grouped <- list()

for (i in 1:length(clusters)) {
  sc.obj = subset(sc.obj.MCells, integrated_snn_res.0.5_grouped==clusters[i])
  sc.obj@active.ident = sc.obj$CTLA4
  MCells_DE_per_cluster_CTLA4vsno_grouped[[paste0("cluster_", clusters[i])]] <-
    FindMarkers(sc.obj,ident.1 = "yes", ident.2 = "no", assay = "RNA")
}

MCells_DE_per_cluster_CTLA4vsno_grouped_filtered <- 
  lapply(MCells_DE_per_cluster_CTLA4vsno_grouped, 
         function(x) x = x %>%
           dplyr::select(p_val_adj, avg_log2FC))


labs <- levels(sc.obj.MCells$integrated_snn_res.0.5_grouped)
clusters <- names(MCells_DE_per_cluster_CTLA4vsno_grouped_filtered)

for (i in 1:length(labs)) {
  
  png(paste0('Macro_DC_2/VolcanoPlots/CTLA4vsno_res05_grouped/',
             clusters[i], "_volcano.png"), width = 800, height = 800)
  print(EnhancedVolcano(MCells_DE_per_cluster_CTLA4vsno_grouped_filtered[[clusters[i]]], 
                        lab = rownames(MCells_DE_per_cluster_CTLA4vsno_grouped_filtered[[clusters[i]]]),
                        x = 'avg_log2FC',
                        y = 'p_val_adj',
                        legendIconSize = 4.0,
                        FCcutoff = 0.5,
                        drawConnectors = TRUE,
                        widthConnectors = 0.75, title = paste0(clusters[i]), 
                        subtitle = "CTLA4 vs. no CTLA4 (irresp. of RT)"))
  
  dev.off()
}
openxlsx::write.xlsx(lapply(MCells_DE_per_cluster_CTLA4vsno_grouped_filtered,
                            function(x) x = x%>%
                              dplyr::mutate(gene =row.names(.))), 
                     file = 'Macro_DC_2/VolcanoPlots/CTLA4vsno_res05_grouped/MCells_DE_grouped_ctla4.xlsx'
                     , asTable = T, overwrite = T, row.names = TRUE)


# Compare TAM groups CTLA4 ------------------------------------------------
MCells_DE_per_cluster_CTLA4vsno_grouped <- list()

sc.obj = subset(sc.obj.MCells, integrated_snn_res.0.5_grouped%in% c("TAM-Group-1", "TAM-Group-2", "TAM-Group-3"))
sc.obj@active.ident = sc.obj$CTLA4
MCells_DE_per_cluster_CTLA4vsno_TAM <-
  FindMarkers(sc.obj,ident.1 = "yes", ident.2 = "no", assay = "RNA")

keep <- stringr::str_subset(row.names(MCells_DE_per_cluster_CTLA4vsno_TAM), "Gm", negate = T) 
MCells_DE_per_cluster_CTLA4vsno_TAM[keep,]

pdf('Macro_DC_2/VolcanoPlots/CTLA4vsno_AllTAMs/Volcano_allTAMS.pdf')
print(EnhancedVolcano(MCells_DE_per_cluster_CTLA4vsno_TAM[keep,], 
                      lab = rownames(MCells_DE_per_cluster_CTLA4vsno_TAM[keep,]),
                      x = 'avg_log2FC',
                      y = 'p_val_adj',
                      legendIconSize = 4.0,
                      FCcutoff = 0.5,
                      drawConnectors = TRUE,
                      widthConnectors = 0.75, title = paste0("All TAMs"), 
                      subtitle = "CTLA4 vs. no CTLA4 (irresp. of RT)"))
dev.off()

openxlsx::write.xlsx(MCells_DE_per_cluster_CTLA4vsno_TAM[keep,]%>%
                              dplyr::mutate(gene =row.names(.)), 
                     file = 'Macro_DC_2/VolcanoPlots/CTLA4vsno_AllTAMs/MCells_DE_AllTAMs_ctla4.xlsx'
                     , asTable = T, overwrite = T, row.names = TRUE)

# Compare T1 vs ctrl with res05 ------------------------------------------------
clusters <- levels(sc.obj.MCells$integrated_snn_res.0.5_named)

MCells_DE_per_cluster_t1vsctrl <- list()

for (i in 1:length(clusters)) {
  sc.obj = subset(sc.obj.MCells, integrated_snn_res.0.5_named==clusters[i])
  sc.obj@active.ident = sc.obj$Radiation
  MCells_DE_per_cluster_t1vsctrl[[paste0("cluster_", clusters[i])]] <-
    FindMarkers(sc.obj,ident.1 = "radiated.T1", ident.2 = "unrad.Ctrl", assay = "RNA")
}



clusters <- names(MCells_DE_per_cluster_t1vsctrl)
labs <- levels(sc.obj.MCells$integrated_snn_res.0.5_named)

for (i in 1:length(labs)) {
  
  png(paste0('Macro_DC_2/VolcanoPlots/Tumor1_vs_Ctrl_res05/',
             clusters[i], "_volcano.png"), width = 800, height = 800)
  print(EnhancedVolcano(MCells_DE_per_cluster_t1vsctrl[[clusters[i]]], 
                  lab = rownames(MCells_DE_per_cluster_t1vsctrl[[clusters[i]]]),
                  x = 'avg_log2FC',
                  y = 'p_val_adj',
                  legendIconSize = 4.0,
                  FCcutoff = 0.5,
                  drawConnectors = TRUE,
                  widthConnectors = 0.75, title = paste0(clusters[i]), 
                  subtitle = "C12-iF and C12/aCF-iF vs. Ctrl"))
  
  dev.off()
}

openxlsx::write.xlsx(lapply(MCells_DE_per_cluster_t1vsctrl,
                            function(x) x = x%>%
                              dplyr::mutate(gene =row.names(.))), 
                     file = 'Macro_DC_2/VolcanoPlots/Tumor1_vs_Ctrl_res05/MCells_DE_per_cluster_t1vsctrl_res05.xlsx'
                     , asTable = T, overwrite = T, row.names = TRUE)

MCells_DE_All_t1vsctrl <- 
  FindMarkers(sc.obj.MCells,ident.1 = "radiated.T1", ident.2 = "unrad.Ctrl", assay = "RNA"
              , group.by = "Radiation")

png(paste0('Macro_DC_2/VolcanoPlots/Tumor1_vs_Ctrl_res05/', "All_MacroDC_volcano.png"), width = 800, height = 800)
EnhancedVolcano(MCells_DE_All_t1vsctrl, 
                lab = rownames(MCells_DE_All_t1vsctrl),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75, title = "All Macro and DC", 
                subtitle = "C12-iF and C12/aCF-iF vs. Ctrl")

dev.off()

openxlsx::write.xlsx(MCells_DE_All_t1vsctrl%>%
                       dplyr::mutate(gene =row.names(.)), 
                     file = 'Macro_DC_2/VolcanoPlots/Tumor1_vs_Ctrl_res05/MCells_DE_All_t1vsctrl.xlsx'
                     , asTable = T, overwrite = T, row.names = TRUE)




