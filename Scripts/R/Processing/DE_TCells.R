## Diff Exp of cell types T Cells
##------ Sat Oct 30 16:39:00 2021 ------##

## Oliver Eichmueller

# load required libraries ------------------------------------------------------
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


sc.obj.TCells <- 
  readRDS(file = paste0(proc.Dir, "processed_RTC4_TCells_integration_211028.Rds"))

sc.obj.MCells <- 
  readRDS(file = paste0(proc.Dir, "processed_RTC4_MCells_integration_scores_211029.Rds"))



# regress out cell cycle -------------------------------------------------------
s.genes.mouse <- convertHumanGeneList(cc.genes$s.genes)
g2m.genes.mouse <- convertHumanGeneList(cc.genes$g2m.genes)
sc.obj.TCells <- CellCycleScoring(sc.obj.TCells, s.features = s.genes.mouse,
                                  g2m.features = g2m.genes.mouse)

DimPlot(sc.obj.TCells, group.by = "Phase")
sc.obj.TCells@active.assay <- 'integrated'

sc.obj.TCells_regr <- ScaleData(sc.obj.TCells, vars.to.regress = c("S.Score", "G2M.Score"),
                                features = row.names(sc.obj.TCells))

sc.obj.TCells_regr <- RunPCA(sc.obj.TCells_regr)

sc.obj.TCells_regr <- RunUMAP(sc.obj.TCells_regr, dims = 1:p$n.PC, 
                              n.components = 2)

pl1 <- DimPlot(sc.obj.TCells, group.by = c('integrated_snn_res.0.25', 'Phase'))
pl2 <- DimPlot(sc.obj.TCells_regr, group.by = c('integrated_snn_res.0.25', 'Phase'))
pl1+pl2
sc.obj.TCells_regr@active.assay <- "RNA"
FeaturePlot(sc.obj.TCells_regr, features = c("Klrb1c", "Cd8b1"), split.by = 'orig.ident')

saveRDS(sc.obj.TCells_regr, 'processed.DS_2/CC_Regress_processed_RTC4_TCells_integration_211031.Rds')

# DE of groups per identity ----------------------------------------------------

sc.obj.TCells$orig.ident <- factor(sc.obj.TCells$orig.ident)
sc.obj.TCells@active.ident <- factor(sc.obj.TCells$orig.ident)
clusters <- levels(sc.obj.TCells$integrated_snn_res.0.25)

DE_per_cluster <- list()

for (i in 1:length(clusters)) {
  sc.obj = subset(sc.obj.TCells, integrated_snn_res.0.25==clusters[i])
  
  DE_per_cluster[[paste0("cluster_", clusters[i])]] <-
    FindAllMarkers(sc.obj, assay = "RNA", only.pos = TRUE)
}

DE_per_cluster_filtered <- lapply(DE_per_cluster, function(x) x = x %>% filter(p_val_adj <1e-5) %>%
         dplyr::select(p_val_adj, cluster, gene) %>%
         pivot_wider(names_from = "cluster", values_from = "p_val_adj"))


DE_per_cluster_filtered_df <- plyr::ldply(DE_per_cluster_filtered, data.frame)

DE_per_cluster_filtered_df %>% filter(.id == "cluster_1")
DE_per_cluster_filtered_df %>% filter(gene %in% c("Bhlhe40"))

DimPlot(sc.obj.TCells, group.by = 'integrated_snn_res.0.25', label = TRUE)
sc.obj.TCells@active.assay <- "RNA"

imm_NK <- c("Ctla2a", "Ccr2", "Emb", "Cd28", "Thy1", "Sell", "Klra7")
mat_NK <- c("Klrg1", "Zeb2", "Cd69", "Klra9", "Klra8")
Eff_NK <- c("Gzma", "Icam1", "Ifng", "Irf1", "Prf1", "Il18r1", "Il18rap")
antiInfl_NK <- c("Dusp", "Zfp36", "Sgk1", "Klf2", "Klf6", "Fos", "Jun", "Ier2", "Ier5")

FeaturePlot(sc.obj.TCells, features = antiInfl_NK, split.by = 'orig.ident')
VlnPlot(sc.obj.TCells, features = antiInfl_NK, split.by = 'orig.ident',
        group.by = 'integrated_snn_res.0.25', idents = c(0,1))

sc.obj.TCells <- ScaleData(sc.obj.TCells)
goi <- intersect(row.names(sc.obj.TCells), c(imm_NK, mat_NK, Eff_NK, antiInfl_NK))
DoHeatmap(
  sc.obj.TCells, cells = names(sc.obj.TCells$integrated_snn_res.0.25)[sc.obj.TCells$integrated_snn_res.0.25=="1"]
          , features = goi, group.by = "orig.ident", group.bar = T, slot = 'counts') +
  scale_fill_continuous(type = 'viridis')

avg_NK_genes <- AverageExpression(sc.obj.TCells[,names(sc.obj.TCells$integrated_snn_res.0.25)[sc.obj.TCells$integrated_snn_res.0.25=="1"]]
                                  , assays = 'RNA', features = goi, group.by = 'orig.ident')
norM <- function(x){ x =x/max *10}

avg_NK_genes$RNA %>% as.data.frame() %>%
  mutate(max = rowMax(avg_NK_genes$RNA)) %>%
  mutate(across(1:5, ~ .x/max*10))
pheatmap(avg_NK_genes$RNA %>% as.data.frame() %>%
           mutate(max = rowMax(avg_NK_genes$RNA)) %>%
           mutate(across(1:5, ~ .x/max*10)) %>% select(-max), cluster_cols = F)


FeaturePlot(sc.obj.TCells, features = c("Tgfb1", "Smad7"), split.by = 'orig.ident')

VlnPlot(sc.obj.TCells, features = c("Ctsd"), group.by = 'orig.ident')
VlnPlot(sc.obj.TCells, features = c("Tgfb1", "Smad7", "Cotl1", "Car2", "Nfkb1"), split.by = 'orig.ident',
        group.by = 'integrated_snn_res.0.25', idents = c(1), ncol = 2, assay = 'RNA')
VlnPlot(sc.obj.TCells, features = c("Cxcr6", "Ccl4"),ncol = 1,
        group.by = 'integrated_snn_res.0.25')
row.names(sc.obj.TCells)[stringr::str_detect(row.names(sc.obj.TCells), "Rps")]

sc.obj.TCells[["percent.ribo"]] <- PercentageFeatureSet(sc.obj.TCells, pattern = "Rps")

FeaturePlot(sc.obj.MCells, features = c("Irf8"), split.by = 'orig.ident')

sc.obj.MCells@active.ident <- sc.obj.MCells$integrated_snn_res.0.3
VlnPlot(sc.obj.MCells, features = c("Tlr2", "Tlr4"), split.by = 'orig.ident',
        group.by = 'integrated_snn_res.0.3', idents = 3:5)


TCell_markers <- xlsx::read.xlsx('TCells_new/diff_exp/TCell_markers_seurat.xlsx', sheetIndex = 1)
DE_per_cluster_filtered_df <- xlsx::read.xlsx('TCells_new/diff_exp/DE_per_cluster_filtered_df_seurat.xlsx', sheetIndex = 1)

# Test DE with monocle ---------------------------------------------------------
library(monocle3)

cds.TCells <- readRDS('TCells_new/monocle/cds_TCells_211028_filtered.Rds')

plot_cells(cds.TCells, color_cells_by = 'integrated_snn_res.0.25')


model_fit <- fit_models(cds.TCells, model_formula_str = "~orig.ident * integrated_snn_res.0.25", 
                        cores = 4, verbose = TRUE)
fit_coefs <- coefficient_table(model_fit)

fit_coefs_clean <- fit_coefs %>% filter(q_value <0.05, estimate>0) %>%
  dplyr::select(gene_short_name, term, q_value) %>%
  pivot_wider(names_from = "term", values_from = "q_value")


fit_coefs_clean %>% 
  filter(gene_short_name=="Bhlhe40") %>% t()

fit_coefs %>% filter(q_value <0.05, estimate>0) %>%
  dplyr::select(gene_short_name, term, q_value) %>%
  filter(term %in% "orig.identrt-2:integrated_snn_res.0.258") %>%
  dplyr::arrange((q_value)) %>% head(20)

colnames(fit_coefs_clean)
id="integrated_snn_res.0.251"
fit_coefs_clean %>% filter(is.na(orig.identrtc4.1.integrated_snn_res.0.251)==FALSE) %>%
  dplyr::arrange(orig.identrtc4.1.integrated_snn_res.0.251)

TCell_markers <- read.csv(file = "TCells_new/all_markers_res025_TCells.csv")

TCell_markers %>% filter(p_val_adj < 1e-20, avg_log2FC>0, gene %in% c("Ccl3"))


fit_coefs_clean <- xlsx::read.xlsx('TCells_new/diff_exp/Model_fit_monocle.xlsx', sheetIndex = 1)
xlsx::write.xlsx(fit_coefs_clean, 'TCells_new/diff_exp/Model_fit_monocle.xlsx')
xlsx::write.xlsx(TCell_markers %>% filter(p_val_adj < 1e-20, avg_log2FC>0), 'TCells_new/diff_exp/TCell_markers_seurat.xlsx')
xlsx::write.xlsx(DE_per_cluster_filtered_df, 'TCells_new/diff_exp/DE_per_cluster_filtered_df_seurat.xlsx')

pl <- DimPlot(sc.obj.integration, group.by = "integrated_snn_res.0.2", combine = FALSE)

png('processed.plots_2/res02_transparent.png', bg="transparent")
pl[[1]] + theme(
  panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
  panel.grid.minor = element_blank(), 
  panel.grid.major = element_blank(),
  plot.background = element_rect(fill = "transparent",colour = NA))
dev.off()
ggsave(plot = pl, filename = 'TCells_new/res025_transparent.png', bg = "transparent", device = 'png')


# Test out Progeny ---------

CellsClusters <- data.frame(Cell = names(sc.obj.TCells$metacell_group), 
                            CellType = as.character(sc.obj.TCells$metacell_group),
                            stringsAsFactors = FALSE)
## We compute the Progeny activity scores and add them to our Seurat object
## as a new assay called Progeny. 
sc.obj.TCells <- progeny(sc.obj.TCells, scale=FALSE, organism="Mouse", top=500, perm=1, 
                  return_assay = TRUE)

## We can now directly apply Seurat functions in our Progeny scores. 
## For instance, we scale the pathway activity scores. 
sc.obj.TCells <- Seurat::ScaleData(sc.obj.TCells, assay = "progeny") 

## We transform Progeny scores into a data frame to better handling the results
progeny_scores_df <- 
  as.data.frame(t(GetAssayData(sc.obj.TCells, slot = "scale.data", 
                               assay = "progeny"))) %>%
  tibble::rownames_to_column("Cell") %>%
  gather(Pathway, Activity, -Cell) 

## We match Progeny scores with the cell clusters.
progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)

## We summarize the Progeny scores by cellpopulation
summarized_progeny_scores <- progeny_scores_df %>% 
  group_by(Pathway, CellType) %>%
  summarise(avg = mean(Activity), std = sd(Activity))

## We prepare the data for the plot
summarized_progeny_scores_df <- summarized_progeny_scores %>%
  dplyr::select(-std) %>%   
  spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) 


paletteLength = 100
myColor = colorRampPalette(c("Darkblue", "white","red"))(paletteLength)

progenyBreaks = c(seq(min(summarized_progeny_scores_df), 0, 
                      length.out=ceiling(paletteLength/2) + 1),
                  seq(max(summarized_progeny_scores_df)/paletteLength, 
                      max(summarized_progeny_scores_df), 
                      length.out=floor(paletteLength/2)))
progeny_hmap = pheatmap(t(summarized_progeny_scores_df[,-1]),fontsize=14, 
                        fontsize_row = 10, 
                        color=myColor, breaks = progenyBreaks, 
                        main = "PROGENy (500)", angle_col = 45,
                        treeheight_col = 0,  border_color = NA)

VlnPlot(sc.obj.TCells, features = "MAPK", split.by = 'orig.ident')


clusters <- levels(sc.obj.TCells$integrated_snn_res.0.25)

DE_per_cluster_progeny <- list()

for (i in 1:length(clusters)) {
  sc.obj = subset(sc.obj.TCells, integrated_snn_res.0.25==clusters[i])
  sc.obj@active.ident = factor(sc.obj$orig.ident)
  DE_per_cluster_progeny[[paste0("cluster_", clusters[i])]] <-
    FindAllMarkers(sc.obj, assay = "progeny", only.pos = TRUE)
}

DE_per_cluster_progeny_filtered <- 
  lapply(DE_per_cluster_progeny[c("cluster_0","cluster_1", "cluster_2",
                                  "cluster_3","cluster_4","cluster_5",
                                  "cluster_7","cluster_8","cluster_9")], 
         function(x) x = x %>% filter(p_val_adj <0.05) %>%
           dplyr::select(p_val_adj, cluster, gene) %>%
           pivot_wider(names_from = "cluster", values_from = "p_val_adj"))

DE_per_cluster_progeny_filtered <- plyr::ldply(DE_per_cluster_progeny_filtered, data.frame)


# DE Testing for tumor 1, 2 and untreated --------------------------------------

sc.obj.TCells$Radiation <- factor(sc.obj.TCells$orig.ident, labels = c("radiated", "unrad.T2", "radiated", "unrad.T2", "unrad.Ctrl"))

DimPlot(sc.obj.TCells, group.by = "Radiation")



clusters <- levels(sc.obj.TCells$integrated_snn_res.0.25)

DE_per_cluster_rad <- list()

for (i in 1:length(clusters)) {
  sc.obj = subset(sc.obj.TCells, integrated_snn_res.0.25==clusters[i])
  sc.obj@active.ident = factor(sc.obj$Radiation)
  DE_per_cluster_rad[[paste0("cluster_", clusters[i])]] <-
    FindAllMarkers(sc.obj, assay = "RNA", only.pos = TRUE)
}

DE_per_cluster_rad_filtered <- 
  lapply(DE_per_cluster_rad, 
         function(x) x = x %>% filter(p_val_adj <0.05) %>%
           dplyr::select(p_val_adj, cluster, gene) %>%
           pivot_wider(names_from = "cluster", values_from = "p_val_adj"))

DE_per_cluster_rad_filtered <- plyr::ldply(DE_per_cluster_rad_filtered, data.frame)
DE_per_cluster_rad_filtered %>%
  filter(.id == "cluster_1")


VlnPlot(sc.obj.TCells, features = "Gzmc", group.by = 'integrated_snn_res.0.25',
        split.by = "orig.ident", slot = 'data', idents = 1)
avg_NK <- AverageExpression(sc.obj.TCells[,names(sc.obj.TCells$integrated_snn_res.0.25)[sc.obj.TCells$integrated_snn_res.0.25=="1"]]
                            , assays = 'RNA', features = VariableFeatures(sc.obj.TCells, assay = 'RNA'),
                            group.by = 'orig.ident', slot = 'data')

avg_NK_df <- avg_NK$RNA

Enrichment_plot <- function(AvgExprDF, Title = "NK-cells",
                            cutoff = 5,
                            Xhigh = "rtc4-1", Xlow = "rt-2", 
                            Yhigh = "rt-1", Ylow = "rtc4-2",
                            UR_quadrant = "RT_induced",
                            BL_quadrant = "Tumor2",
                            BR_quadrant = "CTLA4_induced",
                            UL_quadrant = "FALSE"){
  
  plot_avg = data.frame(gene = row.names(AvgExprDF), 
                            X = AvgExprDF[,Xhigh] - AvgExprDF[,Xlow],
                            Y = AvgExprDF[,Yhigh] - AvgExprDF[,Ylow]) %>%
    mutate(label = 
             factor(ifelse((X>cutoff&Y>cutoff), UR_quadrant,
                           ifelse((X< -cutoff&Y< -cutoff), BL_quadrant, 
                                 ifelse((X> cutoff&Y< -cutoff),BR_quadrant, 
                                        ifelse((X< -cutoff&Y>cutoff), UL_quadrant, "FALSE")))),
                          levels = unique(c("FALSE", UR_quadrant, BL_quadrant, BR_quadrant, UL_quadrant))))
  
  label = plot_avg %>% filter(label != "FALSE") %>%
    magrittr::use_series(gene)
  
  pl = ggplot(plot_avg, aes(x=(X), y=(Y), color = label)) + geom_point() +
    geom_vline(xintercept = c(-cutoff, cutoff)) +
    geom_hline(yintercept = c(-cutoff, cutoff)) + 
    scale_color_manual("Variable Genes", values = c("grey", pals::brewer.set1((length(levels(plot_avg$label))-1))))+
    coord_equal() + 
    xlab(paste0(Xlow, " vs. ", Xhigh)) +
    ylab(paste0(Ylow, " vs. ", Yhigh)) +
    theme_classic() + ggtitle(Title)
  
  LabelPoints(pl, points = label, repel = T, xnudge = 0, ynudge = 0)
  
}
Enrichment_plot(AvgExprDF = avg_NK_df)


avg_cl4 <- AverageExpression(sc.obj.TCells[,names(sc.obj.TCells$integrated_snn_res.0.25)[sc.obj.TCells$integrated_snn_res.0.25=="4"]]
                            , assays = 'RNA', features = var_goi,
                            group.by = 'orig.ident', slot = 'data')

avg_cl4_df <- avg_cl4$RNA
Enrichment_plot(avg_cl4_df, "Tregs")

VlnPlot(sc.obj.TCells, features = "Tmsb4x", group.by = 'integrated_snn_res.0.25',
        split.by = "orig.ident", slot = 'data', idents = 4)

avg_cl8 <- AverageExpression(sc.obj.TCells[,names(sc.obj.TCells$integrated_snn_res.0.25)[sc.obj.TCells$integrated_snn_res.0.25=="8"]]
                            , assays = 'RNA', features = var_goi,
                            group.by = 'orig.ident', slot = 'data')

avg_cl0 <- AverageExpression(sc.obj.TCells[,names(sc.obj.TCells$integrated_snn_res.0.25)[sc.obj.TCells$integrated_snn_res.0.25=="0"]]
                            , assays = 'RNA', features = var_goi,
                            group.by = 'orig.ident', slot = 'data')
avg_cl3 <- AverageExpression(sc.obj.TCells[,names(sc.obj.TCells$integrated_snn_res.0.25)[sc.obj.TCells$integrated_snn_res.0.25=="3"]]
                            , assays = 'RNA', features = var_goi,
                            group.by = 'orig.ident', slot = 'data')
avg_cl7 <- AverageExpression(sc.obj.TCells[,names(sc.obj.TCells$integrated_snn_res.0.25)[sc.obj.TCells$integrated_snn_res.0.25=="7"]]
                            , assays = 'RNA', features = var_goi,
                            group.by = 'orig.ident', slot = 'data')

var_goi <- VariableFeatures(sc.obj.TCells, assay = 'RNA')[VariableFeatures(sc.obj.TCells, assay = 'RNA')!="Gm42418"]


avg_cl8_df <- avg_cl8$RNA
avg_cl3_df <- avg_cl3$RNA
avg_cl7_df <- avg_cl7$RNA
avg_cl3_df['Bhlhe40',]
avg_cl7_df['Bhlhe40',]
Enrichment_plot(avg_cl7_df, "imm Cd8", cutoff = 3)
VlnPlot(sc.obj.TCells, features = "Junb", group.by = 'integrated_snn_res.0.25',
        split.by = "orig.ident", slot = 'data', idents = 7)


# new DE just comparing Rad with no rad ----------------------------------------
clusters <- levels(sc.obj.TCells$integrated_snn_res.0.25)

DE_per_cluster_radvsunrad <- list()

for (i in 1:length(clusters)) {
  sc.obj = subset(sc.obj.TCells, integrated_snn_res.0.25==clusters[i])
  sc.obj@active.ident = factor(sc.obj$Radiation)
  DE_per_cluster_radvsunrad[[paste0("cluster_", clusters[i])]] <-
    FindMarkers(sc.obj,ident.1 = "radiated", ident.2 = "unrad.T2", assay = "RNA")
}

DE_per_cluster_radvsunrad_filtered <- 
  lapply(DE_per_cluster_radvsunrad, 
         function(x) x = x %>%
           dplyr::select(p_val_adj, avg_log2FC))

openxlsx::write.xlsx(lapply(DE_per_cluster_radvsunrad_filtered,
                            function(x) x = x%>%
                              dplyr::mutate(gene =row.names(.))), file = 'TCells_new/diff_exp/VolcanoPlots/Radiation_vs_Unradiated/DE_per_cluster_radvsunrad_filtered.xlsx'
                     , asTable = T, overwrite = T, row.names = TRUE)


# Compare CTLA4 ----------------------------------------------------------------
sc.obj.TCells$CTLA4 <- factor(sc.obj.TCells$orig.ident, labels = c("no", "no", "yes", "yes", "Ctrl"))

clusters <- levels(sc.obj.TCells$integrated_snn_res.0.25)

DE_per_cluster_CTLA4vsno <- list()

for (i in 1:length(clusters)) {
  sc.obj = subset(sc.obj.TCells, integrated_snn_res.0.25==clusters[i])
  sc.obj@active.ident = sc.obj$CTLA4
  DE_per_cluster_CTLA4vsno[[paste0("cluster_", clusters[i])]] <-
    FindMarkers(sc.obj,ident.1 = "yes", ident.2 = "no", assay = "RNA")
}

DE_per_cluster_CTLA4vsno_filtered <- 
  lapply(DE_per_cluster_CTLA4vsno, 
         function(x) x = x %>%
           dplyr::select(p_val_adj, avg_log2FC))

openxlsx::write.xlsx(lapply(DE_per_cluster_CTLA4vsno_filtered,
                            function(x) x = x%>%
                              dplyr::mutate(gene =row.names(.))), file = 'TCells_new/diff_exp/VolcanoPlots/CTLA4_vs_no/DE_per_cluster_CTLA4vsno_filtered.xlsx'
                     , asTable = T, overwrite = T, row.names = TRUE)




clusters <- names(DE_per_cluster_CTLA4vsno_filtered)
labs <- c("Dividing", "NK", "Cd8-EM", "naive Cd4, Cd8", "Treg", "Cd8-Eff", "Cd8 mito-low", 
          "naive Cd8", "Th1", "Dividing Cd8 mito-low")

EnhancedVolcano(DE_per_cluster_CTLA4vsno_filtered[[clusters[1]]], 
                lab = rownames(DE_per_cluster_CTLA4vsno_filtered[[clusters[1]]]),
                x = 'avg_log2FC',
                y = 'p_val_adj')

i=10
png(paste0('TCells_new/diff_exp/VolcanoPlots/CTLA4_vs_no/', labs[i], "_",
           clusters[i], "_volcano.png"), width = 800, height = 800)
EnhancedVolcano(DE_per_cluster_CTLA4vsno_filtered[[clusters[i]]], 
                lab = rownames(DE_per_cluster_CTLA4vsno_filtered[[clusters[i]]]),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75, title = paste0(clusters[i], " ", labs[i]), 
                subtitle = "CTLA4 vs. no CTLA4 (irresp. of radiation)")

dev.off()


VlnPlot(sc.obj.TCells, features = c("Cd52"), group.by = 'integrated_snn_res.0.25',
        split.by = "orig.ident", slot = 'data')

FeaturePlot(sc.obj.TCells, features = "Bhlhe40", split.by = 'orig.ident')
FeaturePlot(sc.obj.TCells, features = "Bhlhe40")
VlnPlot(sc.obj.TCells, features = "Bhlhe40", group.by = 'orig.ident',
         slot = 'data', idents = c(7,8))
VlnPlot(sc.obj.TCells, features = "Bhlhe40",group.by = 'integrated_snn_res.0.25',
        split.by = "orig.ident",
         slot = 'data', idents = c(7,8))

# for whole dataset
all_t_ctla4_de <- FindMarkers(sc.obj.TCells, ident.1 = "yes", ident.2 = "no", group.by = "CTLA4")
png(paste0('TCells_new/diff_exp/VolcanoPlots/CTLA4_vs_no/AllTCells_volcano_Bhlhe40.png'), 
    width = 800, height = 800)
EnhancedVolcano(all_t_ctla4_de, 
                lab = rownames(all_t_ctla4_de),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                legendIconSize = 4.0,selectLab = 'Bhlhe40',
                FCcutoff = .5,
                drawConnectors = TRUE,
                widthConnectors = 0.75, title = "CTLA4 vs. no CTLA4", 
                subtitle = "irresp. of radiation")
dev.off()



avg_metacell_group <- AverageExpression(sc.obj.TCells
                  , assays = 'RNA',
                  group.by = 'metacell_group', slot = 'data')

openxlsx::write.xlsx(as.data.frame(all_t_ctla4_de)%>% mutate(gene=row.names(.)), 
                     file = 'TCells_new/diff_exp/VolcanoPlots/CTLA4_vs_no/DE_allTcells_CTLA4vsno_filtered.xlsx'
                     , overwrite = TRUE)
# compare just within tumor 1 --------------------------------------------------

clusters <- levels(sc.obj.TCells$integrated_snn_res.0.25)

DE_per_cluster_t1 <- list()

for (i in 1:length(clusters)) {
  sc.obj = subset(sc.obj.TCells, integrated_snn_res.0.25==clusters[i])
  sc.obj@active.ident = factor(sc.obj$orig.ident)
  DE_per_cluster_t1[[paste0("cluster_", clusters[i])]] <-
    FindMarkers(sc.obj,ident.1 = "rt-1", ident.2 = "rtc4-1", assay = "RNA")
}



clusters <- names(DE_per_cluster_CTLA4vsno_filtered)
labs <- c("Dividing", "NK", "Cd8-EM", "naive Cd4, Cd8", "Treg", "Cd8-Eff", "Cd8 mito-low", 
          "naive Cd8", "Th1", "Dividing Cd8 mito-low")

EnhancedVolcano(DE_per_cluster_t1[[clusters[1]]], 
                lab = rownames(DE_per_cluster_t1[[clusters[1]]]),
                x = 'avg_log2FC',
                y = 'p_val_adj')

i=10
png(paste0('TCells_new/diff_exp/VolcanoPlots/Tumor1_rtvscombi/', labs[i], "_",
           clusters[i], "_volcano.png"), width = 800, height = 800)
EnhancedVolcano(DE_per_cluster_t1[[clusters[i]]], 
                lab = rownames(DE_per_cluster_t1[[clusters[i]]]),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75, title = paste0(clusters[i], " ", labs[i]), 
                subtitle = "RT vs. Combination in Tumor 1")

dev.off()

openxlsx::write.xlsx(lapply(DE_per_cluster_t1,
                            function(x) x = x%>%
                              dplyr::mutate(gene =row.names(.))), 
                     file = 'TCells_new/diff_exp/VolcanoPlots/Tumor1_rtvscombi/DE_per_cluster_t1.xlsx'
                     , asTable = T, overwrite = T, row.names = TRUE)

# compare just within tumor 2 --------------------------------------------------

clusters <- levels(sc.obj.TCells$integrated_snn_res.0.25)

DE_per_cluster_t2 <- list()

for (i in 1:length(clusters)) {
  sc.obj = subset(sc.obj.TCells, integrated_snn_res.0.25==clusters[i])
  sc.obj@active.ident = factor(sc.obj$orig.ident)
  DE_per_cluster_t2[[paste0("cluster_", clusters[i])]] <-
    FindMarkers(sc.obj,ident.1 = "rt-2", ident.2 = "rtc4-2", assay = "RNA")
}



clusters <- names(DE_per_cluster_t2)
labs <- c("Dividing", "NK", "Cd8-EM", "naive Cd4, Cd8", "Treg", "Cd8-Eff", "Cd8 mito-low", 
          "naive Cd8", "Th1", "Dividing Cd8 mito-low")

EnhancedVolcano(DE_per_cluster_t2[[clusters[1]]], 
                lab = rownames(DE_per_cluster_t2[[clusters[1]]]),
                x = 'avg_log2FC',
                y = 'p_val_adj')

i=10
png(paste0('TCells_new/diff_exp/VolcanoPlots/Tumor2_rtvscombi/', labs[i], "_",
           clusters[i], "_volcano.png"), width = 800, height = 800)
EnhancedVolcano(DE_per_cluster_t2[[clusters[i]]], 
                lab = rownames(DE_per_cluster_t2[[clusters[i]]]),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75, title = paste0(clusters[i], " ", labs[i]), 
                subtitle = "RT vs. Combination in Tumor 2")

dev.off()

openxlsx::write.xlsx(lapply(DE_per_cluster_t2,
                            function(x) x = x%>%
                              dplyr::mutate(gene =row.names(.))), 
                     file = 'TCells_new/diff_exp/VolcanoPlots/Tumor2_rtvscombi/DE_per_cluster_t2.xlsx'
                     , asTable = T, overwrite = T, row.names = TRUE)

# RT vs ctrl -------------------------------------------------------------------

clusters <- levels(sc.obj.TCells$integrated_snn_res.0.25)

DE_per_cluster_t1vsctrl <- list()

for (i in 1:length(clusters)) {
  sc.obj = subset(sc.obj.TCells, integrated_snn_res.0.25==clusters[i])
  sc.obj@active.ident = sc.obj$Radiation
  DE_per_cluster_t1vsctrl[[paste0("cluster_", clusters[i])]] <-
    FindMarkers(sc.obj,ident.1 = "radiated", ident.2 = "unrad.Ctrl", assay = "RNA")
}



clusters <- names(DE_per_cluster_t1vsctrl)
labs <- c("Dividing", "NK", "Cd8-EM", "naive Cd4, Cd8", "Treg", "Cd8-Eff", "Cd8 mito-low", 
          "naive Cd8", "Th1", "Dividing Cd8 mito-low")

i=10
png(paste0('TCells_new/diff_exp/VolcanoPlots/Tumor1_vs_Ctrl/', labs[i], "_",
           clusters[i], "_volcano.png"), width = 800, height = 800)
EnhancedVolcano(DE_per_cluster_t1vsctrl[[clusters[i]]], 
                lab = rownames(DE_per_cluster_t1vsctrl[[clusters[i]]]),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75, title = paste0(clusters[i], " ", labs[i]), 
                subtitle = "RT vs. Ctrl")

dev.off()

openxlsx::write.xlsx(lapply(DE_per_cluster_t1vsctrl,
                            function(x) x = x%>%
                              dplyr::mutate(gene =row.names(.))), 
                     file = 'TCells_new/diff_exp/VolcanoPlots/Tumor1_vs_Ctrl/DE_per_cluster_t1vsctrl.xlsx'
                     , asTable = T, overwrite = T, row.names = TRUE)

DE_All_t1vsctrl <- 
  FindMarkers(sc.obj.TCells,ident.1 = "radiated", ident.2 = "unrad.Ctrl", assay = "RNA"
              , group.by = "Radiation")

png(paste0('TCells_new/diff_exp/VolcanoPlots/Tumor1_vs_Ctrl/', "All_MacroDC_volcano.png"), width = 800, height = 800)
EnhancedVolcano(DE_All_t1vsctrl, 
                lab = rownames(DE_All_t1vsctrl),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75, title = "All T-cells and NK", 
                subtitle = "RT vs. Ctrl in Tumor 1")

dev.off()

openxlsx::write.xlsx(DE_All_t1vsctrl%>%
                       dplyr::mutate(gene =row.names(.)), 
                     file = 'TCells_new/diff_exp/VolcanoPlots/Tumor1_vs_Ctrl/DE_All_t1vsctrl.xlsx'
                     , asTable = T, overwrite = T, row.names = TRUE)

# T2 vs ctrl -------------------------------------------------------------------

clusters <- levels(sc.obj.TCells$integrated_snn_res.0.25)

DE_per_cluster_t2vsctrl <- list()

for (i in 1:length(clusters)) {
  sc.obj = subset(sc.obj.TCells, integrated_snn_res.0.25==clusters[i])
  sc.obj@active.ident = sc.obj$Radiation
  DE_per_cluster_t2vsctrl[[paste0("cluster_", clusters[i])]] <-
    FindMarkers(sc.obj,ident.1 = "unrad.T2", ident.2 = "unrad.Ctrl", assay = "RNA")
}



clusters <- names(DE_per_cluster_t2vsctrl)
labs <- c("Dividing", "NK", "Cd8-EM", "naive Cd4, Cd8", "Treg", "Cd8-Eff", "Cd8 mito-low", 
          "naive Cd8", "Th1", "Dividing Cd8 mito-low")

i=10
png(paste0('TCells_new/diff_exp/VolcanoPlots/Tumor2_vs_Ctrl/', labs[i], "_",
           clusters[i], "_volcano.png"), width = 800, height = 800)
EnhancedVolcano(DE_per_cluster_t2vsctrl[[clusters[i]]], 
                lab = rownames(DE_per_cluster_t2vsctrl[[clusters[i]]]),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75, title = paste0(clusters[i], " ", labs[i]), 
                subtitle = "T2 vs. Ctrl")

dev.off()

openxlsx::write.xlsx(lapply(DE_per_cluster_t2vsctrl,
                            function(x) x = x%>%
                              dplyr::mutate(gene =row.names(.))), 
                     file = 'TCells_new/diff_exp/VolcanoPlots/Tumor2_vs_Ctrl/DE_per_cluster_t2vsctrl.xlsx'
                     , asTable = T, overwrite = T, row.names = TRUE)

DE_All_t2vsctrl <- 
  FindMarkers(sc.obj.TCells,ident.1 = "unrad.T2", ident.2 = "unrad.Ctrl", assay = "RNA"
              , group.by = "Radiation")

png(paste0('TCells_new/diff_exp/VolcanoPlots/Tumor2_vs_Ctrl/', "All_MacroDC_volcano.png"), width = 800, height = 800)
EnhancedVolcano(DE_All_t2vsctrl, 
                lab = rownames(DE_All_t2vsctrl),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75, title = "All T-cells and NK", 
                subtitle = "T2 vs. Ctrl")

dev.off()

openxlsx::write.xlsx(DE_All_t2vsctrl%>%
                       dplyr::mutate(gene =row.names(.)), 
                     file = 'TCells_new/diff_exp/VolcanoPlots/Tumor2_vs_Ctrl/DE_All_t2vsctrl.xlsx'
                     , asTable = T, overwrite = T, row.names = TRUE)


# Combi vs ctrl -------------------------------------------------------------------

clusters <- levels(sc.obj.TCells$integrated_snn_res.0.25)

DE_per_cluster_combivsctrl <- list()

for (i in 1:length(clusters)) {
  sc.obj = subset(sc.obj.TCells, integrated_snn_res.0.25==clusters[i])
  sc.obj@active.ident = factor(sc.obj$orig.ident)
  DE_per_cluster_combivsctrl[[paste0("cluster_", clusters[i])]] <-
    FindMarkers(sc.obj,ident.1 = "rtc4-2", ident.2 = "untreated-1", assay = "RNA")
}



clusters <- names(DE_per_cluster_combivsctrl)
labs <- c("Dividing", "NK", "Cd8-EM", "naive Cd4, Cd8", "Treg", "Cd8-Eff", "Cd8 mito-low", 
          "naive Cd8", "Th1", "Dividing Cd8 mito-low")
 
i=10
png(paste0('TCells_new/diff_exp/VolcanoPlots/CombiT2_vs_Ctrl/', labs[i], "_",
           clusters[i], "_volcano.png"), width = 800, height = 800)
EnhancedVolcano(DE_per_cluster_combivsctrl[[clusters[i]]], 
                lab = rownames(DE_per_cluster_combivsctrl[[clusters[i]]]),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75, title = paste0(clusters[i], " ", labs[i]), 
                subtitle = "Combination Tumor 2 vs. Ctrl")

dev.off()


openxlsx::write.xlsx(lapply(DE_per_cluster_combivsctrl,
                            function(x) x = x%>%
                              dplyr::mutate(gene =row.names(.))), 
                     file = 'TCells_new/diff_exp/VolcanoPlots/CombiT2_vs_Ctrl/DE_per_cluster_combivsctrl.xlsx'
                     , asTable = T, overwrite = T, row.names = TRUE)

DE_All_t2vsctrl <- 
  FindMarkers(sc.obj.TCells,ident.1 = "unrad.T2", ident.2 = "unrad.Ctrl", assay = "RNA"
              , group.by = "Radiation")



png(paste0('TCells_new/diff_exp/VolcanoPlots/Tumor2_vs_Ctrl/', "All_MacroDC_volcano.png"), width = 800, height = 800)
EnhancedVolcano(DE_All_t2vsctrl, 
                lab = rownames(DE_All_t2vsctrl),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75, title = "All T-cells and NK", 
                subtitle = "T2 vs. Ctrl")

dev.off()

openxlsx::write.xlsx(DE_All_t2vsctrl%>%
                       dplyr::mutate(gene =row.names(.)), 
                     file = 'TCells_new/diff_exp/VolcanoPlots/Tumor2_vs_Ctrl/DE_All_t2vsctrl.xlsx'
                     , asTable = T, overwrite = T, row.names = TRUE)

# Plot Heatmap on T-cells ------------------------------------------------------



goi_unique <- TCell_markers %>% filter(avg_log2FC>0.6) %>%
  dplyr::group_by(cluster)  %>% 
  filter(p_val_adj<1e-50) %>% 
  ungroup() %>%
  dplyr::count(gene) %>%
  # filter(n<=1) %>% 
  magrittr::use_series(gene)
# goi_unique <- goi_unique[!(stringr::str_detect(goi_unique, "^mt\\-|Rps|Rpl"))]

goi_new <- TCell_markers %>% 
  dplyr::group_by(cluster)  %>% 
  filter(avg_log2FC>0.6, gene %in% goi_unique, 
         p_val_adj<1e-50) %>% 
  ungroup() %>%
  magrittr::use_series(gene) %>%
  unique()

labs <- c("Dividing Eff", "NK", "Cd8-EM", "naive Cd4, Cd8", "Treg", "Cd8 Act/Exh", "Cd8 mt-high", "naive Cd8", "Th1", "Dividing Cd8 mt-high")

sc.obj.TCells$integrated_snn_res.0.25_named <-
  factor(sc.obj.TCells$integrated_snn_res.0.25, labels = labs)

DimPlot(sc.obj.TCells, group.by = 'integrated_snn_res.0.25_named', label = T)

avg_goi_T <- AverageExpression(sc.obj.TCells, features = goi_new, 
                                 group.by = 'integrated_snn_res.0.25_named')

gene_order <- TCell_markers %>% 
   dplyr::group_by(cluster)  %>% 
   filter(avg_log2FC>0.6, gene %in% goi_unique, p_val_adj<1e-50) %>% 
  dplyr::arrange(p_val_adj) %>%
   ungroup() %>%
   select(gene,order) %>%
     filter(!duplicated(gene)==TRUE) %>%
   arrange(order) %>% magrittr::use_series(gene) %>% unique()


data.frame(real = 0:9, new = c(0,9,3,5,7,4,2,6,8, 1))

TCell_markers$order <- NA
vec <- 0:9
order_vec <- as.numeric( c(0,9,3,5,7,4,2,6,8, 1))
for (i in 1:10) {
  TCell_markers$order[TCell_markers$cluster==vec[i]] <- order_vec[i]
}

data.frame(0:9, cl = levels(sc.obj.TCells$integrated_snn_res.0.25_named))

avg_goi_T_df <- avg_goi_T$RNA[gene_order,c("Dividing Eff","Dividing Cd8 mt-high","Cd8 mt-high","Cd8-EM","Cd8 Act/Exh",
                                           "naive Cd4, Cd8",'naive Cd8',"Treg","Th1","NK")]

vec <- (alpha(c(as.vector(pals::tol(10))), .8))
names(vec) <- levels(sc.obj.TCells$integrated_snn_res.0.25_named)
annot_col <- list()
annot_col[['cluster']] = vec
ph1 <- pheatmap(avg_goi_T_df, 
                cluster_rows = F,
                annotation_col = data.frame(row.names = levels(sc.obj.TCells$integrated_snn_res.0.25_named), 
                                            cluster = as.factor(levels(sc.obj.TCells$integrated_snn_res.0.25_named))),
                color = rev(pals::brewer.rdbu(100)),
                annotation_colors = annot_col, main = "Clustering Heatmap T- and NK-cells",
                show_rownames = T, scale = 'row', cluster_cols = F)


# annotate only selected genes
pl1 <- add.flag(ph1, c('Foxp3', 'Gzma', 'Mki67', 'Top2a', 'mt-Atp6', "Klrb1c", "S1pr1",
                       'Klrd1', 'Pdcd1', 'Lag3',"Havcr2", 'Cd4', 'Cd8a', "Rpl12"
                       ), repel.degree = .1)

pl2 <- DimPlot(sc.obj.TCells, group.by = 'integrated_snn_res.0.25_named',
               cols = (alpha(c(as.vector(pals::tol(10))), .4)), label = TRUE, label.size = 5)
png('TCells_new/UMAP_Heatmap_Temp_res03_TCells.png', width = 1500, height = 750)
cowplot::plot_grid(pl2, pl1,rel_widths = c(1,1))
dev.off()


TCell_markers$cluster_named <- factor(TCell_markers$cluster,
                                      labels = levels(sc.obj.TCells$integrated_snn_res.0.25_named))


openxlsx::write.xlsx(TCell_markers %>% 
                       dplyr::group_by(cluster)  %>% 
                       filter(avg_log2FC>0.6, gene %in% goi_unique, 
                              p_val_adj<1e-50) %>% select(-order), 
                     file = 'TCells_new/Markers_4Heatmap_Tcells_211106.xlsx'
                     , asTable = T, overwrite = T, row.names = TRUE)

openxlsx::write.xlsx(data.frame(avg_goi_T_df) %>% mutate(gene = row.names(.)),
                     'TCells_new/heatmap_df_TCells_211106.xlsx',
                     overwrite = T)
write.csv(avg_goi_T_df,
          'TCells_new/heatmap_df_TCells_211106.csv')



# Directly compare mito high vs other effector cells ---------------------------

mito_vs_EM <- FindMarkers(sc.obj.TCells, group.by = 'integrated_snn_res.0.25_named',
                          ident.1 = 'Cd8 mt-high', ident.2 = 'Cd8-EM', assay = 'RNA')

goi <- mito_vs_EM %>% filter(p_val_adj <1e-20, avg_log2FC>0.5|avg_log2FC< -0.5) %>% 
  mutate(gene = row.names(.)) %>%
  magrittr::use_series(gene)



goi_others <- TCell_markers %>% filter(p_val_adj <1e-50, avg_log2FC>0.5, cluster_named != "Cd8 mt-high",
                         cluster_named !="Cd8-EM")  %>% magrittr::use_series(gene)
setdiff(goi,goi_others)

avg_T_mito <- AverageExpression(sc.obj.TCells, features = setdiff(goi,goi_others), 
                               group.by = 'integrated_snn_res.0.25_named')

genes_mito_ordered <- row.names(mito_vs_EM[setdiff(goi,goi_others),] %>% arrange(desc(avg_log2FC)))

avg_T_mito_df <- avg_T_mito$RNA[genes_mito_ordered,c("Dividing Eff","Dividing Cd8 mt-high","Cd8 mt-high","Cd8-EM","Cd8 Act/Exh",
                                           "naive Cd4, Cd8",'naive Cd8',"Treg","Th1","NK")]

vec <- (alpha(c(as.vector(pals::tol(10))), .8))
names(vec) <- levels(sc.obj.TCells$integrated_snn_res.0.25_named)
annot_col <- list()
annot_col[['cluster']] = vec
ph1 <- pheatmap(avg_T_mito_df, 
                cluster_rows = F,
                annotation_col = data.frame(row.names = levels(sc.obj.TCells$integrated_snn_res.0.25_named), 
                                            cluster = as.factor(levels(sc.obj.TCells$integrated_snn_res.0.25_named))),
                color = rev(pals::brewer.rdbu(100)),
                annotation_colors = annot_col, main = "Mt-high vs mt-low Cd8",
                show_rownames = T, scale = 'row', cluster_cols = F)
png('TCells_new/Heatmap_mthighvslow_21106.png')
pl1 <- add.flag(ph1, c('Il18r1', 'Il18rap', 'Macf1', 'Itgal', 'Kmt2a', "Huwe1", "Hmha1",
                       'Usp34', 'Son', 'Gmfg',"Btf3", 'Ubb', 'S100a10', "Cd52", "Fxyd5", "Tomm7"), repel.degree = .1)
dev.off()


openxlsx::write.xlsx(mito_vs_EM[genes_mito_ordered,] %>% 
                       mutate(gene = row.names(.)) %>%
                      arrange(desc(avg_log2FC)), 
                     file = 'TCells_new/Markers_mitohighvslow_211106.xlsx'
                     , asTable = T, overwrite = T, row.names = TRUE)

# Directly compare mito high lineage vs mit-low lineage ------------------------

cells_use <- row.names(sc.obj.TCells@meta.data %>% filter(orig.ident %in% c("rt-2", "rtc4-2")))

mito_vs_EM_lineage <- FindMarkers(sc.obj.TCells[,cells_use], group.by = 'integrated_snn_res.0.25_named',
                          ident.1 = c('Cd8 mt-high', 'Dividing Cd8 mt-high'), 
                          ident.2 = c('Cd8-EM', 'Dividing Eff'), assay = 'RNA')

goi2 <- mito_vs_EM_lineage %>% filter(p_val_adj <1e-20, avg_log2FC>0.5|avg_log2FC< -0.5) %>% 
  mutate(gene = row.names(.)) %>%
  magrittr::use_series(gene)

goi_others <- TCell_markers %>% filter(p_val_adj <1e-50, avg_log2FC>0.5, cluster_named != "Cd8 mt-high",
                         cluster_named !="Cd8-EM",
                         cluster_named !="Dividing Cd8 mt-high",
                         cluster_named !="Dividing Eff"
                         )  %>% magrittr::use_series(gene)
setdiff(setdiff(goi2,goi_others), setdiff(goi,goi_others))

avg_T_mito_lineage <- AverageExpression(sc.obj.TCells[,cells_use], features = setdiff(goi2,goi_others), 
                               group.by = 'integrated_snn_res.0.25_named')

genes_mito_ordered <- row.names(mito_vs_EM_lineage[setdiff(goi2,goi_others),] %>% arrange(desc(avg_log2FC)))

avg_T_mito_df_lineage <- avg_T_mito_lineage$RNA[genes_mito_ordered,c("Dividing Eff","Dividing Cd8 mt-high","Cd8 mt-high","Cd8-EM","Cd8 Act/Exh",
                                           "naive Cd4, Cd8",'naive Cd8',"Treg","Th1","NK")]

vec <- (alpha(c(as.vector(pals::tol(10))), .8))
names(vec) <- levels(sc.obj.TCells$integrated_snn_res.0.25_named)
annot_col <- list()
annot_col[['cluster']] = vec
ph1 <- pheatmap(avg_T_mito_df, 
                cluster_rows = F,
                annotation_col = data.frame(row.names = levels(sc.obj.TCells$integrated_snn_res.0.25_named), 
                                            cluster = as.factor(levels(sc.obj.TCells$integrated_snn_res.0.25_named))),
                color = rev(pals::brewer.rdbu(100)),
                annotation_colors = annot_col, main = "Mt-high vs mt-low Cd8",
                show_rownames = T, scale = 'row', cluster_cols = F)

png('TCells_new/Heatmap_mthighvslow_lineage_21106.png')
pl1 <- add.flag(ph1, c('Il18r1', 'Il18rap', 'Macf1', 'Itgal', 'Kmt2a', "Huwe1", "Hmha1", "Xist", "Macf1", "Phip", "Ikzf1",
                       'Usp34', 'Son', 'Gmfg',"Btf3", 'Ubb', 'S100a10', "Cd52", "Fxyd5", "Tomm7"), repel.degree = .1)
dev.off()


openxlsx::write.xlsx(mito_vs_EM[genes_mito_ordered,] %>% 
                       mutate(gene = row.names(.)) %>%
                      arrange(desc(avg_log2FC)), 
                     file = 'TCells_new/Markers_mitohighvslow_lineage_211106.xlsx'
                     , asTable = T, overwrite = T, row.names = TRUE)


