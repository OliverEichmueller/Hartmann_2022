library(pheatmap)
library(ggplot2)
library(CodeAndRoll2)
library(monocle3)
library(ComplexHeatmap)
library(EnhancedVolcano)
library(Seurat)
library(UCell)
library(dplyr)
library(ggplot2)

color_paper <- list(
  `rt-1` = "#ffc080",
  `rt-2` = "#c06000",
  `rtc4-1` = "#90bff9",
  `rtc4-2` = "#0000c0",
  `untreated-1` = "#a0a0a4"
)



score_df_M_all <- openxlsx::read.xlsx(xlsxFile = 'Macro_DC_2/DiffExp/GO_heatmaps/Heatmap_SourceData_211118.xlsx',
                                  rowNames = T)
score_df_M_sel <- openxlsx::read.xlsx(xlsxFile = 'Macro_DC_2/DiffExp/GO_heatmaps/Heatmap_SourceData_211118.xlsx',
                                  rowNames = T, sheet = 2)
score_df_M_oth <- openxlsx::read.xlsx(xlsxFile = 'Macro_DC_2/DiffExp/GO_heatmaps/Heatmap_SourceData_211118.xlsx',
                                  rowNames = T, sheet = 3)
score_df_gubin <- openxlsx::read.xlsx(xlsxFile = 'Macro_DC_2/DiffExp/GO_heatmaps/Heatmap_SourceData_211118.xlsx',
                                      sheet = 4, rowNames = T)
score_df_gubin <- score_df_gubin[,c(1:4,9,5:8)]


annot <- data.frame(row.names = colnames(score_df_gubin),
                    sample = c("aCTLA4", "aCTLA4-PD1", "aPD1", "ctrl",
                               "untreated-1", "C12-iF", "C12-oF",
                               "C12/aC4-iF", "C12/aC4-oF"),
                    source = c(rep("Gubin et al.", 4), rep("Hartmann et al.", 5)))
colnames(annot) <- c("sample", "source")


sample_col <- as.vector(c(pals::tol(4), as.named.vector(color_paper)[c(5,1:4)]))
names(sample_col) <- annot$sample
source_col <- as.vector(pals::tol(2))
names(source_col) <- annot$source %>% unique()

annot_col <- list(sample = sample_col, source = source_col)

pdf('Macro_DC_2/DiffExp/GO_heatmaps/GO_heatmap_comparison_Gubin_rearrange.pdf', width = 10, height = 4)
print(pheatmap::pheatmap(score_df_gubin, color = pals::coolwarm(100), cluster_cols = F,
                           gaps_col = 4, 
                           annotation_col = annot, cluster_rows = F,
                           show_colnames = F, annotation_colors = annot_col))
dev.off()

# Selection of clusters
score_df_M_sel <- score_df_M_sel[,c(5,1:4,10,6:9,15,11:14,20,16:19,25,21:24,30,26:29,35,31:34,40,36:39,45,41:44)]

annot <- data.frame(row.names = colnames(score_df_M_sel),
                    sample = rep(c("untreated-1", "C12-iF", "C12-oF",
                                 "C12/aC4-iF", "C12/aC4-oF"),9),
                    celltype = stringr::str_extract(colnames(score_df_M_sel), "(?<=_).{1,}"))
colnames(annot) <- c("sample", "cluster")


celltype_col <-c(as.vector(pals::tol(12)), "grey")
names(celltype_col) <- stringr::str_extract(colnames(score_df_M_all), "(?<=_).{1,}") %>% unique()
sample_col <- as.named.vector(color_paper)[c(5,1:4)]
names(sample_col) <- c("untreated-1", "C12-iF", "C12-oF",
                       "C12/aC4-iF", "C12/aC4-oF")
annot_col <- list(cluster = celltype_col, sample = sample_col)

pdf('Macro_DC_2/DiffExp/GO_heatmaps/GO_heatmap_selection_rearrange.pdf', width = 15, height = 4)
print(pheatmap(as.matrix(score_df_M_sel), color = pals::coolwarm(100), cluster_cols = F,
                           gaps_col = seq(0,45, by = 5), 
                           annotation_col = annot,
                           show_colnames = F, annotation_colors = annot_col))
dev.off()


# All samples ------
score_df_M_all <-score_df_M_all[,c(5,1:4,10,6:9,15,11:14,20,16:19,25,21:24,30,26:29,35,31:34,40,36:39,45,41:44,50,46:49,55,51:54,60,56:59,65,61:64)]


annot <- data.frame(row.names = colnames(score_df_M_all),
                    sample = rep(c("untreated-1", "C12-iF", "C12-oF",
                                   "C12/aC4-iF", "C12/aC4-oF"),13),
                    celltype = stringr::str_extract(colnames(score_df_M_all), "(?<=_).{1,}"))
colnames(annot) <- c("sample", "cluster")

annot$cluster <- factor(as.vector(annot$cluster), levels = unique(annot$cluster))

celltype_col <-c(as.vector(pals::tol(12)), "grey")
names(celltype_col) <- stringr::str_extract(colnames(score_df_M_all), "(?<=_).{1,}") %>% unique()
sample_col <- as.named.vector(color_paper)[c(5,1:4)]
names(sample_col) <- c("untreated-1", "C12-iF", "C12-oF",
                       "C12/aC4-iF", "C12/aC4-oF")
annot_col <- list(cluster = celltype_col, sample = sample_col)


pdf('Macro_DC_2/DiffExp/GO_heatmaps/GO_heatmap_all_rearrange.pdf', width = 15, height = 4)
print(pheatmap(as.matrix(score_df_M_all), color = pals::coolwarm(100), 
               cluster_cols = F,
               gaps_col = seq(0,65, by = 5), 
               annotation_col = annot,
               show_colnames = F, annotation_colors = annot_col))
dev.off()




# Do on whole integrated dataset for comparison

score_df_M_oth <- score_df_M_oth[,c(5,1:4,10,6:9,15,11:14,20,16:19,25,21:24,30,26:29,35,31:34)]

annot <- data.frame(row.names = colnames(score_df_M_oth),
                    sample = rep(c("untreated-1", "C12-iF", "C12-oF",
                                   "C12/aC4-iF", "C12/aC4-oF"),7),
                    celltype = stringr::str_replace(stringr::str_extract(colnames(score_df_M_oth), "(?<=_).{1,}"), "\\.", " "))
colnames(annot) <- c("sample", "cluster")


celltype_col <-c(as.vector(pals::tol(12)), "grey")[c(3,6:9,11:12)]
names(celltype_col) <- stringr::str_replace(stringr::str_extract(colnames(score_df_M_oth), "(?<=_).{1,}"), "\\.", " ")%>%
  unique()
sample_col <- as.named.vector(color_paper)[c(5,1:4)]
names(sample_col) <- c("untreated-1", "C12-iF", "C12-oF",
                       "C12/aC4-iF", "C12/aC4-oF")
annot_col <- list(cluster = celltype_col, sample = sample_col)

pdf('Macro_DC_2/DiffExp/GO_heatmaps/GO_heatmap_OtherCellTypes_rearranged.pdf', width = 15, height = 4)
print(pheatmap(score_df_M_oth, color = pals::coolwarm(100), cluster_cols = F,
                           gaps_col = seq(0,30, by = 5), 
                           annotation_col = annot,
                           show_colnames = F, annotation_colors = annot_col))
dev.off()


openxlsx::write.xlsx(list(`All Clusters` = score_df_M_all, `Selection Clusters` = score_df_M_sel,
                          `Other Cell Types` = score_df_M_oth, 
                          `Comparison Gubin Macrophages` = score_df_gubin),
                     row.names = TRUE, overwrite = T,
                     file = 'Macro_DC_2/DiffExp/GO_heatmaps/Heatmap_SourceData_rearrange_211209.xlsx')

