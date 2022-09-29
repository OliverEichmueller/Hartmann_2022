## Heatmaps for figure
## Mon Dec  6 11:33:22 2021 ----------------------------------------------------
## Oliver Eichmueller
library(gridExtra)
library(openxlsx)
library(grid)
library(zoo)
library(ggplot2)
library(Seurat)
library(tidyr)
library(SeuratDisk)
library(dplyr)
require(MarkdownReports)
require(pheatmap)


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


sc.obj.all <- readRDS('processed.DS_2/processed_RTC4_integration_211115_scores.Rds')

## T Cells ---------------------------------------------------------------------
sc.obj.TCells <- readRDS('processed.DS_2/processed_RTC4_TCells_integration_scores_211121.Rds')


TCell_markers <- openxlsx::read.xlsx('/Users/Oliver.Eichmueller/Downloads/Markers_4Heatmap_Tcells_211106.xlsx', rowNames = T)


goi_unique <- TCell_markers %>% filter(avg_log2FC>0.6) %>%
  dplyr::group_by(cluster)  %>% 
  filter(p_val_adj<1e-50) %>% 
  ungroup() %>%
  dplyr::count(gene) %>%
  # filter(n<=1) %>% 
  magrittr::use_series(gene)

goi_new <- TCell_markers %>% 
  dplyr::group_by(cluster)  %>% 
  filter(avg_log2FC>0.6, gene %in% goi_unique, 
         p_val_adj<1e-50) %>% 
  ungroup() %>%
  magrittr::use_series(gene) %>%
  unique()

labs <- c("Dividing Eff", "NK", "Cd8-EM", "naive Cd4, Cd8", "Treg", "Cd8 Act/Exh", 
          "Cd8 mt-high", "naive Cd8", "Th1", "Dividing Cd8 mt-high")

DimPlot(sc.obj.TCells, group.by = 'integrated_snn_res.0.25_Final', label = T)

avg_goi_T <- AverageExpression(sc.obj.TCells, features = goi_new, 
                               group.by = 'integrated_snn_res.0.25_Final')

data.frame(0:9, cl = levels(sc.obj.TCells$integrated_snn_res.0.25_Final))
data.frame(real = 0:9, new = c(0,9,6,2,5,3,7,4,8, 1), cl = levels(sc.obj.TCells$integrated_snn_res.0.25_Final),
           gene =as.numeric( c(0,9,3,5,7,4,2,6,8, 1)))

TCell_markers$order <- NA
vec <- 0:9
order_vec <- as.numeric( c(0,9,3,5,7,4,2,6,8, 1))
for (i in 1:10) {
  TCell_markers$order[TCell_markers$cluster==vec[i]] <- order_vec[i]
}

gene_order <- TCell_markers %>% 
  dplyr::group_by(cluster)  %>% 
  filter(avg_log2FC>0.6, gene %in% goi_unique, p_val_adj<1e-50) %>% 
  dplyr::arrange(p_val_adj, desc(avg_log2FC)) %>%
  ungroup() %>%
  select(gene,order) %>%
  filter(!duplicated(gene)==TRUE) %>%
  arrange(order) %>% magrittr::use_series(gene) %>% unique()


ordered <- c(0,9,6,2,5,3,7,4,8, 1)+1

avg_goi_T_df <- avg_goi_T$RNA[gene_order,levels(sc.obj.TCells$integrated_snn_res.0.25_Final)[ordered]]

vec <- (alpha(c(as.vector(pals::tol(10))), .8))
names(vec) <- levels(sc.obj.TCells$integrated_snn_res.0.25_Final)
vec <- vec[ordered]
annot_col <- list()
annot_col[['cluster']] = vec

ph1 <- pheatmap(avg_goi_T_df[,1:10], 
                cluster_rows = F,
                annotation_col = data.frame(row.names = levels(sc.obj.TCells$integrated_snn_res.0.25_Final), 
                                            cluster = as.factor(levels(sc.obj.TCells$integrated_snn_res.0.25_Final))),
                color = rev(pals::brewer.rdbu(100)),
                annotation_colors = annot_col, main = "Clustering Heatmap T- and NK-cells",
                show_rownames = T, show_colnames = F, scale = 'row', cluster_cols = F)


# annotate only selected genes
annot_genes <- TCell_markers %>% filter(is.na(CellType) == F) %>% magrittr::use_series(gene)

annot_genes2 <- c("Mki67", "Top2a", "Aurkb", "Cxcr6", "Il18r1", "Il18rap","Lag3"
                  ,"Havcr2", "Ifng", "Cd160", "Pdcd1", "Bhlhe40", "Litaf"
                  , "Lef1", "Il6ra", "Il7r","S1pr1", "Ccr7", "Foxp3"
                  ,"Itgb8", "Ctla4", "Il2ra"
                  ,"Tnfsf11", "Tnfsf8", "Gzma", "Klrb1c", "Prf1")
  



pl1 <- add.flag(ph1, annot_genes2, repel.degree = .1)
dev.off()

goi2 <- TCell_markers %>% 
  dplyr::group_by(cluster)  %>% 
  filter(avg_log2FC>0.6, gene %in% goi_unique, p_val_adj<1e-50) %>% 
  filter(gene %in% c(annot_genes2)) %>%
  dplyr::arrange(p_val_adj, desc(avg_log2FC)) %>%
  ungroup() %>%
  select(gene,order) %>%
  filter(!duplicated(gene)==TRUE) %>%
  arrange(order) %>% magrittr::use_series(gene) %>% unique()

color_new_row <- TCell_markers %>% 
  dplyr::group_by(cluster)  %>% 
  filter(avg_log2FC>0.6, gene %in% goi_unique, p_val_adj<1e-50) %>% 
  dplyr::arrange(p_val_adj, desc(avg_log2FC)) %>%
  ungroup() %>%
  select(gene,order) %>%
  filter(!duplicated(gene)==TRUE) %>%
  arrange(order) %>% left_join(data.frame(names = levels(sc.obj.TCells$integrated_snn_res.0.25_Final)[ordered], order = 0:9),
                               by = "order") %>%
  left_join(data.frame(vec, names = names(vec)), by = "names") %>%
  magrittr::use_series(vec)

color_new_lines <- TCell_markers %>% 
  dplyr::group_by(cluster)  %>% 
  filter(avg_log2FC>0.6, gene %in% goi_unique, p_val_adj<1e-50) %>% 
  filter(gene %in% c(annot_genes2)) %>%
  dplyr::arrange(p_val_adj, desc(avg_log2FC)) %>%
  ungroup() %>%
  select(gene,order) %>%
  filter(!duplicated(gene)==TRUE) %>%
  arrange(order) %>% left_join(data.frame(names = levels(sc.obj.TCells$integrated_snn_res.0.25_Final)[ordered], order = 0:9),
                               by = "order") %>%
  left_join(data.frame(vec, names = names(vec)), by = "names") %>%
  magrittr::use_series(vec)


pl1$grobs[[1]]$gp=gpar(fontsize=20, fontface = "bold")

pl1$grobs[[3]]$gp=gpar(col=color_new_row, fontsize=20, fontface = "bold")
pl1$grobs[[6]]$children$GRID.text.2339$gp=gpar(fontsize=20, fontface = "bold", vjust = -1)
pl1$grobs[[6]]$children$GRID.text.2339$vjust=0.3
pl1$grobs[[6]]$children$GRID.rect.2340$width=unit(1, "cm")
pl1$grobs[[6]]$children$GRID.rect.2340$height=unit(1, "cm")
pl1$grobs[[6]]$children$GRID.rect.2340$y=rev(unit((1:10)+5, "cm"))
pl1$grobs[[6]]$children$GRID.text.2341$y=rev(unit((1:10)+4.5, "cm"))
pl1$grobs[[6]]$children$GRID.text.2341$x=unit(1.1, "cm")
pl1$grobs[[6]]$children$GRID.text.2341$gp=gpar(fontsize=15, fontface = "bold", vjust = -1)
pl1$grobs[[8]]$x0=unit(0, "cm")
pl1$grobs[[8]]$gp=gpar(col=color_new_lines, lwd=2)
dev.off()

png('Heatmaps_paper_final/TCells_Heatmap_colored.png', width = 800, height = 800)
grid.arrange(pl1)
dev.off()
pdf('Heatmaps_paper_final/TCells_Heatmap_colored.pdf', width = 10, height = 10)
grid.arrange(pl1)
dev.off()
group <- TCell_markers %>% 
  dplyr::group_by(cluster)  %>% 
  filter(avg_log2FC>0.6, gene %in% goi_unique, p_val_adj<1e-50) %>% 
  dplyr::arrange(p_val_adj) %>%
  ungroup() %>%
  select(gene,order) %>%
  filter(!duplicated(gene)==TRUE) %>%
  arrange(order) %>% left_join(data.frame(names = levels(sc.obj.TCells$integrated_snn_res.0.25_Final)[ordered], order = 0:9),
                               by = "order") %>% magrittr::use_series(names)
avg_goi_T_df <- as.data.frame(avg_goi_T_df)
avg_goi_T_df$maxPval <- group

TCell_markers_2 <- TCell_markers

TCell_markers_2$cluster_named2 <- factor(TCell_markers_2$cluster,
                                         labels = levels(sc.obj.TCells$integrated_snn_res.0.25_Final))

write.xlsx(list(`T Cells Heatmap unscaled` = avg_goi_T_df, 
                `T Cells All DEGs` = TCell_markers_2[,c(1:8, 14)]), overwrite = T,
                     file = 'Heatmaps_paper_final/supp_tables/TCells_Supp_heatmap.xlsx',
                    rowNames = T)

# All Cells --------------------------------------------------------------------

All_markers <- openxlsx::read.xlsx('/Users/Oliver.Eichmueller/Downloads/Markers_4Heatmap_AllCells_211106 (1).xlsx')

sc.obj.all$integrated_snn_res.0.2_named <-
  factor(sc.obj.all$integrated_snn_res.0.2_named, 
         labels = c("TAM-A", "TAM-B", levels(sc.obj.all$integrated_snn_res.0.2_named)[3:13]))

goi_unique <- All_markers %>% filter(avg_log2FC>0.6) %>%
  dplyr::group_by(cluster)  %>% 
  filter(p_val_adj<1e-50) %>% 
  ungroup() %>%
  dplyr::count(gene) %>%
  # filter(n<=1) %>% 
  magrittr::use_series(gene)

goi_unique <- goi_unique[!(stringr::str_detect(goi_unique, "^mt\\-|Rps|Rpl"))]

All_markers$order <- NA
vec <- 0:12
order_vec <- data.frame(orig = 0:12, labels = colnames(avg_goi_all$RNA)) %>%
  left_join(data.frame(new = 0:12, labels = colnames(avg_goi_all_df2)), by = "labels") %>%
  magrittr::use_series(new)
for (i in 1:13) {
  All_markers$order[All_markers$cluster==vec[i]] <- order_vec[i]
}

goi_new <- All_markers %>% 
  dplyr::group_by(cluster)  %>% 
  filter(avg_log2FC>0.6, gene %in% goi_unique, p_val_adj<1e-100) %>% 
  dplyr::arrange(p_val_adj, desc(avg_log2FC)) %>%
  ungroup() %>%
  select(gene,order) %>%
  filter(!duplicated(gene)==TRUE) %>%
  arrange(order) %>% magrittr::use_series(gene) %>% unique()

avg_goi_all <- AverageExpression(sc.obj.all, features = goi_new, assays = "RNA",
                                 group.by = 'integrated_snn_res.0.2_named')


order_all <- c(0,1,3,11,4,9,2,8,10,5,6,7,12)+1

avg_goi_all_df2 <- avg_goi_all$RNA[goi_new,order_all]

vec <- (alpha(c(as.vector(pals::tol(12)), "dark grey"), .8))
names(vec) <- levels(sc.obj.all$integrated_snn_res.0.2_named)
vec <- vec[order_all]
annot_col <- list()
annot_col[['cluster']] = vec
ph1 <- pheatmap(avg_goi_all_df2, 
                cluster_rows = F,
                # annotation_row = data.frame(row.names = goi_df$gene, cluster = as.factor(goi_df$cluster)),
                annotation_col = data.frame(row.names = levels(sc.obj.all$integrated_snn_res.0.2_named), 
                                            cluster = as.factor(levels(sc.obj.all$integrated_snn_res.0.2_named))),
                color = rev(pals::brewer.rdbu(100)),
                annotation_colors = annot_col, main = "Clustering Heatmap all cells",
                show_rownames = T, scale = 'row', cluster_cols = F)


# annotate only selected genes
annot_genes_all <-  All_markers %>% 
  dplyr::group_by(cluster)  %>% 
  filter(avg_log2FC>0.6, gene %in% goi_unique, p_val_adj<1e-100) %>% 
  dplyr::arrange(p_val_adj, desc(avg_log2FC)) %>%
  filter(is.na(Celltype) == FALSE) %>% 
  ungroup() %>%
  select(gene, order) %>%
  filter(!duplicated(gene)==TRUE) %>%
  ungroup() %>% 
  arrange(order) %>%magrittr::use_series(gene)

annot_genes_all_new <- c("Cxcl10","Csf1r","Itgam","Arg1",
"Apoe","Mrc1",
"Havcr2","Lcn2","Mmp9",
"Ccl22","Il12b","Xcr1","Cd8b1",
"Ifng","Cd8a", "Il18r1",
"Il2ra","Foxp3","S1pr1", "Tcf7","Gzma","Klrb1c",
"Pax5","Cd22","Cd19","Siglech",
"Cd209d","Tex2")


ph1 <-pheatmap(avg_goi_all_df2[,], 
         cluster_rows = F,
         # annotation_row = data.frame(row.names = goi_df$gene, cluster = as.factor(goi_df$cluster)),
         annotation_col = data.frame(row.names = levels(sc.obj.all$integrated_snn_res.0.2_named), 
                                     cluster = as.factor(levels(sc.obj.all$integrated_snn_res.0.2_named))),
         color = rev(pals::brewer.rdbu(100)),
         annotation_colors = annot_col, main = "Clustering Heatmap all cells",
         show_rownames = T, show_colnames = F, scale = 'row', cluster_cols = F)





pl1 <- add.flag(ph1, annot_genes_all_new, repel.degree = .1)

color_new_row <- All_markers %>% 
  dplyr::group_by(cluster)  %>% 
  filter(avg_log2FC>0.6, gene %in% goi_unique, p_val_adj<1e-50) %>% 
  dplyr::arrange(p_val_adj, desc(avg_log2FC)) %>%
  ungroup() %>%
  select(gene,order) %>%
  filter(!duplicated(gene)==TRUE) %>%
  arrange(order) %>% left_join(data.frame(names = levels(sc.obj.all$integrated_snn_res.0.2_named)[order_all], 
                                          order = 0:12),
                               by = "order") %>%
  left_join(data.frame(vec, names = names(vec)), by = "names") %>%
  magrittr::use_series(vec)

color_new_lines <-  All_markers %>% 
  dplyr::group_by(cluster)  %>% 
  filter(avg_log2FC>0.6, gene %in% goi_unique, p_val_adj<1e-100) %>% 
  filter(gene %in% annot_genes_all_new) %>%
  dplyr::arrange(p_val_adj, desc(avg_log2FC)) %>%
  ungroup() %>%
  select(gene,order) %>%
  filter(!duplicated(gene)==TRUE) %>%
  arrange(order) %>% left_join(data.frame(names = levels(sc.obj.all$integrated_snn_res.0.2_named)[order_all], 
                                          order = 0:12),
                               by = "order") %>%
  left_join(data.frame(vec, names = names(vec)), by = "names") %>%
  magrittr::use_series(vec)



pl1$grobs[[1]]$gp=gpar(fontsize=20, fontface = "bold")

pl1$grobs[[3]]$gp=gpar(col=color_new_row, fontsize=20, fontface = "bold")
pl1$grobs[[6]]$children$GRID.text.2997$gp=gpar(fontsize=20, fontface = "bold", vjust = -1)
pl1$grobs[[6]]$children$GRID.text.2997$vjust=0.3
pl1$grobs[[6]]$children$GRID.rect.2998$width=unit(1, "cm")
pl1$grobs[[6]]$children$GRID.rect.2998$height=unit(1, "cm")
pl1$grobs[[6]]$children$GRID.rect.2998$y=rev(unit((1:13)+5, "cm"))
pl1$grobs[[6]]$children$GRID.text.2999$y=rev(unit((1:13)+4.5, "cm"))
pl1$grobs[[6]]$children$GRID.text.2999$x=unit(1.1, "cm")
pl1$grobs[[6]]$children$GRID.text.2999$gp=gpar(fontsize=15, fontface = "bold", vjust = -1)
pl1$grobs[[8]]$x0=unit(0, "cm")
pl1$grobs[[8]]$gp=gpar(col=color_new_lines, lwd=2)
dev.off()
grid.arrange(pl1)



png('Heatmaps_paper_final/AllCells_Heatmap_colored.png', width = 800, height = 800)
grid.arrange(pl1)
dev.off()
pdf('Heatmaps_paper_final/AllCells_Heatmap_colored.pdf', width = 10, height = 10)
grid.arrange(pl1)
dev.off()




All_markers$cluster_named <-
  factor(All_markers$cluster, labels = levels(sc.obj.all$integrated_snn_res.0.2_named))



avg_goi_all_df3 <- as.data.frame(avg_goi_all_df2)
avg_goi_all_df3$clusters <- All_markers %>% 
  dplyr::group_by(cluster)  %>% 
  filter(avg_log2FC>0.6, gene %in% goi_unique, p_val_adj<1e-100) %>% 
  dplyr::arrange(p_val_adj, desc(avg_log2FC)) %>%
  ungroup() %>%
  filter(!duplicated(gene)==TRUE) %>%
  arrange(order) %>% magrittr::use_series(cluster_named)



write.xlsx(list(`All Cells Heatmap unscaled` = avg_goi_all_df3, 
                `All Cells All DEGs` = All_markers[,c(1:9)]), overwrite = T,
           file = 'Heatmaps_paper_final/supp_tables/AllCells_Supp_heatmap.xlsx',
           rowNames = T)

# M Cells --------------------------------------------------------------------

M_markers <- openxlsx::read.xlsx('/Users/Oliver.Eichmueller/Downloads/Markers_4Heatmap_Mcells_res05_211117.xlsx')

sc.obj.MCells <- readRDS('processed.DS_2/processed_RTC4_MCells_211115_scores.Rds')

goi_unique <- M_markers %>% filter(avg_log2FC>0.6) %>%
  dplyr::group_by(cluster)  %>% 
  filter(p_val_adj<1e-50) %>% 
  ungroup() %>%
  dplyr::count(gene) %>%
  # filter(n<=1) %>% 
  magrittr::use_series(gene)

goi_new <- M_markers %>% 
  dplyr::group_by(cluster)  %>% 
  filter(avg_log2FC>0.6, gene %in% goi_unique, 
         p_val_adj<1e-50) %>% 
  ungroup() %>%
  magrittr::use_series(gene) %>%
  unique()


DimPlot(sc.obj.MCells, group.by = 'integrated_snn_res.0.5_named', label = T)

avg_goi_M <- AverageExpression(sc.obj.MCells, features = goi_new, assays = "RNA",
                               group.by = 'integrated_snn_res.0.5_named')

data.frame(0:12, cl = levels(sc.obj.MCells$integrated_snn_res.0.5_named))


reorder_df <- data.frame(real = 0:12, new = c(9,5,2,3,4,6,1,0,8,11,7,10,12), cl = levels(sc.obj.MCells$integrated_snn_res.0.5_named)) %>%
  left_join(data.frame(cl=c("Monocytes", "NC-Monocytes", "TAM-3", "TAM-4", "TAM-5", "TAM-6", "TAM-2", "TAM-1", 
                                   "TAM-7", "TRM", "cDC2", "cDC1-1", "cDC1-2"), gene = 0:12),
            by = "cl")

M_markers$order <- NA
vec <- reorder_df$cl
order_vec <- reorder_df$gene
for (i in 1:13) {
  M_markers$order[M_markers$cluster==vec[i]] <- order_vec[i]
}

gene_order <- M_markers %>% 
  dplyr::group_by(cluster)  %>% 
  filter(avg_log2FC>0.6, gene %in% goi_unique, p_val_adj<1e-50) %>% 
  dplyr::arrange(p_val_adj, desc(avg_log2FC)) %>%
  ungroup() %>%
  select(gene,order) %>%
  filter(!duplicated(gene)==TRUE) %>%
  arrange(order) %>% magrittr::use_series(gene) %>% unique()


ordered <- reorder_df$new+1

avg_goi_M_df <- avg_goi_M$RNA[gene_order,levels(sc.obj.MCells$integrated_snn_res.0.5_named)[ordered]]

vec <- (alpha(c(as.vector(pals::tol(12)), "dark grey"), .8))
names(vec) <- levels(sc.obj.MCells$integrated_snn_res.0.5_named)
vec <- vec[ordered]
annot_col <- list()
annot_col[['cluster']] = vec
ph1 <- pheatmap(avg_goi_M_df, 
                cluster_rows = F,
                annotation_col = data.frame(row.names = levels(sc.obj.MCells$integrated_snn_res.0.5_named), 
                                            cluster = as.factor(levels(sc.obj.MCells$integrated_snn_res.0.5_named))),
                color = rev(pals::brewer.rdbu(100)),
                annotation_colors = annot_col, main = "Clustering Heatmap Myeloid cells",
                show_rownames = T, scale = 'row', cluster_cols = F)


# annotate only selected genes
goi <- M_markers %>% 
  filter(is.na(Column1) == FALSE) %>% 
  magrittr::use_series(gene) %>% unique()
annot_genes_M <-  M_markers %>% 
  dplyr::group_by(cluster)  %>% 
  filter(avg_log2FC>0.6, gene %in% goi, p_val_adj<1e-50) %>% 
  dplyr::arrange(p_val_adj, desc(avg_log2FC)) %>%
  ungroup() %>%
  select(gene,order) %>%
  filter(!duplicated(gene)==TRUE) %>%
  arrange(order) %>% magrittr::use_series(gene) %>% unique()

annot_genes_M_2 <- 
  c("Lyz2", "Sell", "Chil3",
    "Cx3cr1", "Nr4a1","Itgal",
    "Cxcl10","Ccr2",
    "Ly6c2","Vegfa","Il1b",
    "Arg1","Ccl9","Ccl2",
    "Ptgs2", "Klf6", "Tnf",
    "Apoe","Folr2","Mrc1",
    "Fcgr4", "Maf", "Ms4a7",
    "Mgl2", "Mki67",
    "Fcer1g","Pparg","Rac2",
    "Clec10a","Cd209a",
    "Zbtb46","Irf8","Id2",
    "Xcr1","Clec9a", "Btla")


ph1 <- pheatmap(avg_goi_M_df, 
                cluster_rows = F,
                annotation_col = data.frame(row.names = levels(sc.obj.MCells$integrated_snn_res.0.5_named), 
                                            cluster = as.factor(levels(sc.obj.MCells$integrated_snn_res.0.5_named))),
                color = rev(pals::brewer.rdbu(100)),
                annotation_colors = annot_col, main = "Clustering Heatmap Myeloid cells",
                show_rownames = T, show_colnames = F, scale = 'row', cluster_cols = F)





pl1 <- add.flag(ph1, annot_genes_M_2, repel.degree = .1)

color_new_row <- M_markers %>% 
  dplyr::group_by(cluster)  %>% 
  filter(avg_log2FC>0.6, gene %in% goi_unique, p_val_adj<1e-50) %>% 
  dplyr::arrange(p_val_adj, desc(avg_log2FC)) %>%
  ungroup() %>%
  select(gene,order) %>%
  filter(!duplicated(gene)==TRUE) %>%
  arrange(order) %>% left_join(data.frame(names = levels(sc.obj.MCells$integrated_snn_res.0.5_named)[ordered], 
                                          order = 0:12),
                               by = "order") %>%
  left_join(data.frame(vec, names = names(vec)), by = "names") %>%
  magrittr::use_series(vec)

color_new_lines <-  M_markers %>% 
  dplyr::group_by(cluster)  %>% 
  filter(avg_log2FC>0.6, gene %in% goi_unique, p_val_adj<1e-50) %>% 
  filter(gene %in% annot_genes_M_2) %>%
  dplyr::arrange(p_val_adj, desc(avg_log2FC)) %>%
  ungroup() %>%
  select(gene,order) %>%
  filter(!duplicated(gene)==TRUE) %>%
  arrange(order) %>% left_join(data.frame(names = levels(sc.obj.MCells$integrated_snn_res.0.5_named)[ordered], 
                                          order = 0:12),
                               by = "order") %>%
  left_join(data.frame(vec, names = names(vec)), by = "names") %>%
  magrittr::use_series(vec)


pl1$grobs[[1]]$gp=gpar(fontsize=20, fontface = "bold")

pl1$grobs[[3]]$gp=gpar(col=color_new_row, fontsize=20, fontface = "bold")
pl1$grobs[[6]]$children$GRID.text.3038$gp=gpar(fontsize=20, fontface = "bold", vjust = -1)
pl1$grobs[[6]]$children$GRID.text.3038$vjust=0.3
pl1$grobs[[6]]$children$GRID.rect.3039$width=unit(1, "cm")
pl1$grobs[[6]]$children$GRID.rect.3039$height=unit(1, "cm")
pl1$grobs[[6]]$children$GRID.rect.3039$y=rev(unit((1:13)+5, "cm"))
pl1$grobs[[6]]$children$GRID.text.3040$y=rev(unit((1:13)+4.5, "cm"))
pl1$grobs[[6]]$children$GRID.text.3040$x=unit(1.1, "cm")
pl1$grobs[[6]]$children$GRID.text.3040$gp=gpar(fontsize=15, fontface = "bold", vjust = -1)
pl1$grobs[[8]]$x0=unit(0, "cm")
pl1$grobs[[8]]$gp=gpar(col=color_new_lines, lwd=2)
dev.off()
grid.arrange(pl1)



png('Heatmaps_paper_final/MyeloidCells_Heatmap_colored.png', width = 800, height = 800)
grid.arrange(pl1)
dev.off()
pdf('Heatmaps_paper_final/MyeloidCells_Heatmap_colored.pdf', width = 10, height = 10)
grid.arrange(pl1)
dev.off()


avg_goi_M_df2 <- as.data.frame(avg_goi_M_df)
avg_goi_M_df2$clusters <- M_markers %>% 
  dplyr::group_by(cluster)  %>% 
  filter(avg_log2FC>0.6, gene %in% goi_unique, p_val_adj<1e-50) %>% 
  dplyr::arrange(p_val_adj, desc(avg_log2FC)) %>%
  ungroup() %>%
  filter(!duplicated(gene)==TRUE)%>%
  arrange(order) %>% magrittr::use_series(cluster)



write.xlsx(list(`Myeloid Cells Heatmap unscaled` = avg_goi_M_df2, 
                `Myeloid Cells All DEGs` = M_markers[,c(1:7)]), overwrite = T,
           file = 'Heatmaps_paper_final/supp_tables/MyeloidCells_Supp_heatmap.xlsx',
           rowNames = T)

sc.obj.MCells <- UCell::AddModuleScore_UCell(sc.obj.MCells,features = list(DC_new = c("H2-Ab1+", "Itgax+", "Adgre1-")))

table(sc.obj.MCells$DC_new_UCell >0.66)
bc <- names(sc.obj.MCells$DC_new_UCell)[sc.obj.MCells$DC_new_UCell >0.75]
1810/26055

DimPlot(sc.obj.MCells, cells.highlight = bc)



