## Sub-Clustering T-cells
##------ Tue Oct 28 13:12:56 2021 ------##
## Oliver Eichmueller

# load required libraries ------------------------------------------------------
library(Seurat)
library(SeuratDisk)
library(dplyr)
require(MarkdownReports)

source('~/R/Git_Projects/RTC4_2021/Scripts/Parameters.RTC4.processing.R')
source('~/Documents/GitHub/Seurat.utils/Functions/Plotting.dim.reduction.3D.R')
source('~/Documents/GitHub/Seurat.utils/00.Load.Seurat.Utils.LOCAL.R')
source("https://raw.githubusercontent.com/vertesy/CodeAndRoll/master/CodeAndRoll.R")

proc.Dir <- "~/R/Git_Projects/RTC4_2021/processed.DS_2/"

sc.obj.integration <- 
  readRDS(file = paste0(proc.Dir, "2processed_RTC4_integration_211028.Rds"))

TCellDir <- "~/R/Git_Projects/RTC4_2021/TCells_new/"
dir.create(TCellDir)

# Version #1: Subcluster using existing CCA and integration --------------------

sc.obj.TCells <- subset(sc.obj.integration, integrated_snn_res.0.2 %in% c(2,5,6,8,10))

# calculate 2D UMAP
sc.obj.TCells <- RunUMAP(sc.obj.TCells, dims = 1:p$n.PC, 
                              n.components = 2)

# store 2D UMAP in misc
sc.obj.TCells@misc$reductions.backup[["umap_2D"]] <- 
  sc.obj.TCells@reductions[["umap"]]

# calculate 3D UMAP
sc.obj.TCells <- RunUMAP(sc.obj.TCells, dims = 1:p$n.PC, 
                              n.components = 3)

# store 3D UMAP in misc
sc.obj.TCells@misc$reductions.backup[["umap_3D"]] <- 
  sc.obj.TCells@reductions[["umap"]]

# Load 2D UMAP for plotting
sc.obj.TCells@reductions[["umap"]] <- 
  sc.obj.TCells@misc$reductions.backup$umap_2D

DimPlot(sc.obj.TCells, group.by = "orig.ident")


# Calculate clusters -----------------------------------------------------------

sc.obj.TCells@reductions[["umap"]] <- 
  sc.obj.TCells@misc$reductions.backup$umap_2D

sc.obj.TCells@active.assay <- "integrated"
sc.obj.TCells <- FindNeighbors(sc.obj.TCells, reduction = "pca", dims = 1:p$'n.PC')
sc.obj.TCells <- FindClusters(object = sc.obj.TCells, resolution = p$'snn_res')
sc.obj.TCells <- FindClusters(object = sc.obj.TCells, resolution = 0.25)

DimPlot(sc.obj.TCells, group.by = "integrated_snn_res.0.25", label = T)

pl = DimPlot(sc.obj.TCells, group.by = "orig.ident")
ggsave(paste0(TCellDir, "2D_TCells_UMAP_origident.png"), device = "png", plot = pl)

pl = DimPlot(sc.obj.TCells, group.by = "integrated_snn_res.0.25", label = T)
ggsave(paste0(TCellDir, "2D_TCells_UMAP_res025.png"), device = "png", plot = pl)

res = stringr::str_extract(colnames(sc.obj.TCells@meta.data), "integrated_snn_res.{1,}")
res = res[!is.na(res)]

for (i in 1:length(res)) {
  pl = DimPlot(sc.obj.TCells, group.by = res[i], label = T)
  
  ggsave(paste0(TCellDir, "2D_UMAP_", res[i], ".png"), device = "png", plot = pl)
}

saveRDS(sc.obj.TCells, file = paste0(proc.Dir, "processed_RTC4_TCells_integration_211028.Rds"))

# Differential expression ------------------------------------------------------

sc.obj.TCells@active.ident <- sc.obj.TCells$integrated_snn_res.0.25
sc.obj.TCells@active.assay <- "RNA"

TCell_markers <- FindAllMarkers(sc.obj.TCells)

goi_TCell <- TCell_markers %>% filter(avg_log2FC>0.3, p_val_adj<1e-100) %>%
  magrittr::use_series(gene) %>% unique()


average_TCell <- AverageExpression(sc.obj.TCells, assays = "RNA", features = goi_TCell,
                                 group.by = "integrated_snn_res.0.3")

pheatmap(average_TCell$RNA, show_rownames = F, scale = "row")


sc.obj.TCells@meta.data
sum_all <- sc.obj.integration@meta.data$orig.ident %>% table() %>% as.data.frame()
colnames(sum_all) = c("orig.ident", "sum")

metadata_TCells <- sc.obj.TCells@meta.data
metadata_TCells %>%
  dplyr::group_by(orig.ident) %>%
  dplyr::count(integrated_snn_res.0.25) %>%
  left_join(sum_all, by = "orig.ident") %>%
  mutate(relative = round(n / sum*100, 2)) %>%
  # filter(integrated_snn_res.0.25 %in% c(2,3,6)) %>%
  ggplot(aes(x = integrated_snn_res.0.25, y = relative, fill = orig.ident)) +
  geom_col(position = "dodge")

metadata_TCells %>%
  dplyr::group_by(orig.ident) %>%
  dplyr::count(integrated_snn_res.0.25) %>%
  mutate(relative = round(n / sum(n)*100, 2)) %>%
  ggplot(aes(x = integrated_snn_res.0.25, y = relative, fill = orig.ident)) +
  geom_col(position = "dodge")


# Export diff expression table

write.csv(TCell_markers, file = paste0(TCellDir, "all_markers_res025_TCells.csv"), row.names = T)
write.csv(TCell_markers %>% filter(p_val_adj < 1e-100), file = paste0(TCellDir, "all_markers_res025_TCells_filtered.csv"), row.names = T)
xlsx::write.xlsx(TCell_markers %>% filter(p_val_adj < 1e-100), file = paste0(TCellDir, "all_markers_res025_TCells_filtered.xlsx"), row.names = T)

FeaturePlot(sc.obj.TCells, features = c("Cd27", "Sell", "Ly6c1", "Il2rb", "Itgal", "Cd44"))
FeaturePlot(sc.obj.TCells, features = c("Havcr2"))
TCell_markers %>% filter(cluster == 0)



# TCells all --------------



sc.obj.TCells@meta.data$UMAP_1 <- 
  sc.obj.TCells@reductions$umap@cell.embeddings[row.names(sc.obj.TCells@meta.data),"UMAP_1"]
sc.obj.TCells@meta.data$UMAP_2 <- 
  sc.obj.TCells@reductions$umap@cell.embeddings[row.names(sc.obj.TCells@meta.data),"UMAP_2"]

loom_subset <- as.loom(sc.obj.TCells, filename = "~/velocyto_test/seurat_TCells_211028.loom")

loom_subset$close_all()




