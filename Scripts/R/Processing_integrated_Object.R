## Processing integrated object 
##------ Mon Oct 28 21:15:55 2021 ------##
## Oliver Eichmueller


# load required libraries ------------------------------------------------------
library(Seurat)
library(dplyr)
require(MarkdownReports)

source('~/R/Git_Projects/RTC4_2021/Scripts/Parameters.RTC4.processing.R')
source('~/Documents/GitHub/Seurat.utils/Functions/Plotting.dim.reduction.3D.R')
source('~/Documents/GitHub/Seurat.utils/00.Load.Seurat.Utils.LOCAL.R')
source("https://raw.githubusercontent.com/vertesy/CodeAndRoll/master/CodeAndRoll.R")

proc.Dir <- "~/R/Git_Projects/RTC4_2021/processed.DS_2/"

sc.obj.integration <- 
  readRDS(file = paste0(OutDir, "unprocessed_RTC4_integration_211028.Rds"))

OutDir <- "~/R/Git_Projects/RTC4_2021/processed.plots_2/"
dir.create(OutDir)
# Perform PCA, calculate UMAP --------------------------------------------------

sc.obj.integration <- ScaleData(sc.obj.integration)
sc.obj.integration <- RunPCA(sc.obj.integration)

# calculate 2D UMAP
sc.obj.integration <- RunUMAP(sc.obj.integration, dims = 1:p$n.PC, 
                              n.components = 2)

# store 2D UMAP in misc
sc.obj.integration@misc$reductions.backup <- list()

sc.obj.integration@misc$reductions.backup[["umap_2D"]] <- 
  sc.obj.integration@reductions[["umap"]]

# calculate 3D UMAP
sc.obj.integration <- RunUMAP(sc.obj.integration, dims = 1:p$n.PC, 
                              n.components = 3)

# store 3D UMAP in misc
sc.obj.integration@misc$reductions.backup[["umap_3D"]] <- 
  sc.obj.integration@reductions[["umap"]]

# Load 2D UMAP for plotting
sc.obj.integration@reductions[["umap"]] <- 
  sc.obj.integration@misc$reductions.backup$umap_2D

DimPlot(sc.obj.integration, group.by = "orig.ident")

# Load 3D UMAP for plotting
sc.obj.integration@reductions[["umap"]] <- 
  sc.obj.integration@misc$reductions.backup$umap_3D

plot3D.umap(obj = sc.obj.integration, category = "orig.ident", AutoAnnotBy = FALSE)


# Calculate clusters -----------------------------------------------------------

sc.obj.integration@reductions[["umap"]] <- 
  sc.obj.integration@misc$reductions.backup$umap_2D

sc.obj.integration <- FindNeighbors(sc.obj.integration, reduction = "pca", dims = 1:p$'n.PC')
sc.obj.integration <- FindClusters(sc.obj.integration, resolution = p$'snn_res')

pl = DimPlot(sc.obj.integration, group.by = "orig.ident")

ggsave(paste0(OutDir, "2D_UMAP_origident.png"), device = "png", plot = pl)


res = stringr::str_extract(colnames(sc.obj.integration@meta.data), "integrated_snn_res.{1,}")
res = res[!is.na(res)]

for (i in 1:length(res)) {
  pl = DimPlot(sc.obj.integration, group.by = res[i], label = T)
  
  ggsave(paste0(OutDir, "2D_UMAP_", res[i], ".png"), device = "png", plot = pl)
}

saveRDS(sc.obj.integration, file = paste0(proc.Dir, "processed_RTC4_integration_211028.Rds"))


# Calculate Cluster markers ----------------------------------------------------
sc.obj.integration@active.assay <- "RNA"
sc.obj.integration@active.ident <- sc.obj.integration$integrated_snn_res.0.2

all_markers <- FindAllMarkers(sc.obj.integration)

goi <- all_markers %>%
  filter(p_val_adj<1e-50, avg_log2FC>0.3) %>%
  magrittr::use_series(gene) %>% unique()

average_goi <- AverageExpression(sc.obj.integration, assays = "RNA", features = goi,
                  group.by = "integrated_snn_res.0.2")

# Export diff expression table

write.csv(all_markers, file = paste0(proc.Dir, "all_markers_res02.csv"), row.names = T)

saveRDS(sc.obj.integration, file = paste0(proc.Dir, "2processed_RTC4_integration_211028.Rds"))
