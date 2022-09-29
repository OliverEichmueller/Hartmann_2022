## Sub-Clustering Macrophages and DCs on mito filtered dataset
##------ Fri Oct 29 07:31:15 2021 ------##
## Oliver Eichmueller


# load required libraries ------------------------------------------------------
library(Seurat)
library(dplyr)
require(MarkdownReports)
require(ggplot2)

source('~/R/Git_Projects/RTC4_2021/Scripts/Parameters.RTC4.processing.R')
source('~/Documents/GitHub/Seurat.utils/Functions/Plotting.dim.reduction.3D.R')
source('~/Documents/GitHub/Seurat.utils/00.Load.Seurat.Utils.LOCAL.R')
source("https://raw.githubusercontent.com/vertesy/CodeAndRoll/master/CodeAndRoll.R")

proc.Dir <- "~/R/Git_Projects/RTC4_2021/processed.DS_2/"

sc.obj.integration <- 
  readRDS(file = paste0(proc.Dir, "2processed_RTC4_integration_211028.Rds"))

MCellDir <- "~/R/Git_Projects/RTC4_2021/Macro_DC_2/"
dir.create(MCellDir)

# Version #1: Subcluster using existing CCA and integration --------------------
DimPlot(sc.obj.integration, group.by = "integrated_snn_res.0.2", label = T)
sc.obj.MCells <- subset(sc.obj.integration, integrated_snn_res.0.2 %in% c(0,1,3,4,9))

# calculate 2D UMAP
sc.obj.MCells <- RunUMAP(sc.obj.MCells, dims = 1:p$n.PC, 
                         n.components = 2)

# store 2D UMAP in misc
sc.obj.MCells@misc$reductions.backup[["umap_2D"]] <- 
  sc.obj.MCells@reductions[["umap"]]

# calculate 3D UMAP
sc.obj.MCells <- RunUMAP(sc.obj.MCells, dims = 1:p$n.PC, 
                         n.components = 3)

# store 3D UMAP in misc
sc.obj.MCells@misc$reductions.backup[["umap_3D"]] <- 
  sc.obj.MCells@reductions[["umap"]]

# Load 2D UMAP for plotting
sc.obj.MCells@reductions[["umap"]] <- 
  sc.obj.MCells@misc$reductions.backup$umap_2D


DimPlot(sc.obj.MCells, group.by = "orig.ident")


# Calculate clusters -----------------------------------------------------------

sc.obj.MCells@reductions[["umap"]] <- 
  sc.obj.MCells@misc$reductions.backup$umap_2D

sc.obj.MCells@active.assay <- "integrated"
sc.obj.MCells <- FindNeighbors(sc.obj.MCells, reduction = "pca", dims = 1:p$'n.PC')
sc.obj.MCells <- FindClusters(object = sc.obj.MCells, resolution = p$'snn_res')
sc.obj.MCells <- FindClusters(object = sc.obj.MCells, resolution = 0.25)

pl = DimPlot(sc.obj.MCells, group.by = "orig.ident")
ggsave(paste0(MCellDir, "2D_MCells_UMAP_origident.png"), device = "png", plot = pl)


res = stringr::str_extract(colnames(sc.obj.MCells@meta.data), "integrated_snn_res.{1,}")
res = res[!is.na(res)]

for (i in 1:length(res)) {
  pl = DimPlot(sc.obj.MCells, group.by = res[i], label = T)
  
  ggsave(paste0(MCellDir, "2D_UMAP_", res[i], ".png"), device = "png", plot = pl)
}

saveRDS(sc.obj.MCells, file = paste0(proc.Dir, "processed_RTC4_MCells_integration_211028.Rds"))

# sc.obj.MCells <- readRDS(paste0(proc.Dir, "processed_RTC4_MCells_integration_211028.Rds"))

# Compare groups ---------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(UCell)

pl2 <- DimPlot(sc.obj.MCells, group.by = 'integrated_snn_res.0.3', label = T)

pl1 <- sc.obj.MCells@meta.data %>%
  dplyr::select(orig.ident, integrated_snn_res.0.3) %>%
  dplyr::group_by(orig.ident) %>%
  dplyr::count(integrated_snn_res.0.3) %>%
  left_join(sc.obj.integration@meta.data %>% 
              dplyr::count(orig.ident), by = "orig.ident") %>%
  mutate(rel = round(n.x/n.y*100,2)) %>%
  ggplot(aes(x=orig.ident, y=rel, fill = orig.ident)) +
  geom_col() + facet_wrap(~integrated_snn_res.0.3, ncol = 10) +
  theme(axis.text.x = element_text(angle=90,hjust=01, vjust=.5)) + ylab('Pct of whole dataset') + xlab('Condition')
pl1+pl2
sc.obj.MCells@active.ident <- sc.obj.MCells$integrated_snn_res.0.3
allmarkers_myeloid <- FindAllMarkers(sc.obj.MCells, only.pos = T)

allmarkers_myeloid %>% filter(p_val_adj<1e-50) %>%
  dplyr::arrange(cluster, desc(p_val_adj))
allmarkers_myeloid_filtered <- allmarkers_myeloid %>% filter(p_val_adj<1e-50) %>%
    dplyr::arrange(cluster, (p_val_adj))
write.csv(allmarkers_myeloid_filtered, 'Macro_DC_2/DiffExp_res03_filtered.csv')

allmarkers_myeloid %>% filter(gene %in%  (markers %>% filter(annot == "M2") %>% magrittr::use_series(marker))) %>%
  filter(p_val_adj <1e-100)

library(STRINGdb)
string_db <- STRINGdb$new(species = 10090)

example_string <- string_db$map(allmarkers_myeloid %>% filter(p_val_adj<1e-100) %>%
                                  filter(cluster==0),
                                "gene")
string_hits <- example_string %>% magrittr::use_series(STRING_id)


string_db$plot_network(string_hits[!is.na(string_hits)])
enrichment <- string_db$get_enrichment(string_hits[!is.na(string_hits)])

enrichment$category %>% table
enrichment %>%
  dplyr::select(category, term, description, p_value, fdr, number_of_genes) %>%
  filter(category %in% "Process") %>%
  arrange((p_value)) %>% dplyr::slice(1:20) %>%
  grid.table()
dev.off()
FeaturePlot(sc.obj.MCells, features = c("Cxcl3"), col = c(alpha("grey", .5),"red"))

saveRDS(sc.obj.MCells, 'processed.DS_2/processed_RTC4_MCells_integration_211029.Rds')





# MCells all --------------



sc.obj.MCells@meta.data$UMAP_1 <- 
  sc.obj.MCells@reductions$umap@cell.embeddings[row.names(sc.obj.MCells@meta.data),"UMAP_1"]
sc.obj.MCells@meta.data$UMAP_2 <- 
  sc.obj.MCells@reductions$umap@cell.embeddings[row.names(sc.obj.MCells@meta.data),"UMAP_2"]

loom_subset <- as.loom(sc.obj.MCells, filename = "~/velocyto_test/seurat_MCells_211028.loom")

loom_subset$close_all()

