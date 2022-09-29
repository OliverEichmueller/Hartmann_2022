## Pre-processing objects for RTC4
##------ Mon Oct 18 19:14:43 2021 ------##
## Oliver Eichmueller

# load required libraries ------------------------------------------------------
library(Seurat)
library(dplyr)
source('~/R/Git_Projects/RTC4_2021/Scripts/Parameters.RTC4.processing.R')


# define output directory ------------------------------------------------------
OutDir <- "pre.processed/"
dir.create(OutDir)

# define data directory --------------------------------------------------------
DataDir <- "data/OTP/"
datafiles <- list.files(DataDir)


# create seurat object for data dir --------------------------------------------

for (i in 1:length(datafiles)) {
  cond = datafiles[i]
  
  sc.data = Read10X(
    data.dir = paste0(DataDir, cond, "/filtered_feature_bc_matrix/"))
  
  sc.obj = CreateSeuratObject(sc.data, min.cells = p$min.cells, 
                              min.features = p$min.features,
                              project = cond)
  
  sc.obj[["percent.mt"]] = PercentageFeatureSet(sc.obj, pattern = "^mt\\-")
  
  saveRDS(sc.obj, file = 
            paste0(OutDir, "preprocessed_", 
                   cond,"_min.cells.10.min.features.200.Rds"))
  
}

# getting all metadata for QC plotting -----------------------------------------

preproc.Dir <- "~/R/Git_Projects/RTC4_2021/pre.processed/"

metadata.Out <- "~/R/Git_Projects/RTC4_2021/QC_table/"
dir.create(metadata.Out)

datasets <- list.files(preproc.Dir)
dataset.list <- list()

for (i in 1:length(datasets)) {
  dataset.list[[datasets[i]]] <- readRDS(paste0(preproc.Dir, datasets[i]))
}

nFeature_df_1 <- 
  data.frame(nFeature_RNA = dataset.list[[datasets[1]]]$nFeature_RNA,
             nCount_RNA = dataset.list[[datasets[1]]]$nCount_RNA,
             percent.mt = dataset.list[[datasets[1]]]$percent.mt,
             condition = unique(dataset.list[[datasets[1]]]$orig.ident),
             row.names = NULL)
nFeature_df_2 <- 
  data.frame(nFeature_RNA = dataset.list[[datasets[2]]]$nFeature_RNA,
             nCount_RNA = dataset.list[[datasets[2]]]$nCount_RNA,
             percent.mt = dataset.list[[datasets[2]]]$percent.mt,
             condition = unique(dataset.list[[datasets[2]]]$orig.ident),
             row.names = NULL)
nFeature_df_3 <- 
  data.frame(nFeature_RNA = dataset.list[[datasets[3]]]$nFeature_RNA,
             nCount_RNA = dataset.list[[datasets[3]]]$nCount_RNA,
             percent.mt = dataset.list[[datasets[3]]]$percent.mt,
             condition = unique(dataset.list[[datasets[3]]]$orig.ident),
             row.names = NULL)
nFeature_df_4 <- 
  data.frame(nFeature_RNA = dataset.list[[datasets[4]]]$nFeature_RNA,
             nCount_RNA = dataset.list[[datasets[4]]]$nCount_RNA,
             percent.mt = dataset.list[[datasets[4]]]$percent.mt,
             condition = unique(dataset.list[[datasets[4]]]$orig.ident),
             row.names = NULL)
nFeature_df_5 <- 
  data.frame(nFeature_RNA = dataset.list[[datasets[5]]]$nFeature_RNA,
             nCount_RNA = dataset.list[[datasets[5]]]$nCount_RNA,
             percent.mt = dataset.list[[datasets[5]]]$percent.mt,
             condition = unique(dataset.list[[datasets[5]]]$orig.ident),
             row.names = NULL)

nFeature_df <- 
  rbind(nFeature_df_1, nFeature_df_2, nFeature_df_3, nFeature_df_4, 
        nFeature_df_5)

write.csv(nFeature_df, paste0(metadata.Out, "QC_preprocessed.csv"))
