## Integrating objects for RTC4
##------ Mon Oct 18 21:15:55 2021 ------##
## Oliver Eichmueller

# load required libraries ------------------------------------------------------
library(Seurat)
library(dplyr)
source('~/R/Git_Projects/RTC4_2021/Scripts/Parameters.RTC4.processing.R')


# Integrate 10X data ===========================================================
# this was performed to combine the datasets

# read in pre-processed objects ------------------------------------------------
preproc.Dir <- "~/R/Git_Projects/RTC4_2021/pre.processed/"

proc.Dir <- "~/R/Git_Projects/RTC4_2021/processed.DS_2/"
dir.create(proc.Dir)

datasets <- list.files(preproc.Dir, pattern = "preprocessed")
object.list <- list()

for (i in 1:length(datasets)) {
  object.list[[datasets[i]]] <- readRDS(paste0(preproc.Dir, datasets[i]))
}


# preprocess for integration ---------------------------------------------------
# store all datasets in a list

# re-do pre-processing
for (i in 1: length(object.list)) {
  object.list[[i]][["percent.mt"]] <- 
    PercentageFeatureSet(object.list[[i]], pattern = "^mt\\-")
  
}

# subset based on parameters of p
for (i in 1: length(object.list)) {
  object.list[[i]] <- subset(object.list[[i]], 
                             subset = nFeature_RNA > p$thr.hp.nFeature_RNA)
  object.list[[i]] <- subset(object.list[[i]], 
                             subset = percent.mt < p$thr.lp.mito)
  
}

# normalise and calculate variable features
for (i in 1: length(object.list)) {
  object.list[[i]] <- NormalizeData(object.list[[i]], 
                                    normalization.method = "LogNormalize", 
                                    scale.factor = 10000)
  
  object.list[[i]] <- FindVariableFeatures(object.list[[i]], 
                                           selection.method = "vst", 
                                           nfeatures = 2000)
}


# perform integration ----------------------------------------------------------
object.list_int_anchors <- FindIntegrationAnchors(object.list = object.list, 
                                                  dims = 1:p$n.PC)
saveRDS(object.list_int_anchors, 
        file = paste0(proc.Dir, "integration_anchors_211028.Rds"))

sc.obj.integration <- IntegrateData(anchorset = object.list_int_anchors, 
                                    dims = 1:p$n.PC)

saveRDS(sc.obj.integration, 
        file = paste0(proc.Dir, "unprocessed_RTC4_integration_211028.Rds"))
