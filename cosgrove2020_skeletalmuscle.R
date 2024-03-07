library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(ComplexHeatmap)
library(Matrix)

##### read it in ####

cosgrove2020_muscle <- readRDS("cosgrove2020_muscle.RDS")

unique(colnames(cosgrove2020_muscle@meta.data))
dim(cosgrove2020_muscle@assays$RNA@data)

QCplots <- VlnPlot(cosgrove2020_muscle, features=c("nFeature_RNA.x", "nCount_RNA.x", "percent_mito"), ncol=3)

#### _ ####

