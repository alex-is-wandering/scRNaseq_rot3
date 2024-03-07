remotes::install_github("mojaveazure/seurat-disk")

install.packages("stringi")
  
library("Seurat")
library("SeuratObject")
library("hdf5r")
library("Matrix")
library("R6")
library("cli")
library("crayon")
library("rlang")
#library("stringi")
library("withr")
library("SeuratDisk")

seurat_obj <- readH5AD("Full_obj_raw_counts_nosoupx_v2.h5ad")
saveRDS(seurat_obj, file="raw_gut_atlas_data.RDS")