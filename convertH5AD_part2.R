library("Seurat")
library("SeuratObject")
library("Matrix")
setwd("/workdir/data/crosstissue_data/human/gut/gutcellatlas")

##### load everything and label appropriately for the matrix #####

obj=readMM("gutatlas_processed_v.mtx") #some 20GB, takes a while
obj=t(obj)

obs=read.csv("gutatlas_processed_obs.csv")
var=read.csv("gutatlas_processed_var.csv")

head(obs,5)

rownames(obj)=var$X #"X" is the not-ensemble gene ID format
colnames(obj)=obs$X #sure looks like barcodes

#### make the object, confirm it's real; may take a while ####

processed_gut_atlas=CreateSeuratObject(counts=obj,project="processed_gut_atlas")

head(processed_gut_atlas@meta.data,5)

saveRDS(processed_gut_atlas, file="processed_gut_atlas_data.RDS")


#### ADD THE REST OF THE METADATA FROM OBS ####

#OBSERVATIONS
#for loop is slower, don't do this one
#for (col_name in colnames(obs)) {
#processed_gut_atlas@meta.data[[col_name]] <- obs[[col_name]]}

raw_gut_atlas <- readRDS("raw_gut_atlas_data.RDS")

library(dplyr)

obs=read.csv("gutatlas_raw_obs.csv")

meta_df <- as.data.frame(rawgutatlas@meta.data)
meta_df <- mutate(meta_df, row.names=row.names(meta_df))
head(meta_df,5)
head(obs,5)
#column "X" in obs contains barcodes
meta_df <- left_join(meta_df,obs, by=c("row.names" = "X")) #X in obs, not matched in metadata (rownames) -> name it
rownames(meta_df) <- meta_df$'row.names'
meta_df$'row.names' <- NULL
rawgutatlas@meta.data <- meta_df

head(rawgutatlas@meta.data,5) #this works the same
str(rawgutatlas)

#SAVE THE OBJECT

saveRDS(rawgutatlas, file="raw_gut_atlas_data.RDS")

#uns=read.csv("gutatlas_processed_uns.csv")
#obsm=read.csv("gutatlas_processed_obsm.csv")