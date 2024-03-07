#skeletalsample <- readMM("matrix.mtx")
#dim(skeletalsample)
#var=read.table("features.tsv") #gene IDs
#cosgrovemuscle=read.table("barcodes.tsv.gz")
#dim(var)
#head(cosgrovemuscle)
#rownames(skeletalsample)=var$V2
#colnames(skeletalsample)=cosgrovemuscle$V1
#skeletalsample[0:10,0:3]

library(Seurat)
library(dplyr)
library(readr)
library(Matrix)

#loading and making Seurat

cosgrovemuscle <- read.table("GSE143704_DeMicheli_HumanMuscleAtlas_rawdata.txt.gz", header=T, row.names=1)
#cosgrovemuscle <- as.matrix(cosgrovemuscle)

dim(cosgrovemuscle)
cosgrovemuscle[0:10,1:3]

cosgrove2020_muscle = CreateSeuratObject(counts=cosgrovemuscle,project="cosgrove2020_muscle")
head(cosgrove2020_muscle@meta.data,5)
#head(cosgrove2020_muscle...
cosgrove2020_muscle@assays$RNA@counts[5,5]

#appending metadata

cosgrovemeta <- read.delim("GSE143704_DeMicheli_HumanMuscleAtlas_metadata.txt")
#unique(colnames(cosgrovemeta))
#head(cosgrovemeta$cell_annotation.1)
#head(cosgrovemeta$cell_annotation)

meta_df <- as.data.frame(cosgrove2020_muscle@meta.data)
meta_df <- mutate(meta_df, X=row.names(meta_df))

#head(meta_df,5)
#intersect(meta_df$X, cosgrovemeta$X) #they all match

#column "X" in cosgrovemeta contains barcodes; "X" in the meta.data...
meta_df <- left_join(meta_df,cosgrovemeta, by="X")
rownames(meta_df) <- meta_df$'X'
meta_df$'X' <- NULL #cleaning up
cosgrove2020_muscle@meta.data <- meta_df #back on in

head(cosgrove2020_muscle@meta.data,5)
str(cosgrove2020_muscle@meta.data)

QCplots <- VlnPlot(cosgrove2020_muscle, features=c("nFeature_RNA.x", "nCount_RNA.x", "percent_mito"), ncol=3)

#adding in additional metadata from the series matrix file
#after separating out what i wanted from the series matrix - initially it was all text that messed up read.delim

cosgroveseries <- t(read.delim("seriesmatrixtruncated.txt"))
colnames(cosgroveseries) <- unlist(cosgroveseries[1, ])
cosgroveseries <- cosgroveseries[-1, ]
cosgroveseries <- as.data.frame(cosgroveseries)
cosgroveseries$sampleID <- rownames(cosgroveseries)
rownames(cosgroveseries) <- NULL
dim(cosgroveseries)

meta_df2 <- as.data.frame(cosgrove2020_muscle@meta.data)
meta_df2 <- mutate(meta_df2, X=row.names(meta_df2)) #lost those somehow, protecting them
meta_df2 <- left_join(meta_df2,cosgroveseries, by=c("sampleID" = "sampleID"))
rownames(meta_df2) <- meta_df2$'X'
meta_df2$'X' <- NULL #cleaning up

dim(meta_df)
dim(meta_df2)
str(meta_df)
str(meta_df2)
head(meta_df,2)
head(meta_df2,2)

cosgrove2020_muscle@meta.data <- meta_df2

QCplots <- VlnPlot(cosgrove2020_muscle, features=c("nFeature_RNA.x", "nCount_RNA.x", "percent_mito"), ncol=3) 

saveRDS(cosgrove2020_muscle, file="cosgrove2020_muscle.RDS")