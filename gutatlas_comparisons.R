#LOADING

library("Seurat")
library("Matrix")
processedgutatlas <- readRDS("processed_gut_atlas_data.RDS")
rawgutatlas <- readRDS("raw_gut_atlas_data.RDS")

##### DIMENSIONS #####

dim(processedgutatlas) # or dim(processedgutatlas@assays$RNA@counts)
dim(rawgutatlas) # or dim(rawgutatlas@assays$RNA@counts)
#3000 fewer genes in the processed - why?

sum(GetAssayData(rawgutatlas)) #2736448287
sum(GetAssayData(processedgutatlas)) #986821345

##### CONFIRM QC OR LACK THEREOF ##### 
#visualize QC metrics: 
head(processedgutatlas@meta.data)
#length(unique(rawgutatlas@meta.data$Diagnosis=="Healthy adult"))

raw_QC <- VlnPlot(rawgutatlas, features=c("nFeature_RNA","nCount_RNA","pct_counts_mt"))
pro_QC <- VlnPlot(processedgutatlas, features=c("nFeature_RNA","nCount_RNA","pct_counts_mt"))

#find max and min Counts for slots: COMPARE
min(rawgutatlas@assays$RNA@counts) # min =0 
max(rawgutatlas@assays$RNA@counts) # max =41487
min(processedgutatlas@assays$RNA@counts) # min =0
max(processedgutatlas@assays$RNA@counts) # max =8.97

min(rawgutatlas@assays$RNA@data) # min =0 
max(rawgutatlas@assays$RNA@data) # max =41487
min(processedgutatlas@assays$RNA@data) # min =0
max(processedgutatlas@assays$RNA@data) # min =8.97

dim(rawgutatlas@assays$RNA@scale.data) #no scale.data
dim(processedgutatlas@assays$RNA@scale.data) #no scale.data

min(rawgutatlas@meta.data$nFeature_RNA) # min = 501
max(rawgutatlas@meta.data$nFeature_RNA) # max = 8236
min(processedgutatlas@meta.data$nFeature_RNA) # min = 501
max(processedgutatlas@meta.data$nFeature_RNA) # max = 7999

min(rawgutatlas@meta.data$nCount_RNA) # min = 608
max(rawgutatlas@meta.data$nCount_RNA) # max = 133696
min(processedgutatlas@meta.data$nCount_RNA) # min = 325.8225 
max(processedgutatlas@meta.data$nCount_RNA) # min = 4772.657 

#checking for doublet filtering - nope
max(processedgutatlas@meta.data$doublet_scores) # max = 0.6518 so it looks like they DIDN'T subset there yet

##### SAMPLE INFO #####
#umbrella cell types - 9
unique(rawgutatlas@meta.data$category)
unique(processedgutatlas@meta.data$category)
#granular cell types - 134
length(unique(rawgutatlas@meta.data$Integrated_05))
length(unique(processedgutatlas@meta.data$Integrated_05)) #134 in each 

#how many donors - ISSUE: WHAT IS SAMPLE.NAME (X2)
length(unique(rawgutatlas@meta.data$donor.name))
length(unique(rawgutatlas@meta.data$sample.name)) # >1 sample from same
length(unique(processedgutatlas@meta.data$donor.name))
length(unique(processedgutatlas@meta.data$sample.name))

#rename Sample.name -> donor.name: 
#colnames(rawgutatlas@meta.data)[colnames(rawgutatlas@meta.data)=="Sample.name"] <- "donor.name"
#colnames(processedgutatlas@meta.data)[colnames(processedgutatlas@meta.data)=="Sample.name"] <- "donor.name"
#head(rawgutatlas@meta.data,5)
#head(processedgutatlas@meta.data,5)

#how many adult donors? 7 -> matches in both datasets now that colnames corrected
length(unique(rawgutatlas@meta.data$donor.name[rawgutatlas@meta.data$Age_group=="Adult"])) #7
length(unique(processedgutatlas@meta.data$donor.name[processedgutatlas@meta.data$Age_group=="Adult"])) #7
length(unique(rawgutatlas@meta.data$donor.name[rawgutatlas@meta.data$Diagnosis=="Healthy adult"])) #7
length(unique(processedgutatlas@meta.data$donor.name[processedgutatlas@meta.data$Diagnosis=="Healthy adult"])) #7


##### FEATURE MISMATCH: 3K FEWER IN PROCESSED ##### 

unique_genes_raw <- setdiff(rownames(rawgutatlas@assays$RNA@data), rownames(processedgutatlas@assays$RNA@data))
unique_genes_raw <- as.data.frame(unique_genes_raw)
write.table(unique_genes_raw,"discarded3k.txt",row.names=F,col.names=F)
test <- read.table("discarded3k.txt")
#followed by this in command line: sed -i 's/"//g' discarded3k.txt

#unique_genes_pro <- setdiff(rownames(processedgutatlas@assays$RNA@data), rownames(rawgutatlas@assays$RNA@data))
#length(unique_genes_raw) #3022
#length(unique_genes_pro) #19 - looks like just numbers appended for genes appearing multiple times?


###### VISUALIZATION: CHANGE IDENTS #####

head(Idents(rawgutatlas),5)
head(Idents(processedgutatlas),5)

Idents(raw_gut_atlas) <- raw_gut_atlas@meta.data$Integrated_05
table(Idents(raw_gut_atlas))


#looking for fibroblast:
library(stringr)
sum(str_detect(raw_gut_atlas@meta.data$Integrated_05, 'fibroblast') >0) #number of cells there granular cell type contains word "fibroblast"

length(unique(raw_gut_atlas@meta.data$Sample.name[raw_gut_atlas@meta.data$Diagnosis=="Healthy adult"])) #7 samples from healthy adults
length(unique(raw_gut_atlas@meta.data$Sample.name[raw_gut_atlas@meta.data$Age_group=="Adult"])) #just confirming; also 7
sum(raw_gut_atlas@meta.data$Diagnosis=="Healthy adult") #and there's 124367 cells originating from those adults total
table(rawgutatlas@meta.data$Sample.name[rawgutatlas@meta.data$Diagnosis=="Healthy adult"]) #per-sample total cell counts

table(raw_gut_atlas@meta.data$Sample.name[raw_gut_atlas@meta.data$Integrated_05=="myofibroblast" & raw_gut_atlas@meta.data$Diagnosis=="Healthy adult"])

sum(raw_gut_atlas@meta.data$Integrated_05=="myofibroblast")

