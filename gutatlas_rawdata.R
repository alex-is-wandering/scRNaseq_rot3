library("Seurat")
library("Matrix")
library("ggplot2")
setwd("/workdir/data/crosstissue_data/human/gut/elmentaite2021_gutcellatlas")
raw_gut_atlas <- readRDS("raw_gut_atlas_data.RDS")
#rawgut_healthyadult_mesenchymal <- readRDS("/workdir/data/crosstissue_data/human/kidney/subramanian2021/2020-06-07.human_baseline_paren_v3.RDS")

##### checking for annotations and info #####
head(raw_gut_atlas@meta.data,5)

#get numbers of cells per-cell type
table(Idents(raw_gut_atlas))

#categorical info: what metadata exists
colnames(raw_gut_atlas@meta.data)
#umbrella cell types - 9
unique(raw_gut_atlas@meta.data$category)
unique(raw_gut_atlas@meta.data$Region)
#granular cell types - 144
unique(raw_gut_atlas@meta.data$Integrated_05)
#looking for fibroblast:
library(stringr)
sum(str_detect(raw_gut_atlas@meta.data$Integrated_05, 'fibroblast') >0) #number of cells in granular cell type list that contain word "fibroblast"

#subsets by age - the only adult is "healthy adult"
#unique(raw_gut_atlas@meta.data$Diagnosis)
sum(unique(raw_gut_atlas@meta.data$Sample.name))

##### process #####
#*QC metrics & filtering based on them
raw_QC <- VlnPlot(raw_gut_atlas, features=c("nFeature_RNA","nCount_RNA","pct_counts_mt"))
#*cells have already been filtered based on %mito (<50%) and nFeature (>500)
plot1 <- FeatureScatter(raw_gut_atlas,feature1="nCount_RNA",feature2="pct_counts_mt") #=-0.06
plot2 <- FeatureScatter(raw_gut_atlas,feature1="nCount_RNA",feature2="nFeature_RNA") #=0.89
correlates <- plot1+plot2
ggsave(filename = "rawQCcorrelations.png",plot=correlates)
#*as in AS's paper, don't filter based on an upper limit of nFeatures

#*NORMALIZING:
raw_gut_atlas <- NormalizeData(raw_gut_atlas)

##### finding adult group & isolating #####
#*confirming coherency between labels:
length(unique(raw_gut_atlas@meta.data$Sample.name[raw_gut_atlas@meta.data$Diagnosis=="Healthy adult"])) #7 samples from healthy adults
length(unique(raw_gut_atlas@meta.data$Sample.name[raw_gut_atlas@meta.data$Age_group=="Adult"])) #just confirming; also 7
#sum(raw_gut_atlas@meta.data$Diagnosis=="Healthy adult") #and there's 124367 cells originating from those adults total
table(raw_gut_atlas@meta.data$Sample.name[raw_gut_atlas@meta.data$Diagnosis=="Healthy adult"]) #per-sample total cell counts
table(raw_gut_atlas@meta.data$Sample.name[raw_gut_atlas@meta.data$Age_group=="Adult"]) #per-sample total cell counts

#goal: get a small table with proportions of cell types per donor
table(raw_gut_atlas@meta.data$Sample.name[raw_gut_atlas@meta.data$Integrated_05=="myofibroblast" & raw_gut_atlas@meta.data$Diagnosis=="Healthy adult"]) #that's pretty stark - why does A32 have almost all of them?

##### SUBSET for healthy adult; going off the difference in cell# and diagnosis being a more embedded column ####
rawgut_healthyadult <- subset(raw_gut_atlas, subset=Diagnosis=="Healthy adult")
head(rawgut_healthyadult@meta.data)

saveRDS(rawgut_healthyadult, "rawgutatlas_adultsubset.RDS")

rawgut_healthyadult <- readRDS("rawgutatlas_adultsubset.RDS")

##### ECM association #####
matrisomelist <- read.delim("/local/workdir/agh227/reference_ECM_signatures/Hs_Matrisome_Masterlist_Naba et al_2012.xlsx - Hs_Matrisome_Masterlist.tsv")
reactome_ECMOrg <- read.delim("/workdir/agh227/reference_ECM_signatures/reactome_ECMOrganization.tsv")
#genes at "reactome_ECMOrg$associations.gene.symbol"

#set Idents to relevant cell types so they group right when graphing
#colnames(rawgut_healthyadult@meta.data)
#unique(rawgut_healthyadult@meta.data$category)
Idents(rawgut_healthyadult)
#Idents(rawgut_healthyadult) <- raw_gut_atlas@meta.data$Integrated_05
#Idents(rawgut_healthyadult) <- rawgut_healthyadult@meta.data$category

#COL1A1 as a proxy check
Idents(rawgut_healthyadult) <- rawgut_healthyadult@meta.data$category
VlnPlot(rawgut_healthyadult,features="COL1A1",slot="data") #by larger category, nothing much

#cell type specific, SUBSETTING CELL TYPES OF INTEREST:
rawgut_healthyadult_mesenchymal <- subset(raw_gut_atlas, subset=category=="Mesenchymal")
Idents(rawgut_healthyadult_mesenchymal) <- rawgut_healthyadult_mesenchymal@meta.data$Integrated_05
VlnPlot(rawgut_healthyadult_mesenchymal,features="COL1A1",slot="data") #by larger category, nothing much

##### REACTOME: ECM ORG #####

#looking at collagens:
mesen_geneIDs <- rownames(rawgut_healthyadult_mesenchymal@assays$RNA@data)
COL_genes <- mesen_geneIDs[grepl("COL",mesen_geneIDs,ignore.case=FALSE)]
selected_genes <- reactome_ECMOrg$associations.gene.symbol
selected_genes <- selected_genes[selected_genes %in% COL_genes]
dotplot_mesenchymal_collagens <- DotPlot(rawgut_healthyadult_mesenchymal,features=selected_genes) + theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave("/workdir/agh227/elmentaite2021_gutatlas/Reactome_ECMOrg/dotplot_mesenchymal_collagens.png",plot=dotplot_mesenchymal_collagens,width=16,height=7)
dev.off()

#getting genes in common, module scoring the term against the healthy adult/mesenchymal subset:
selected_genes <- reactome_ECMOrg$associations.gene.symbol
selected_genes <- selected_genes[selected_genes %in% rownames(rawgut_healthyadult_mesenchymal@assays$RNA@data)] 

#module scoring for the functional term:
rawgut_healthyadult_mesenchymal <- AddModuleScore(rawgut_healthyadult_mesenchymal, features = list(selected_genes), name = "ECMOrg")
VlnPlot(rawgut_healthyadult_mesenchymal, features="ECMOrg1")
ggsave("/workdir/agh227/elmentaite2021_gutatlas/Reactome_ECMOrg/mesenchymal_ViolinECMOrg_AddModuleScored.png",width=20,height=10)

###### Seurat DoHeatmap: gene expression vs ECM organization ###### 
#*naming which categories to map
#*specifying geneID location
selected_genes <- reactome_ECMOrg$associations.gene.symbol
selected_genes <- selected_genes[selected_genes %in% rownames(rawgut_healthyadult_mesenchymal@assays$RNA@data)]

#*making heatmaps: ALL CORE MATRISOME COMPONENTS (3) categories, all cell types
DoHeatmap(rawgut_healthyadult_mesenchymal, features=selected_genes, slot="data")
ggsave("/workdir/agh227/elmentaite2021_gutatlas/Reactome_ECMOrg/mesen_heatmap_ECMOrg.png",width=20,height=10)

dev.off()

##### MATRISOME: CATEGORIES #####

for (Matrisome.Category in unique(matrisomelist$Matrisome.Category)) {
  
  # Remove spaces from the category name
  cleaned_category <- gsub(" ", "_", gsub("-", "_", Matrisome.Category))
  
  selected_genes <- matrisomelist$Gene.Symbol[matrisomelist$Matrisome.Category == Matrisome.Category]
  
  rawgut_healthyadult_mesenchymal <- AddModuleScore(rawgut_healthyadult_mesenchymal, features = list(selected_genes), name = cleaned_category) #uses "data" unless "scale=TRUE"
  
}

unique_categories <- unique(matrisomelist$Matrisome.Category)  
VlnPlot(rawgut_healthyadult_mesenchymal, features="ECM_Glycoproteins1")
ggsave("/workdir/agh227/elmentaite2021_gutatlas/Matrisome/mesen_ViolinGlycoproteins_AddModuleScored.png",width=20,height=10)
dev.off()

VlnPlot(rawgut_healthyadult_mesenchymal, features="Proteoglycans1")
ggsave("/workdir/agh227/elmentaite2021_gutatlas/Matrisome/mesen_ViolinProteoglycans_AddModuleScored.png",width=20,height=10)
dev.off()

VlnPlot(rawgut_healthyadult_mesenchymal, features="Collagens1")
ggsave("/workdir/agh227/elmentaite2021_gutatlas/Matrisome/mesen_ViolinCollagens_AddModuleScored.png",width=20,height=10)
dev.off()

VlnPlot(rawgut_healthyadult_mesenchymal, features="ECM_affiliated_Proteins1")
ggsave("/workdir/agh227/elmentaite2021_gutatlas/Matrisome/mesen_ViolinECM_affiliated_AddModuleScored.png",width=20,height=10)
dev.off()

VlnPlot(rawgut_healthyadult_mesenchymal, features="Secreted_Factors1")
ggsave("/workdir/agh227/elmentaite2021_gutatlas/Matrisome/mesen_ViolinSecretedFactors_AddModuleScored.png",width=20,height=10)
dev.off()

VlnPlot(rawgut_healthyadult_mesenchymal, features="ECM_Regulators1")
ggsave("/workdir/agh227/elmentaite2021_gutatlas/Matrisome/mesen_ViolinECMRegulators_AddModuleScored.png",width=20,height=10)
dev.off()

##### checking emptiness of matrix #####

#dim(raw_gut_atlas@assays$RNA@counts)
#pulling specific entries
checkmatrix <- raw_gut_atlas@assays$RNA@counts[33530:33538, 1930:1945]
View(checkmatrix) #Error in validityMethod(as(object, superClass)) : object 'Csparse_validate' not found -- can't use "print" still. with checkmatrix@i and @p I can tell that this contains 0s and 1s

#class(GetAssayData(raw_gut_atlas))
sum(GetAssayData(raw_gut_atlas)) #sum is nonzero, matrix isn't empty
head(raw_gut_atlas@assays$RNA@counts) #doesn't work
dim(raw_gut_atlas@assays$RNA@counts)
