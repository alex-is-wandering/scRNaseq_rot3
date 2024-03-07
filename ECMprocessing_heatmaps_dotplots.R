library(Seurat)
library(dplyr)
library(patchwork)
library(ComplexHeatmap)
library(ggplot2)

###### READING IN DATA ###### 
#*pre-processed data from AS's bioxriv paper

#loading matrix
#/local1/workdir/agh227/423-Vastus83M/filtered_feature_bc_matrix

#immune<-readRDS("2021-01-22.human_baseline_immune_v3.RDS")
parenchymal<-readRDS("2020-06-07.human_baseline_paren_v3.RDS")
#endothelial<-readRDS("2021-01-22.human_baseline_endo_v3.RDS")
#endo_granular<-readRDS("2021-01-22.human_baseline_endo_granular_v3.RDS")

#ENDOTHELIAL
#*extract metadata from immune - nephrectomy & age column - and map to orig.ident

#sorting the endothelial data by region of the kidney
endo_cellnames <- colnames(endothelial)

endothelial$cortex <- grepl("cortex",endo_cellnames,ignore.case=T)
endothelial$medulla <- grepl("medulla", endo_cellnames,ignore.case=T)
endothelial$hilum <-grepl("pelvis", endo_cellnames,ignore.case=T)

endo_cortex_subset <- subset(endothelial, subset=cortex==TRUE)
endo_medulla_subset <- subset(endothelial, subset=medulla==TRUE)
endo_hilum_subset <- subset(endothelial, subset=hilum==TRUE)

DimPlot(endo_cortex_subset,reduction="umap")
DimPlot(endo_medulla_subset,reduction="umap")
DimPlot(endo_hilum_subset,reduction="umap")

#PARENCHYMAL
head(parenchymal@meta.data,5) 
#set ident to celltype, for granular vs high-level
Idents(object = parenchymal) <- parenchymal$celltype
DimPlot(parenchymal,reduction="umap")

paren_cortex_subset <- subset(parenchymal, subset=location=="C")
paren_medulla_subset <- subset(parenchymal, subset=location=="M")
paren_hilum_subset <- subset(parenchymal, subset=location=="P")

DimPlot(paren_cortex_subset,reduction="umap")
DimPlot(paren_medulla_subset,reduction="umap")
DimPlot(paren_hilum_subset,reduction="umap")

#IMMUNE

imm_cortex_subset <- subset(immune, subset=location=="C")
imm_medulla_subset <- subset(immune, subset=location=="M")
imm_hilum_subset <- subset(immune, subset=location=="P")

DimPlot(imm_cortex_subset,reduction="umap",group.by="celltype")
DimPlot(imm_medulla_subset,reduction="umap",group.by="celltype")
DimPlot(imm_hilum_subset,reduction="umap",group.by="celltype")


###### ECM EXPLORATION ###### 

parenchymal<-readRDS("2020-06-07.human_baseline_paren_v3.RDS")
matrisomelist <- read.delim("/local/workdir/agh227/Hs_Matrisome_Masterlist_Naba et al_2012.xlsx - Hs_Matrisome_Masterlist.tsv")

paren_geneIDs <- rownames(parenchymal@assays$RNA@data)
COL_genes <- paren_geneIDs[grepl("COL",paren_geneIDs,ignore.case=FALSE)]

###### looking at celltypes vs individual collagen subunits ###### 
#come back to later
#VIOLIN PLOT PER-GENE
VlnPlot(parenchymal,features="COL10A1",slot="data")

#DOT PLOTS FOR CATEGORIES, ILLUSTRATING PER-GENE
categories = c("Collagens")
selected_genes <- matrisomelist$Gene.Symbol[matrisomelist$Matrisome.Category %in% categories]
#*we only want to map those genes that are named in the matrisomelist
selected_genes <- selected_genes[selected_genes %in% rownames(parenchymal@assays$RNA@data)]

DotPlot(parenchymal,features=selected_genes) + theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave("/workdir/agh227/kidneydata_figures/dotplot_alltypes_collagens.png",width=20,height=10)
dev.off()

###### SCORING: looking at celltypes vs ECM component categories ######

# Loop through each unique category in the gene set
for (Matrisome.Category in unique(matrisomelist$Matrisome.Category)) {
  
  # Remove spaces from the category name
  cleaned_category <- gsub(" ", "_", gsub("-", "_", Matrisome.Category))
  
  selected_genes <- matrisomelist$Gene.Symbol[matrisomelist$Matrisome.Category == Matrisome.Category]
  
  parenchymal <- AddModuleScore(parenchymal, features = list(selected_genes), name = Matrisome.Category)
  
  }

# next:

cat_counts <- table(matrisomelist$Matrisome.Category)

head(parenchymal@meta.data,5)

unique_categories <- unique(matrisomelist$Matrisome.Category)  

paren_Glycopro <- VlnPlot(parenchymal, features="ECM_Glycoproteins1")
paren_Proteogl <- VlnPlot(parenchymal, features="Proteoglycans1")
paren_Collagen <- VlnPlot(parenchymal, features="Collagens1")
paren_Affiliat <- VlnPlot(parenchymal, features="ECM_affiliated_Proteins1")
paren_Secreted <- VlnPlot(parenchymal, features="Secreted_Factors1")
paren_Regulato <- VlnPlot(parenchymal, features="ECM_Regulators1")

###### ComplexHeatmap: gene expression vs ECM categories ###### 
unique_celltypes <- unique(parenchymal$celltype)
unique_categories <- unique(matrisomelist$Matrisome.Category)

matrisomelist <- read.delim("/local/workdir/agh227/Hs_Matrisome_Masterlist_Naba et al_2012.xlsx - Hs_Matrisome_Masterlist.tsv")

categories = c("ECM Regulators")
#matrisome-affiliated: "ECM-affiliated Proteins","ECM Regulators","Secreted Factors"
#core matrisome: "Collagens","ECM-affiliated Proteins","ECM Glycoproteins"

selected_genes <- matrisomelist$Gene.Symbol[matrisomelist$Matrisome.Category %in% categories]
#*we only want to map those genes that are named in the matrisomelist
selected_genes <- selected_genes[selected_genes %in% rownames(parenchymal@assays$RNA@data)]




###### Seurat DoHeatmap: gene expression vs ECM categories ###### 
#GetAssayData(parenchymal,slot='scale.data')[1:10,1:5]
#expression data: parenchymal@assays$RNA@scale.data

unique_celltypes <- unique(parenchymal$celltype)
unique_categories <- unique(matrisomelist$Matrisome.Category)

#if needed, reload matrisome reference
matrisomelist <- read.delim("/local/workdir/agh227/Hs_Matrisome_Masterlist_Naba et al_2012.xlsx - Hs_Matrisome_Masterlist.tsv")


#*naming which categories to map
categories = c("Collagens","ECM-affiliated Proteins","ECM Glycoproteins")
#matrisome-affiliated: "ECM-affiliated Proteins","ECM Regulators","Secreted Factors"
#core matrisome: "Collagens","ECM-affiliated Proteins","ECM Glycoproteins"

#*specifying geneID location
selected_genes <- matrisomelist$Gene.Symbol[matrisomelist$Matrisome.Category %in% categories]
#*we only want to map those genes that are named in the matrisomelist
selected_genes <- selected_genes[selected_genes %in% rownames(parenchymal@assays$RNA@data)]

#*making heatmaps
heatmap <- DoHeatmap(parenchymal, features=selected_genes, slot="data")
ggsave("/workdir/agh227/kidneydata_figures/heatmap_alltypes_corematrisomea.png",width=20,height=10)
dev.off()

###*what if we only want to look at fibroblasts?
selected_celltypes = c("Fibroblasts-1","Fibroblasts-2","Fibroblasts-3")
selected_celltypes <- selected_celltypes[selected_celltypes %in% parenchymal@meta.data$celltype]

subset <- subset(parenchymal, subset = celltype %in% selected_celltypes)

#*making heatmaps
heatmap <- DoHeatmap(subset,features=selected_genes, slot="data")
ggsave("/workdir/agh227/kidneydata_figures/heatmap_fibroblasts_corematrisome.png",width=20,height=10)
dev.off()