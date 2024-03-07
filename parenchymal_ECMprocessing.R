library(Seurat)
library(dplyr)
library(patchwork)
library(ComplexHeatmap)
library(ggplot2)

###### READING IN DATA ###### 
setwd("/workdir/data/crosstissue_data/kidney")
parenchymal<-readRDS("2020-06-07.human_baseline_paren_v3.RDS")
matrisomelist <- read.delim("/local/workdir/agh227/reference_ECM_signatures/Hs_Matrisome_Masterlist_Naba et al_2012.xlsx - Hs_Matrisome_Masterlist.tsv")

#head(parenchymal@meta.data,5) 
#set ident to celltype, for granular vs high-level
Idents(object = parenchymal) <- parenchymal$celltype
DimPlot(parenchymal,reduction="umap")
ggsave("/workdir/agh227/kidneydata_figures/parenchymal_UMAP.png",width=10,height=10)
dev.off()

#cell origin: region of the kidney
paren_cortex_subset <- subset(parenchymal, subset=location=="C")
paren_medulla_subset <- subset(parenchymal, subset=location=="M")
paren_hilum_subset <- subset(parenchymal, subset=location=="P")

###### ECM EXPLORATION ###### 

paren_geneIDs <- rownames(parenchymal@assays$RNA@data) #not scale.data
COL_genes <- paren_geneIDs[grepl("COL",paren_geneIDs,ignore.case=FALSE)]

###### dotplot & vlnplot: looking at celltypes vs individual collagen subunits ###### 
#VIOLIN PLOT PER-GENE
VlnPlot(parenchymal,features="COL10A1",slot="data")
VlnPlot(parenchymal,features="COL1A1",slot="data")
ggsave("/workdir/agh227/kidneydata_figures/paren_COL1A1_violin.png",width=10,height=5)
dev.off()

#DOT PLOTS FOR CATEGORIES, ILLUSTRATING PER-GENE
categories = c("Collagens")
selected_genes <- matrisomelist$Gene.Symbol[matrisomelist$Matrisome.Category %in% categories]
#*we only want to map those genes that are named in the matrisomelist
selected_genes <- selected_genes[selected_genes %in% rownames(parenchymal@assays$RNA@data)]

DotPlot(parenchymal,features=selected_genes) + theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave("/workdir/agh227/kidneydata_figures/paren_dotplot_alltypes_collagens.png",width=25,height=10) #not sure why background is wonky, but it works
dev.off()

###### Seurat AddModuleScore: looking at celltypes vs ECM component categories ######

# Loop through each unique category in the gene set - looking @ all the data, not just scale.data
for (Matrisome.Category in unique(matrisomelist$Matrisome.Category)) {
  
  # Remove spaces from the category name
  cleaned_category <- gsub(" ", "_", gsub("-", "_", Matrisome.Category))

  selected_genes <- matrisomelist$Gene.Symbol[matrisomelist$Matrisome.Category == Matrisome.Category]
  
  parenchymal <- AddModuleScore(parenchymal, features = list(selected_genes), name = cleaned_category) #uses "data" unless "scale=TRUE"
  
}

saveRDS(parenchymal, file="/workdir/agh227/processeddata/parenchymal_matrisomeScored.RDS")
head(parenchymal@meta.data,5)

# next:

cat_counts <- table(matrisomelist$Matrisome.Category)

head(parenchymal@meta.data,5)

unique_categories <- unique(matrisomelist$Matrisome.Category)  

VlnPlot(parenchymal, features="ECM_Glycoproteins1")
ggsave("/workdir/agh227/kidneydata_figures/paren_ViolinGlycoproteins_AddModuleScored.png",width=20,height=10)
dev.off()

VlnPlot(parenchymal, features="Proteoglycans1")
ggsave("/workdir/agh227/kidneydata_figures/paren_ViolinProteoglycans_AddModuleScored.png",width=20,height=10)
dev.off()

VlnPlot(parenchymal, features="Collagens1")
ggsave("/workdir/agh227/kidneydata_figures/paren_ViolinCollagens_AddModuleScored.png",width=20,height=10)
dev.off()

VlnPlot(parenchymal, features="ECM_affiliated_Proteins1")
ggsave("/workdir/agh227/kidneydata_figures/paren_ViolinECM_affiliated_AddModuleScored.png",width=20,height=10)
dev.off()

VlnPlot(parenchymal, features="Secreted_Factors1")
ggsave("/workdir/agh227/kidneydata_figures/paren_ViolinSecretedFactors_AddModuleScored.png",width=20,height=10)
dev.off()

VlnPlot(parenchymal, features="ECM_Regulators1")
ggsave("/workdir/agh227/kidneydata_figures/paren_ViolinECMRegulators_AddModuleScored.png",width=20,height=10)
dev.off()

###### Seurat DoHeatmap: gene expression vs ECM categories ###### 
#GetAssayData(parenchymal,slot='scale.data')[1:10,1:5]
#expression data: parenchymal@assays$RNA@scale.data

unique_celltypes <- unique(parenchymal$celltype)
unique_categories <- unique(matrisomelist$Matrisome.Category)

#*naming which categories to map
categories = c("Collagens","Proteoglycans","ECM Glycoproteins")
#matrisome-affiliated: "ECM-affiliated Proteins","ECM Regulators","Secreted Factors"
#core matrisome: "Collagens","Proteoglycans","ECM Glycoproteins"
#*specifying geneID location
selected_genes <- matrisomelist$Gene.Symbol[matrisomelist$Matrisome.Category %in% categories]
#*we only want to map those genes that are named in the matrisomelist
selected_genes <- selected_genes[selected_genes %in% rownames(parenchymal@assays$RNA@data)]

#*making heatmaps: ALL CORE MATRISOME COMPONENTS (3) categories, all cell types
DoHeatmap(parenchymal, features=selected_genes, slot="data")
ggsave("/workdir/agh227/kidneydata_figures/paren_heatmap_alltypes_corematrisome.png",width=20,height=10)
dev.off()

#making heatmaps: collagen across all cell types
categories = c("Collagens")
selected_genes <- matrisomelist$Gene.Symbol[matrisomelist$Matrisome.Category %in% categories]
selected_genes <- selected_genes[selected_genes %in% rownames(parenchymal@assays$RNA@data)]
DoHeatmap(parenchymal, features=selected_genes, slot="data")
ggsave("/workdir/agh227/kidneydata_figures/paren_heatmap_alltypes_collagens.png",width=20,height=10)
dev.off()

#making heatmaps: proteoglycans across all cell types
categories = c("Proteoglycans")
selected_genes <- matrisomelist$Gene.Symbol[matrisomelist$Matrisome.Category %in% categories]
selected_genes <- selected_genes[selected_genes %in% rownames(parenchymal@assays$RNA@data)]
DoHeatmap(parenchymal, features=selected_genes, slot="data")
ggsave("/workdir/agh227/kidneydata_figures/paren_heatmap_alltypes_proteoglycans.png",width=20,height=10)
dev.off()

#making heatmaps: glycoproteins across all cell types
categories = c("ECM Glycoproteins")
selected_genes <- matrisomelist$Gene.Symbol[matrisomelist$Matrisome.Category %in% categories]
selected_genes <- selected_genes[selected_genes %in% rownames(parenchymal@assays$RNA@data)]
DoHeatmap(parenchymal, features=selected_genes, slot="data")
ggsave("/workdir/agh227/kidneydata_figures/paren_heatmap_alltypes_glycoproteins.png",width=20,height=10)
dev.off()

#focusing on fibroblasts:
selected_celltypes = c("Fibroblasts-1","Fibroblasts-2","Fibroblasts-3")
selected_celltypes <- selected_celltypes[selected_celltypes %in% parenchymal@meta.data$celltype]

subset <- subset(parenchymal, subset = celltype %in% selected_celltypes)

#*making heatmaps, fibroblasts: all core matrisome
categories = c("Collagens","Proteoglycans","ECM Glycoproteins")
selected_genes <- matrisomelist$Gene.Symbol[matrisomelist$Matrisome.Category %in% categories]
selected_genes <- selected_genes[selected_genes %in% rownames(parenchymal@assays$RNA@data)]
DoHeatmap(subset,features=selected_genes, slot="data")
ggsave("/workdir/agh227/kidneydata_figures/paren_heatmap_fibroblasts_corematrisome.png",width=10,height=40)
dev.off()

#*making heatmaps, other cell types with collagen expression: all core matrisomes

selected_celltypes = c("Podocyte","Mesangial","EGM")
selected_celltypes <- selected_celltypes[selected_celltypes %in% parenchymal@meta.data$celltype]

subset <- subset(parenchymal, subset = celltype %in% selected_celltypes)

categories = c("Collagens","Proteoglycans","ECM Glycoproteins")
selected_genes <- matrisomelist$Gene.Symbol[matrisomelist$Matrisome.Category %in% categories]
selected_genes <- selected_genes[selected_genes %in% rownames(parenchymal@assays$RNA@data)]
DoHeatmap(subset,features=selected_genes, slot="data")

ggsave("/workdir/agh227/kidneydata_figures/paren_heatmap_3othersofinterest_corematrisome.png",width=10,height=40)
dev.off()