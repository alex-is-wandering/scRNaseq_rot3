library(Seurat)
library(dplyr)
library(patchwork)
library(ComplexHeatmap)
library(ggplot2)

#1. "find all" "pa/renchymal" and replace with "im/mune"
#2. "find all" "pa/ren" and replace with "im/m"

###### READING IN DATA ###### 
setwd("/workdir/data/crosstissue_data/kidney")
immune<-readRDS("2021-01-22.human_baseline_immune_v3.RDS")
matrisomelist <- read.delim("/local/workdir/agh227/reference_ECM_signatures/Hs_Matrisome_Masterlist_Naba et al_2012.xlsx - Hs_Matrisome_Masterlist.tsv")

#head(immune@meta.data,5) 
#set ident to celltype, for granular vs high-level
Idents(object = immune) <- immune$celltype
DimPlot(immune,reduction="umap")
ggsave("/workdir/agh227/kidneydata_figures/immume_UMAP.png",width=10,height=10)
dev.off()

#cell origin: region of the kidney
imm_cortex_subset <- subset(immune, subset=location=="C")
imm_medulla_subset <- subset(immune, subset=location=="M")
imm_hilum_subset <- subset(immune, subset=location=="P")

###### ECM EXPLORATION ###### 

imm_geneIDs <- rownames(immune@assays$RNA@data) #not scale.data
COL_genes <- imm_geneIDs[grepl("COL",imm_geneIDs,ignore.case=FALSE)]

###### dotplot & vlnplot: looking at celltypes vs individual collagen subunits ###### 
#VIOLIN PLOT PER-GENE
#VlnPlot(immune,features="COL10A1",slot="data")
VlnPlot(immune,features="COL1A1",slot="data")
ggsave("/workdir/agh227/kidneydata_figures/imm_COL1A1_violin.png",width=10,height=5)
dev.off()

#DOT PLOTS FOR CATEGORIES, ILLUSTRATING PER-GENE
categories = c("Collagens")
selected_genes <- matrisomelist$Gene.Symbol[matrisomelist$Matrisome.Category %in% categories]
#*we only want to map those genes that are named in the matrisomelist
selected_genes <- selected_genes[selected_genes %in% rownames(immune@assays$RNA@data)]

DotPlot(immune,features=selected_genes) + theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave("/workdir/agh227/kidneydata_figures/imm_dotplot_alltypes_collagens.png",width=25,height=10) #not sure why background is wonky, but it works
dev.off()

###### Seurat AddModuleScore: looking at celltypes vs ECM component categories ######

# Loop through each unique category in the gene set - looking @ all the data, not just scale.data
for (Matrisome.Category in unique(matrisomelist$Matrisome.Category)) {
  
  # Remove spaces from the category name
  cleaned_category <- gsub(" ", "_", gsub("-", "_", Matrisome.Category))
  
  selected_genes <- matrisomelist$Gene.Symbol[matrisomelist$Matrisome.Category == Matrisome.Category]
  
  immune <- AddModuleScore(immune, features = list(selected_genes), name = cleaned_category) #uses "data" unless "scale=TRUE"
  
}

saveRDS(immune, file="/workdir/agh227/processeddata/immchymal_matrisomeScored.RDS")
head(immune@meta.data,5)

# next:

cat_counts <- table(matrisomelist$Matrisome.Category)

head(immune@meta.data,5)

unique_categories <- unique(matrisomelist$Matrisome.Category)  

VlnPlot(immune, features="ECM_Glycoproteins1")
ggsave("/workdir/agh227/kidneydata_figures/imm_ViolinGlycoproteins_AddModuleScored.png",width=20,height=10)
dev.off()

VlnPlot(immune, features="Proteoglycans1")
ggsave("/workdir/agh227/kidneydata_figures/imm_ViolinProteoglycans_AddModuleScored.png",width=20,height=10)
dev.off()

VlnPlot(immune, features="Collagens1")
ggsave("/workdir/agh227/kidneydata_figures/imm_ViolinCollagens_AddModuleScored.png",width=20,height=10)
dev.off()

VlnPlot(immune, features="ECM_affiliated_Proteins1")
ggsave("/workdir/agh227/kidneydata_figures/imm_ViolinECM_affiliated_AddModuleScored.png",width=20,height=10)
dev.off()

VlnPlot(immune, features="Secreted_Factors1")
ggsave("/workdir/agh227/kidneydata_figures/imm_ViolinSecretedFactors_AddModuleScored.png",width=20,height=10)
dev.off()

VlnPlot(immune, features="ECM_Regulators1")
ggsave("/workdir/agh227/kidneydata_figures/imm_ViolinECMRegulators_AddModuleScored.png",width=20,height=10)
dev.off()

###### Seurat DoHeatmap: gene expression vs ECM categories ###### 
#GetAssayData(immune,slot='scale.data')[1:10,1:5]
#expression data: immune@assays$RNA@scale.data

unique_celltypes <- unique(immune$celltype)
unique_categories <- unique(matrisomelist$Matrisome.Category)

#*naming which categories to map
categories = c("Collagens","Proteoglycans","ECM Glycoproteins")
#matrisome-affiliated: "ECM-affiliated Proteins","ECM Regulators","Secreted Factors"
#core matrisome: "Collagens","Proteoglycans","ECM Glycoproteins"
#*specifying geneID location
selected_genes <- matrisomelist$Gene.Symbol[matrisomelist$Matrisome.Category %in% categories]
#*we only want to map those genes that are named in the matrisomelist
selected_genes <- selected_genes[selected_genes %in% rownames(immune@assays$RNA@data)]

#*making heatmaps: ALL CORE MATRISOME COMPONENTS (3) categories, all cell types
DoHeatmap(immune, features=selected_genes, slot="data")
ggsave("/workdir/agh227/kidneydata_figures/imm_heatmap_alltypes_corematrisome.png",width=20,height=10)
dev.off()

#making heatmaps: collagen across all cell types
categories = c("Collagens")
selected_genes <- matrisomelist$Gene.Symbol[matrisomelist$Matrisome.Category %in% categories]
selected_genes <- selected_genes[selected_genes %in% rownames(immune@assays$RNA@data)]
DoHeatmap(immune, features=selected_genes, slot="data")
ggsave("/workdir/agh227/kidneydata_figures/imm_heatmap_alltypes_collagens.png",width=20,height=10)
dev.off()

#making heatmaps: proteoglycans across all cell types
categories = c("Proteoglycans")
selected_genes <- matrisomelist$Gene.Symbol[matrisomelist$Matrisome.Category %in% categories]
selected_genes <- selected_genes[selected_genes %in% rownames(immune@assays$RNA@data)]
DoHeatmap(immune, features=selected_genes, slot="data")
ggsave("/workdir/agh227/kidneydata_figures/imm_heatmap_alltypes_proteoglycans.png",width=20,height=10)
dev.off()

#making heatmaps: glycoproteins across all cell types
categories = c("ECM Glycoproteins")
selected_genes <- matrisomelist$Gene.Symbol[matrisomelist$Matrisome.Category %in% categories]
selected_genes <- selected_genes[selected_genes %in% rownames(immune@assays$RNA@data)]
DoHeatmap(immune, features=selected_genes, slot="data")
ggsave("/workdir/agh227/kidneydata_figures/imm_heatmap_alltypes_glycoproteins.png",width=20,height=10)
dev.off()


##### points of interest: ECM regulators and ECM affiliated #####

#focusing on some of these monocytes & macrophages:
selected_celltypes = c("LYVE1+ResMacs-1","LYVE1+ResMacs-2","CD9+ TREM2+","mito-rich macrophage","T cells-1")
selected_celltypes <- selected_celltypes[selected_celltypes %in% immune@meta.data$celltype]

subset <- subset(immune, subset = celltype %in% selected_celltypes)

#*making heatmaps, specifically for ECM regulators and affiliated proteins
categories = c("ECM Regulators","ECM-affiliated Proteins")
selected_genes <- matrisomelist$Gene.Symbol[matrisomelist$Matrisome.Category %in% categories]
selected_genes <- selected_genes[selected_genes %in% rownames(immune@assays$RNA@data)]
DoHeatmap(subset,features=selected_genes, slot="data")
ggsave("/workdir/agh227/kidneydata_figures/imm_heatmap_macros_regulatorsandaffiliate.png",width=15,height=45)
dev.off()

#*making dotplot, covering interesting regulatory expression

categories = c("ECM Regulators")
selected_genes <- matrisomelist$Gene.Symbol[matrisomelist$Matrisome.Category %in% categories]
selected_genes <- selected_genes[selected_genes %in% rownames(immune@assays$RNA@data)]

DotPlot(immune,features=selected_genes) + theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave("/workdir/agh227/kidneydata_figures/imm_dotplot_alltypes_ECMregulators.png",width=25,height=10)
dev.off()