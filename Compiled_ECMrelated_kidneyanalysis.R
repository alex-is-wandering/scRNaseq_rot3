setwd("/local1/workdir/data/crosstissue_data/human/kidney/subramanian2021")
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(ComplexHeatmap)
library(Matrix)

#DATA
parenchymal<-readRDS("/local1/workdir/data/crosstissue_data/human/kidney/subramanian2021/2020-06-07.human_baseline_paren_v3.RDS")
immune<-readRDS("/local1/workdir/data/crosstissue_data/human/kidney/subramanian2021/2021-01-22.human_baseline_immune_v3.RDS")
endothelial<-readRDS("/local1/workdir/data/crosstissue_data/human/kidney/subramanian2021/2021-01-22.human_baseline_endo_v3.RDS")#important:
parenchymal <- ScaleData(parenchymal, features=rownames(parenchymal))
immune <- ScaleData(immune, features=rownames(immune))
endothelial <- ScaleData(endothelial, features=rownames(endothelial))

Idents(object = parenchymal) <- parenchymal$celltype
Idents(object = immune) <- immune$celltype
Idents(object = endothelial) <- endothelial$celltype

Idents(endothelial) <- "all"
Idents(immune) <- "all"
Idents(parenchymal) <- "all"

QCplots_endo <- VlnPlot(endothelial, features=c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol=3)
QCplots_imm <- VlnPlot(immune, features=c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol=3)
QCplots_paren <- VlnPlot(parenchymal, features=c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol=3)
QCplots_endo
QCplots_imm
QCplots_paren

#SIGNATURES
matrisomelist <- read.delim("/local1/workdir/agh227/reference_ECM_signatures/Hs_Matrisome_Masterlist_Naba et al_2012.xlsx - Hs_Matrisome_Masterlist.tsv")
reactome_ECMOrg <- read.delim("/workdir/agh227/reference_ECM_signatures/reactome_ECMOrganization.tsv")
#genes at "reactome_ECMOrg$associations.gene.symbol"
KEGG_ECMrecint <- read.delim("/local1/workdir/agh227/reference_ECM_signatures/KEGG_ECMReceptorInteractions.txt")

#clarifying how it's grouping -> each object needs additional metadata. 
#*NMF = celltype?
head(parenchymal@meta.data,5)
colnames(parenchymal@meta.data) #needs condition, highlevel(?)
colnames(immune@meta.data) # needs NMF
unique(immune@meta.data$res.1)
colnames(endothelial@meta.data) #join by orig.ident (need age, sex, condition, location, nephrectomy, highlevel?)


##### module scoring: Matrisome Project ######

#all individual categories (6) vs parenchymal data - the loop that scores by association with the 6 known matrisome categories and saves them 

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


##### module scoring: Reactome, ECM Org template ######
selected_genes <- reactome_ECMOrg$associations.gene.symbol
selected_genes <- selected_genes[selected_genes %in% rownames(parenchymal@assays$RNA@data)]
#DotPlot(parenchymal,features=selected_genes) + theme(axis.text.x=element_text(angle=45, hjust=1))
#ggsave("/workdir/agh227/kidneydata_figures/parenchymal/Reactome_ECMOrg/paren_dotplot_alltypes_collagens.png",width=25,height=10)
#dev.off()

parenchymal <- AddModuleScore(parenchymal, features = list(selected_genes), name = "ECMOrg")
VlnPlot(parenchymal, features="ECMOrg1")
ggsave("/workdir/agh227/kidneydata_figures/parenchymal/Reactome_ECMOrg/paren_ViolinECMOrg_AddModuleScored.png",width=20,height=10)
dev.off()

##### module scoring: KEGG, ECM receptor interaction ######
#immune
selected_genes <- KEGG_ECMrecint$geneID
selected_genes <- selected_genes[selected_genes %in% rownames(immune@assays$RNA@data)]

Idents(object = immune) <- immune$celltype

immune <- AddModuleScore(immune, features = list(selected_genes), name = "ECMRecIntNarrow")
VlnPlot(immune, features="ECMRecIntNarrow1")
ggsave("/workdir/agh227/subramanian2021_kidney/KEGGPathway_Comparisons/imm_ViolinECMReceptorIntNarrow_AddModuleScored.png",width=10,height=5)
dev.off()

#parenchymal
selected_genes <- KEGG_ECMrecint$geneID
selected_genes <- selected_genes[selected_genes %in% rownames(parenchymal@assays$RNA@data)]

parenchymal <- AddModuleScore(parenchymal, features = list(selected_genes), name = "ECMRecIntNarrow")
VlnPlot(parenchymal, features="ECMRecIntNarrow1")
ggsave("/workdir/agh227/subramanian2021_kidney/KEGGPathway_Comparisons/paren_ViolinECMReceptorIntNarrow_AddModuleScored.png",width=10,height=5)
dev.off()

#endothelial
selected_genes <- KEGG_ECMrecint$geneID
selected_genes <- selected_genes[selected_genes %in% rownames(endothelial@assays$RNA@data)]

endothelial <- AddModuleScore(endothelial, features = list(selected_genes), name = "ECMRecIntNarrow")
VlnPlot(endothelial, features="ECMRecIntNarrow1")
ggsave("/workdir/agh227/subramanian2021_kidney/KEGGPathway_Comparisons/endo_ViolinECMReceptorIntNarrow_AddModuleScored.png",width=10,height=5)
dev.off()

##### ComplexHeatmap for genes of interest #####

unique_celltypes <- unique(immune$celltype)
unique_categories <- unique(matrisomelist$Matrisome.Category)
unique(matrisomelist$Matrisome.Division)

###### parenchymal - core matrisome all #####

categories = c("ECM Glycoproteins", "Proteoglycans", "Collagens")

cluster_anno <- parenchymal@meta.data$celltype
col_fun <- circlize::colorRamp2(c(-1,0,3),c("#FF00FF", "black", "#FFFF00"))

selected_genes <- matrisomelist$Gene.Symbol[matrisomelist$Matrisome.Category %in% categories]
selected_genes <- selected_genes[selected_genes %in% rownames(parenchymal@assays$RNA@scale.data)]

mat <- parenchymal[["RNA"]]@scale.data[selected_genes, ] %>% as.matrix() #moving from sparse to full matrix

#filtering heatmap entries by presence of expression:
nonzero_counts <- rowSums(mat !=0)
nonzero_proportion <- nonzero_counts/ncol(mat)
row_stdev <- apply(mat,1,sd)
row_median <- apply(mat, 1, median)

filtered_genes <- selected_genes[nonzero_proportion >= 0.10 & row_stdev > 0.5 & (row_median > 0.25 | row_median < -0.25)]
print(length(filtered_genes))

mat <- mat[filtered_genes, ]

ht_opt$TITLE_PADDING = unit(c(5,15),"points")
Paren_CoreMatrisome <- Heatmap(mat, 
                     name="Core Components", 
                     cluster_columns = T, 
                     show_column_names = F, 
                     show_column_dend = F,
                     cluster_column_slices = F,
                     column_split=factor(cluster_anno), 
                     column_title_gp = gpar(fontsize=8),
                     column_gap = unit(0.2, "mm"), 
                     cluster_rows = T,
                     row_names_gp = gpar(fontsize = 6),
                     show_row_dend = F,
                     col = col_fun,
                     use_raster = T, 
                     raster_quality = 6,
                     column_title_rot=90)
Paren_CoreMatrisome 

###### immune v Matrisome - collagens #####

categories = c("Collagens")

cluster_anno <- immune@meta.data$celltype
col_fun <- circlize::colorRamp2(c(-1,0,3),c("#FF00FF", "black", "#FFFF00"))

selected_genes <- matrisomelist$Gene.Symbol[matrisomelist$Matrisome.Category %in% categories]
selected_genes <- selected_genes[selected_genes %in% rownames(immune@assays$RNA@scale.data)]

mat <- immune[["RNA"]]@scale.data[selected_genes, ] %>% as.matrix()

nonzero_counts <- rowSums(mat !=0)
nonzero_proportion <- nonzero_counts/ncol(mat)
row_stdev <- apply(mat,1,sd)
row_median <- apply(mat, 1, median)

#filtered_genes <- selected_genes[nonzero_proportion >= 0.10 & (abs(row_median) > 0.1) & row_stdev >0.1]
#print(length(filtered_genes))

mat <- mat[selected_genes, ]

ht_opt$TITLE_PADDING = unit(c(5,35),"points")
Imm_Collagens <- Heatmap(mat, 
                               name="Collagens", 
                               cluster_columns = T, 
                               show_column_names = F, 
                               show_column_dend = F,
                               cluster_column_slices = F,
                               column_split=factor(cluster_anno), 
                               column_title_gp = gpar(fontsize=8),
                               column_gap = unit(0.2, "mm"), 
                               cluster_rows = T,
                               row_names_gp = gpar(fontsize = 6),
                               show_row_dend = F,
                               col = col_fun,
                               use_raster = T, 
                               raster_quality = 6,
                               column_title_rot=90)
Imm_Collagens 

###### Endothelial v Matrisome - Proteoglycans #####
unique(matrisomelist$Matrisome.Category)
categories = c("Proteoglycans")

cluster_anno <- endothelial@meta.data$celltype
col_fun <- circlize::colorRamp2(c(-1,0,3),c("#FF00FF", "black", "#FFFF00"))

selected_genes <- matrisomelist$Gene.Symbol[matrisomelist$Matrisome.Category %in% categories]
selected_genes <- selected_genes[selected_genes %in% rownames(endothelial@assays$RNA@scale.data)]

mat <- endothelial[["RNA"]]@scale.data[selected_genes, ] %>% as.matrix()

nonzero_counts <- rowSums(mat !=0)
nonzero_proportion <- nonzero_counts/ncol(mat)
row_stdev <- apply(mat,1,sd)
row_median <- apply(mat, 1, median)

filtered_genes <- selected_genes[nonzero_proportion >= 0.10 & (abs(row_median) > 0.1) & row_stdev >0.1]
print(length(filtered_genes))

mat <- mat[filtered_genes, ]

ht_opt$TITLE_PADDING = unit(c(5,40),"points")
Endo_Proteoglycans <- Heatmap(mat, 
                           name="Core Proteoglycans", 
                           cluster_columns = T, 
                           show_column_names = F, 
                           show_column_dend = F,
                           cluster_column_slices = F,
                           column_split=factor(cluster_anno), 
                           column_title_gp = gpar(fontsize=6),
                           column_gap = unit(0.2, "mm"), 
                           cluster_rows = T,
                           row_names_gp = gpar(fontsize = 8),
                           show_row_dend = F,
                           col = col_fun,
                           use_raster = T, 
                           raster_quality = 6,
                           column_title_rot=90)
Endo_Proteoglycans 

###### Endothelial v Matrisome - Glycoproteins #####

categories = c("ECM Glycoproteins")

cluster_anno <- endothelial@meta.data$celltype
col_fun <- circlize::colorRamp2(c(-1,0,3),c("#FF00FF", "black", "#FFFF00"))

selected_genes <- matrisomelist$Gene.Symbol[matrisomelist$Matrisome.Category %in% categories]
selected_genes <- selected_genes[selected_genes %in% rownames(endothelial@assays$RNA@scale.data)]

mat <- endothelial[["RNA"]]@scale.data[selected_genes, ] %>% as.matrix()

nonzero_counts <- rowSums(mat !=0)
nonzero_proportion <- nonzero_counts/ncol(mat)
row_stdev <- apply(mat,1,sd)
row_median <- apply(mat, 1, median)

filtered_genes <- selected_genes[nonzero_proportion >= 0.10 & (abs(row_median) > 0.1) & row_stdev >0.1]
print(length(filtered_genes))

mat <- mat[filtered_genes, ]

ht_opt$TITLE_PADDING = unit(c(5,45),"points")
Endo_Glycoproteins <- Heatmap(mat, 
                               name="Core Glycoproteins", 
                               cluster_columns = T, 
                               show_column_names = F, 
                               show_column_dend = F,
                               cluster_column_slices = F,
                               column_split=factor(cluster_anno), 
                               column_title_gp = gpar(fontsize=8),
                               column_gap = unit(0.2, "mm"), 
                               cluster_rows = T,
                               row_names_gp = gpar(fontsize = 6),
                               show_row_dend = F,
                               col = col_fun,
                               use_raster = T, 
                               raster_quality = 4,
                               column_title_rot=90)
Endo_Glycoproteins


###### Endothelial v matrisome - ECM regulators #####

categories = c("ECM Regulators")
cluster_anno <- endothelial@meta.data$celltype
col_fun <- circlize::colorRamp2(c(-1,0,3),c("#FF00FF", "black", "#FFFF00"))

selected_genes <- matrisomelist$Gene.Symbol[matrisomelist$Matrisome.Category %in% categories]
selected_genes <- selected_genes[selected_genes %in% rownames(endothelial@assays$RNA@scale.data)]

  mat <- endothelial[["RNA"]]@scale.data[selected_genes, ] %>% as.matrix()

  nonzero_counts <- rowSums(mat !=0)
  nonzero_proportion <- nonzero_counts/ncol(mat)
  row_stdev <- apply(mat,1,sd)
  row_median <- apply(mat, 1, median)

  filtered_genes <- selected_genes[nonzero_proportion >= 0.10 & (abs(row_median) > 0.1) & row_stdev >0.1]
  print(length(filtered_genes))

mat <- mat[filtered_genes, ]

ht_opt$TITLE_PADDING = unit(c(5,35),"points")
Endo_Regulators <- Heatmap(mat, 
                           name="Expression: ECM Regulators", 
                           cluster_columns = T, 
                           show_column_names = F, 
                           show_column_dend = F,
                           cluster_column_slices = F,
                           column_split=factor(cluster_anno), 
                           column_title_gp = gpar(fontsize=8),
                           column_gap = unit(0.2, "mm"), 
                           cluster_rows = T,
                           row_names_gp = gpar(fontsize = 5),
                           show_row_dend = F,
                           col = col_fun,
                           use_raster = T, 
                           raster_quality = 6,
                           column_title_rot=90)
Endo_Regulators 

###### parenchymal v matrisome - ECM regulators #####

categories = c("ECM Regulators")
cluster_anno <- parenchymal@meta.data$celltype
col_fun <- circlize::colorRamp2(c(-1,0,3),c("#FF00FF", "black", "#FFFF00"))

selected_genes <- matrisomelist$Gene.Symbol[matrisomelist$Matrisome.Category %in% categories]
selected_genes <- selected_genes[selected_genes %in% rownames(parenchymal@assays$RNA@scale.data)]

mat <- parenchymal[["RNA"]]@scale.data[selected_genes, ] %>% as.matrix()

nonzero_counts <- rowSums(mat !=0)
nonzero_proportion <- nonzero_counts/ncol(mat)
row_stdev <- apply(mat,1,sd)
row_median <- apply(mat, 1, median)

filtered_genes <- selected_genes[nonzero_proportion >= 0.30 & (abs(row_median) > 0.2) & row_stdev >0.2]
print(length(filtered_genes))

mat <- mat[filtered_genes, ]

ht_opt$TITLE_PADDING = unit(c(5,25),"points")
Paren_Regulators <- Heatmap(mat, 
                           name="ECM Reg.", 
                           cluster_columns = T, 
                           show_column_names = F, 
                           show_column_dend = F,
                           cluster_column_slices = F,
                           column_split=factor(cluster_anno), 
                           column_title_gp = gpar(fontsize=8),
                           column_gap = unit(0.2, "mm"), 
                           cluster_rows = T,
                           row_names_gp = gpar(fontsize = 4),
                           show_row_dend = F,
                           col = col_fun,
                           use_raster = T, 
                           raster_quality = 6,
                           column_title_rot=90)
Paren_Regulators 

###### Endothelial v Matrisome - Secreted Factors #####
#unique(matrisomelist$Matrisome.Category)
categories = c("Secreted Factors")

cluster_anno <- endothelial@meta.data$celltype
col_fun <- circlize::colorRamp2(c(-1,0,3),c("#FF00FF", "black", "#FFFF00"))

selected_genes <- matrisomelist$Gene.Symbol[matrisomelist$Matrisome.Category %in% categories]
selected_genes <- selected_genes[selected_genes %in% rownames(endothelial@assays$RNA@scale.data)]

mat <- endothelial[["RNA"]]@scale.data[selected_genes, ] %>% as.matrix()

  nonzero_counts <- rowSums(mat !=0)
  nonzero_proportion <- nonzero_counts/ncol(mat)
  row_stdev <- apply(mat,1,sd)
  row_median <- apply(mat, 1, median)

  filtered_genes <- selected_genes[nonzero_proportion >= 0.10 & (abs(row_median) > 0.2) & row_stdev >0.9]
  print(length(filtered_genes))

mat <- mat[filtered_genes, ]

ht_opt$TITLE_PADDING = unit(c(5,35),"points")
Endo_SecretedFactors <- Heatmap(mat, 
                               name="Secreted Factors", 
                               cluster_columns = T, 
                               show_column_names = F, 
                               show_column_dend = F,
                               cluster_column_slices = F,
                               column_split=factor(cluster_anno), 
                               column_title_gp = gpar(fontsize=8),
                               column_gap = unit(0.2, "mm"), 
                               cluster_rows = T,
                               row_names_gp = gpar(fontsize = 8),
                               show_row_dend = F,
                               col = col_fun,
                               use_raster = T, 
                               raster_quality = 2,
                               column_title_rot=90)
Endo_SecretedFactors 

###### endothelial v Matrisome - Affiliated Proteins #####
#unique(matrisomelist$Matrisome.Category)
categories = c("ECM-affiliated Proteins")

cluster_anno <- endothelial@meta.data$celltype
col_fun <- circlize::colorRamp2(c(-1,0,3),c("#FF00FF", "black", "#FFFF00"))

selected_genes <- matrisomelist$Gene.Symbol[matrisomelist$Matrisome.Category %in% categories]
selected_genes <- selected_genes[selected_genes %in% rownames(endothelial@assays$RNA@scale.data)]

mat <- endothelial[["RNA"]]@scale.data[selected_genes, ] %>% as.matrix()

nonzero_counts <- rowSums(mat !=0)
nonzero_proportion <- nonzero_counts/ncol(mat)
row_stdev <- apply(mat,1,sd)
row_median <- apply(mat, 1, median)

filtered_genes <- selected_genes[nonzero_proportion >= 0.10 & (abs(row_median) > 0.1) & row_stdev >0.1]
print(length(filtered_genes))

mat <- mat[filtered_genes, ]

ht_opt$TITLE_PADDING = unit(c(5,30),"points")
Endo_ECMAffil <- Heatmap(mat, 
                                name="ECMAffiliated", 
                                cluster_columns = T, 
                                show_column_names = F, 
                                show_column_dend = F,
                                cluster_column_slices = F,
                                column_split=factor(cluster_anno), 
                                column_title_gp = gpar(fontsize=8),
                                column_gap = unit(0.2, "mm"), 
                                cluster_rows = T,
                                row_names_gp = gpar(fontsize = 8),
                                show_row_dend = F,
                                col = col_fun,
                                use_raster = T, 
                                raster_quality = 2,
                                column_title_rot=90)
Endo_ECMAffil 

##### Heatmapping cluster markers: ####
#ECM Regulators different between immune cell clusters:
categories = c("ECM Regulators")
selected_genes <- matrisomelist$Gene.Symbol[matrisomelist$Matrisome.Category %in% categories]
selected_genes <- selected_genes[selected_genes %in% rownames(immune@assays$RNA@scale.data)]
immune_regulators <- FindAllMarkers(immune, features=selected_genes)
head(immune_regulators)

immune_regulators %>% group_by(cluster) %>% slice_head(n=10) %>% ungroup() -> top10
DoHeatmap(immune, features=top10$gene)
ggsave("/local1/workdir/agh227/subramanian2021_kidney/immune_ECMRegulators.png", width=22, height=10)
dev.off()


##### parenchymal: differences between fibroblasts-1, -2, -3 #####
unique(parenchymal@meta.data$celltype)
fibroblasts <- c("Fibroblasts-1", "Fibroblasts-2","Fibroblasts-3")
fibroblasts_kidney <- subset(parenchymal, subset = celltype %in% fibroblasts)
unique(fibroblasts_kidney@meta.data$celltype)

head(fibroblasts_kidney@meta.data)

fibroblast_markers <- FindAllMarkers(fibroblasts_kidney)

top_markers <- fibroblast_markers %>%
  dplyr::filter(abs(avg_log2FC) >0.5) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n=30) %>%
  ungroup()
top_genes <- top_markers$gene

DoHeatmap(fibroblasts_kidney, features=top_genes)
ggsave("/local1/workdir/agh227/subramanian2021_kidney/heatmap_fibroblastmarkers.png", width = 8, height = 9)
dev.off()

##### immune: T-cells vs macrophages #####

unique(immune@meta.data$celltype)
immune_TandMac <- c("T cells-1", "T cells-2 (stress)", "Activated T cells", "mito-rich macrophage","LYVE1+ResMacs-2", "LYVE1+ResMacs-1", "CD9+ TREM2+")
immune_TMacro_comparison <- subset(immune, subset = celltype %in% immune_TandMac)
unique(immune_TMacro_comparison@meta.data$celltype)
Idents(object = immune_TMacro_comparison) <- immune_TMacro_comparison$celltype

head(immune_TMacro_comparison@meta.data)

TMacro_markers <- FindAllMarkers(immune_TMacro_comparison)

top_markers <- TMacro_markers %>%
  dplyr::filter(abs(avg_log2FC) >0.5) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n=30) %>%
  ungroup()
top_genes <- top_markers$gene

DoHeatmap(immune_TMacro_comparison, features=top_genes)
ggsave("/local1/workdir/agh227/subramanian2021_kidney/heatmap_Tcellmacrophagecomparisons.png", width = 8, height = 9)
dev.off()


##### MERGING ALL KIDNEY DATA ####
#merge for bulk check of gene expression
#by default, this erases previously normalize and scaled data matrices. save the normalized ones with "merge.data=TRUE"
merged_endoparen <- merge (x=endothelial, y=parenchymal, merge.data=T)
merged_kidney <- merge(x=merged_endoparen, y=immune, merge.data=T)
saveRDS(object=merged_kidney, file="/local1/workdir/agh227/subramanian2021_kidney/mergedkidney.RDS")

summary(merged_kidney@assays$RNA@data[1:20, 1:15])
#(merged_kidney@assays$RNA@data[1:10, 1:15]) #doesn't work
colnames(merged_kidney@meta.data)
unique(merged_kidney@meta.data$celltype)
Idents(object = merged_kidney) <- merged_kidney$celltype
merged_kidney <- ScaleData(merged_kidney, features=rownames(merged_kidney))

####*ALL KIDNEY V. COLLAGENS: HEATMAP ####
inquiry <- c(unique(merged_kidney$celltype))
verifkidney_check <- subset(merged_kidney, subset = celltype %in% inquiry)
table(verifkidney_check@meta.data$celltype)

categories = c("ECM Glycoproteins")
#unique(matrisomelist$Matrisome.Category)
cluster_anno <- verifkidney_check@meta.data$celltype
col_fun <- circlize::colorRamp2(c(-1,0,3),c("#FF00FF", "black", "#FFFF00"))

selected_genes <- matrisomelist$Gene.Symbol[matrisomelist$Matrisome.Category %in% categories]
selected_genes <- selected_genes[selected_genes %in% rownames(verifkidney_check@assays$RNA@scale.data)]

mat <- verifkidney_check[["RNA"]]@scale.data[selected_genes, ] %>% as.matrix()

#filtering heatmap entries by presence of expression:
nonzero_counts <- rowSums(mat !=0)
nonzero_proportion <- nonzero_counts/ncol(mat)
row_stdev <- apply(mat,1,sd)
row_median <- apply(mat, 1, median)

filtered_genes <- selected_genes[nonzero_proportion >= 0.1 & row_stdev > 0.1 & (row_median > 0.1 | row_median < -0.1)]
print(length(filtered_genes))

mat <- mat[selected_genes, ]

ht_opt$TITLE_PADDING = unit(c(5,15),"points")
subrakidney_Glycoprot <- Heatmap(mat, 
                                 name="All v. Glycoproteins", 
                                 cluster_columns = T, 
                                 show_column_names = F, 
                                 show_column_dend = F,
                                 cluster_column_slices = F,
                                 column_split=factor(cluster_anno), 
                                 column_title_gp = gpar(fontsize=8),
                                 column_gap = unit(0.6, "mm"), 
                                 cluster_rows = T,
                                 row_names_gp = gpar(fontsize = 8),
                                 show_row_dend = F,
                                 col = col_fun,
                                 use_raster = T, 
                                 raster_quality = 2,
                                 column_title_rot=90)
subrakidney_Glycoprot 
#ggsave("local1/workdir/agh227/kidney_MatrisomeCollagens.png")
dev.off()


##### collagens: dotplot across ALL cells ####
categories = c("Collagens")
selected_genes <- matrisomelist$Gene.Symbol[matrisomelist$Matrisome.Category %in% categories]
selected_genes <- selected_genes[selected_genes %in% rownames(merged_kidney@assays$RNA@data)]

kidney_collagens <- DotPlot(merged_kidney,features=selected_genes) + theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave("/workdir/agh227/subramanian2021_kidney/kidneyDotPlot_MatrisomeCollagenExpression.png", width=20, height=30)
dev.off()


##### looking @ MMP expression ####
paren_geneIDs <- rownames(immune@assays$RNA@data)
MMP_genes <- paren_geneIDs[grepl("^MMP",paren_geneIDs,ignore.case=FALSE)]
cluster_anno <- parenchymal@meta.data$celltype

MMP_genes <- MMP_genes[MMP_genes %in% rownames(parenchymal@assays$RNA@scale.data)]

paren_MMPs <- DotPlot(parenchymal,features=MMP_genes) + theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave("/local1/workdir/agh227/subramanian2021_kidney/parenchymalDotPlot_MMPExpression.png", width=20, height=30)
dev.off()


##### WHICH COMPONENTS ARE ENRICHED IN THESE CELLS? ####
#mass heatmapping celltypes in merged_kidney vs matrisome
#NOT NECESSARY: redirect to filtering_scaledata_bycluster.R and only make heatmaps to verify; not individual

setwd("/local1/workdir/agh227/subramanian2021_kidney/MatrisomeProject_Comparisons/massheatmaps_bycelltype") #all files saved here
table(endothelial@meta.data$celltype)
table(matrisomelist$Matrisome.Category)
Idents(object = merged_kidney) <- merged_kidney$celltype

#START
#saves to working directory!!!

ht_opt$TITLE_PADDING = unit(c(5,25),"points")
col_fun <- circlize::colorRamp2(c(-1,0,3),c("#FF00FF", "black", "#FFFF00"))

categories_celltypes <- unique(merged_kidney@meta.data$celltype)

for (celltype in categories_celltypes) {
  
  #Collagens, ECM Glycoproteins, Proteoglycans
  #ECM-affiliated Proteins, ECM Regulators, Secreted Factors
  categories = ("ECM Regulators")
  selected_cells <- which(Idents(merged_kidney) == celltype)
  
  cluster_anno <- merged_kidney@meta.data$celltype[selected_cells]
  
  selected_genes <- matrisomelist$Gene.Symbol[matrisomelist$Matrisome.Category %in% categories]
  selected_genes <- selected_genes[selected_genes %in% rownames(merged_kidney@assays$RNA@scale.data)]
  
  mat <- merged_kidney[["RNA"]]@scale.data[selected_genes, selected_cells] %>% as.matrix()
  
  nonzero_counts <- rowSums(mat !=0)
  nonzero_proportion <- nonzero_counts/ncol(mat)
  row_median <- apply(mat, 1, median)
  row_stdev <- apply(mat,1,sd)
  
  filtered_genes <- selected_genes[nonzero_proportion >= 0.1 & abs(row_median) > 0.2] # & row_stdev < 0.8
  
  if(length(filtered_genes)<=1) {
    print(paste("filtered_genes is empty or has insufficient entries; skipping heatmap for", celltype))
  } else {
    filtered_genes <- filtered_genes[!is.na(filtered_genes)]
    if(length(filtered_genes)==0) {
      print(paste("filtered_genes only contains NA entries; skipping heatmap for", celltype))
    } else {
      mat <- mat[filtered_genes[!is.na(filtered_genes)], ]  
      
      category_name_clean <- gsub("/","_",celltype)
      heatmap_name <- paste("Heatmap_", categories, "_", category_name_clean, ".pdf",sep="")
      pdf(heatmap_name)
      
      heatmap <- Heatmap(mat, 
                         name=celltype, 
                         cluster_columns = T, 
                         show_column_names = F, 
                         show_column_dend = F,
                         cluster_column_slices = F,
                         column_split=factor(cluster_anno), 
                         column_title_gp = gpar(fontsize=8),
                         column_gap = unit(0.6, "mm"), 
                         cluster_rows = T,
                         row_names_gp = gpar(fontsize = 8),
                         show_row_dend = F,
                         col = col_fun,
                         use_raster = T, 
                         raster_quality = 2,
                         column_title_rot=90)
      
      draw(heatmap)
      dev.off()
    }}}

#DONE
