#seurat object -> genes in common with term of interest (i.e. matrisome categories) -> filter by per-celltype expression, whether or not they're expressed in enough cells -> removing visual bias & replicating visual analysis of heatmaps

library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(Matrix)

##### reference list ####
setwd("/local1/workdir/agh227/lakemenon2023_kidneyverificationset")

matrisomelist <- read.delim("/local1/workdir/agh227/reference_ECM_signatures/Hs_Matrisome_Masterlist_Naba et al_2012.xlsx - Hs_Matrisome_Masterlist.tsv")

##### data ####

lakemenonkidney2023_healthy <- readRDS("/local1/workdir/data/crosstissue_data/human/kidney/LakeMenon2023/scRNAseq/lakemenon2023kidney_healthy.RDS")
#lakemenonkidney2023_healthy <- ScaleData(lakemenonkidney2023_healthy, features=rownames(lakemenonkidney2023_healthy))
#saveRDS(object=lakemenonkidney2023_healthy, file="/local1/workdir/data/crosstissue_data/human/kidney/LakeMenon2023/scRNAseq/lakemenon2023kidney_healthy.RDS")
merged_kidney<- readRDS("/local1/workdir/agh227/subramanian2021_kidney/mergedkidney.RDS")
#merged_kidney <- ScaleData(merged_kidney, features=rownames(merged_kidney))
#saveRDS(object=merged_kidney, file="/local1/workdir/data/crosstissue_data/human/kidney/subramanian2021/mergedkidney.RDS")


Idents(lakemenonkidney2023_healthy) <- lakemenonkidney2023_healthy@meta.data$subclass.l2
Idents(merged_kidney) <- merged_kidney@meta.data$celltype


table(lakemenonkidney2023_healthy@meta.data$subclass.l2)
table(merged_kidney@meta.data$celltype)
table(matrisomelist$Matrisome.Category)

##### START EXTRACTING GENES OF INTEREST: lakemenon2023 #####

Idents(lakemenonkidney2023_healthy) <- lakemenonkidney2023_healthy@meta.data$subclass.l2
categories_l2 <- unique(lakemenonkidney2023_healthy@meta.data$subclass.l2)
filtered_genes_per_cluster <- list()
gene_cluster_table <- list()
gene_categories <- list()

for (category_l2 in categories_l2) {
  
  #Collagens, ECM Glycoproteins, Proteoglycans
  #ECM-affiliated Proteins, ECM Regulators, Secreted Factors
  
categories = c("Collagens", "ECM Glycoproteins", "Proteoglycans", "ECM-affiliated Proteins", "ECM Regulators", "Secreted Factors")
selected_cells <- which(Idents(lakemenonkidney2023_healthy) == category_l2)
  
selected_genes <- matrisomelist$Gene.Symbol[matrisomelist$Matrisome.Category %in% categories]
selected_genes <- selected_genes[selected_genes %in% rownames(lakemenonkidney2023_healthy@assays$RNA@scale.data)]
  
mat <- lakemenonkidney2023_healthy[["RNA"]]@scale.data[selected_genes, selected_cells] %>% as.matrix()

#filter by per-cluster (as this iterates) expression
cells_above_threshold_per_gene <- apply(mat, 1, function(row) sum(row >=2))
prop_above_threshold_per_gene <- cells_above_threshold_per_gene / length(selected_cells)
filtered_genes <- selected_genes[prop_above_threshold_per_gene >=0.30] 

  nonzero_counts <- rowSums(mat !=0)
  nonzero_proportion <- nonzero_counts/ncol(mat)
  #row_median <- apply(mat, 1, median)
  
filtered_genes <- filtered_genes[nonzero_proportion >= 0.1]# & abs(row_median) > 0.1]
filtered_genes <- na.omit(filtered_genes)

filtered_genes_per_cluster[[category_l2]] <- filtered_genes

for (gene in filtered_genes) {
  if(!(gene %in% names(gene_cluster_table))) {
    gene_cluster_table[[gene]] <- list()
    gene_categories[[gene]] <- matrisomelist$Matrisome.Category[match(gene, matrisomelist$Gene.Symbol)]
  }
gene_cluster_table[[gene]] <- c(gene_cluster_table[[gene]], category_l2)}}

gene_cluster_df <- data.frame(Category = unlist(gene_categories),
Gene = names(gene_cluster_table),
Clusters=sapply(gene_cluster_table, function(clusters) paste(clusters, collapse=",")),
row.names=NULL)

print(gene_cluster_df)

gene_cluster_df_sorted <- gene_cluster_df[order(gene_cluster_df$Category, gene_cluster_df$Gene), ]
write.table(gene_cluster_df_sorted, "lakemenon2023_cellspecificgeneexpression.txt", sep="\t", row.names=F, quote=F)

##### START EXTRACTING GENES OF INTEREST: subramanian2021 #####
Idents(merged_kidney) <- merged_kidney@meta.data$celltype

categories_l2 <- unique(merged_kidney@meta.data$celltype)
filtered_genes_per_cluster <- list()
gene_cluster_table <- list()
gene_categories <- list()

for (category_l2 in categories_l2) {
  
  #Collagens, ECM Glycoproteins, Proteoglycans
  #ECM-affiliated Proteins, ECM Regulators, Secreted Factors
  
  categories = c("Collagens", "ECM Glycoproteins", "Proteoglycans", "ECM-affiliated Proteins", "ECM Regulators", "Secreted Factors")
  selected_cells <- which(Idents(merged_kidney) == category_l2)
  
  selected_genes <- matrisomelist$Gene.Symbol[matrisomelist$Matrisome.Category %in% categories]
  selected_genes <- selected_genes[selected_genes %in% rownames(merged_kidney@assays$RNA@scale.data)]
  
  mat <- merged_kidney[["RNA"]]@scale.data[selected_genes, selected_cells] %>% as.matrix()
  
  #filter by per-cluster (as this iterates) expression
  cells_above_threshold_per_gene <- apply(mat, 1, function(row) sum(row >=2))
  prop_above_threshold_per_gene <- cells_above_threshold_per_gene / length(selected_cells)
  filtered_genes <- selected_genes[prop_above_threshold_per_gene >=0.30] 
  
  nonzero_counts <- rowSums(mat !=0)
  nonzero_proportion <- nonzero_counts/ncol(mat)
  #row_median <- apply(mat, 1, median)
  
  filtered_genes <- filtered_genes[nonzero_proportion >= 0.1]# & abs(row_median) > 0.1]
  filtered_genes <- na.omit(filtered_genes)
  
  filtered_genes_per_cluster[[category_l2]] <- filtered_genes
  
  for (gene in filtered_genes) {
    if(!(gene %in% names(gene_cluster_table))) {
      gene_cluster_table[[gene]] <- list()
      gene_categories[[gene]] <- matrisomelist$Matrisome.Category[match(gene, matrisomelist$Gene.Symbol)]
    }
    gene_cluster_table[[gene]] <- c(gene_cluster_table[[gene]], category_l2)}}

gene_cluster_df <- data.frame(Category = unlist(gene_categories),
Gene = names(gene_cluster_table),
Clusters=sapply(gene_cluster_table, function(clusters) paste(clusters, collapse=",")),
row.names=NULL)

print(gene_cluster_df)

gene_cluster_df_sorted <- gene_cluster_df[order(gene_cluster_df$Category, gene_cluster_df$Gene), ]
write.table(gene_cluster_df_sorted, "subramanian2021_cellspecificgeneexpression.txt", sep="\t", row.names=F, quote=F)


##### reading from that table to make heatmaps #####

specificityinfo_lm23 <- read.delim("lakemenon2023_cellspecificgeneexpression.txt")
specificityinfo_sb21 <- read.delim("subramanian2021_cellspecificgeneexpression.txt")

categories=c("Collagens")
selected_genes <- gene_cluster_df_sorted$Gene[gene_cluster_df_sorted$Category %in% categories]

fibroblasts <- c("FIB", "MyoF","aFIB")
verifkidney_Fibroblasts <- subset(lakemenonkidney2023_healthy, subset = subclass.l2 %in% fibroblasts)

cluster_anno <- verifkidney_Fibroblasts@meta.data$subclass.l2
mat <- verifkidney_Fibroblasts[["RNA"]]@scale.data[selected_genes, ] %>% as.matrix()

heatmap_name <- paste("FilteredCollagens.pdf")
ht_opt$TITLE_PADDING = unit(c(5,35),"points")
col_fun <- circlize::colorRamp2(c(-1,0,3),c("#FF00FF", "black", "#FFFF00"))
pdf(heatmap_name)

heatmap <- Heatmap(mat, 
                   name="Filtered Collagens", 
                   cluster_columns = T, 
                   show_column_names = F, 
                   show_column_dend = F,
                   cluster_column_slices = F,
                   column_split=factor(cluster_anno), 
                   column_title_gp = gpar(fontsize=8),
                   column_gap = unit(0.3, "mm"), 
                   cluster_rows = T,
                   row_names_gp = gpar(fontsize = 8),
                   show_row_dend = F,
                   col = col_fun,
                   use_raster = T, 
                   raster_quality = 2,
                   column_title_rot=90)

draw(heatmap)
dev.off()


##### heatmapping some genes in common #####
#10 glycoproteins. cell types from the table previous, all "specific"
table(kidneycheck_LM@meta.data$subclass.l2)
table(kidneycheck_AS@meta.data$celltype)
selected_genes <- c("TGFBI","THBS1","THBS2","THSD4","TINAG","TINAGL1","TNC","TNXB","VWA1","VWF")

selectedcelltypes_LM <- c("aFIB","FIB", "MyoF","VSMC/P","dVSMC","EC-LYM", "MAC-M2","ncMON", "cycMNP", "MC","pDC")
kidneycheck_LM <- subset(lakemenonkidney2023_healthy, subset = subclass.l2 %in% selectedcelltypes_LM)
cluster_anno_LM <- kidneycheck_LM@meta.data$subclass.l2
mat_LM <- kidneycheck_LM[["RNA"]]@scale.data[selected_genes, ] %>% as.matrix()

selectedcelltypes_AS <- c("Fibroblasts-1","Fibroblasts-2","Fibroblasts-3", "vSMC","Pericyte-2","Lymphatic EC","mito-rich macrophage","EGM","LYVE1+ResMacs-2","LYVE1+ResMacs-1","CD9+ TREM2+","cDC-1","cDC-2","Monocyte-2 (C)","pDC")
kidneycheck_AS <- subset(merged_kidney, subset = celltype %in% selectedcelltypes_AS)
cluster_anno_AS <- kidneycheck_AS@meta.data$celltype
mat_AS <- kidneycheck_AS[["RNA"]]@scale.data[selected_genes, ] %>% as.matrix()

ht_opt$TITLE_PADDING = unit(c(5,35),"points")
col_fun <- circlize::colorRamp2(c(-1,0,3),c("#FF00FF", "black", "#FFFF00"))

heatmap_name_LM <- paste("SampleGlycoprot_lake23.pdf")

pdf(heatmap_name_LM)
heatmap <- Heatmap(mat_LM, 
                   name="Sample Glycoprot LM", 
                   cluster_columns = T, 
                   show_column_names = F, 
                   show_column_dend = F,
                   cluster_column_slices = F,
                   column_split=factor(cluster_anno_LM), 
                   column_title_gp = gpar(fontsize=8),
                   column_gap = unit(0.3, "mm"), 
                   cluster_rows = T,
                   row_names_gp = gpar(fontsize = 8),
                   show_row_dend = F,
                   col = col_fun,
                   use_raster = T, 
                   raster_quality = 2,
                   column_title_rot=90)
heatmap@row_names_param$labels <- order(heatmap@row_names_param$labels)
draw(heatmap)
dev.off()

heatmap_name_AS <- paste("SampleGlycoprot_submn21.pdf")
pdf(heatmap_name_AS)
heatmap <- Heatmap(mat_AS, 
                   name="Sample Glycoprot AS", 
                   cluster_columns = T, 
                   show_column_names = F, 
                   show_column_dend = F,
                   cluster_column_slices = F,
                   column_split=factor(cluster_anno_AS), 
                   column_title_gp = gpar(fontsize=8),
                   column_gap = unit(0.3, "mm"), 
                   cluster_rows = T,
                   row_names_gp = gpar(fontsize = 8),
                   show_row_dend = F,
                   col = col_fun,
                   use_raster = T, 
                   raster_quality = 2,
                   column_title_rot=90)
heatmap@row_names_param$labels <- factor(heatmap@row_names_param$labels,
                                         levels=c("TGFBI","THBS1","THBS2","THSD4","TINAG","TINAGL1","TNC","TNXB","VWA1","VWF"),
                                         ordered=TRUE)
draw(heatmap)
dev.off()
