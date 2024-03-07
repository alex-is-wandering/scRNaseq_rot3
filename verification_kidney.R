#library(harmony)
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(ComplexHeatmap)
library(Matrix)

#load("/local1/workdir/data/crosstissue_data/human/kidney/LakeMenon2023/scRNAseq/PREMIERE_V1_ForAS.Robj")
#kpmpld -> lakemenonkidney2023
#rm(kpmpld)
#saveRDS(object=lakemenonkidney2023, file="/local1/workdir/data/crosstissue_data/human/kidney/LakeMenon2023/scRNAseq/lakemenon2023kidney.RDS")

lakemenonkidney2023_healthy <- readRDS("/local1/workdir/data/crosstissue_data/human/kidney/LakeMenon2023/scRNAseq/lakemenon2023kidney_healthy.RDS")

##### loading in lake menon 2023 data and ECM related files: ####

matrisomelist <- read.delim("/local1/workdir/agh227/reference_ECM_signatures/Hs_Matrisome_Masterlist_Naba et al_2012.xlsx - Hs_Matrisome_Masterlist.tsv")
reactome_ECMOrg <- read.delim("/workdir/agh227/reference_ECM_signatures/reactome_ECMOrganization.tsv")
KEGG_ECMrecint <- read.delim("/local1/workdir/agh227/reference_ECM_signatures/KEGG_ECMReceptorInteractions.txt")

#lakemenonkidney2023 <- readRDS("/local1/workdir/data/crosstissue_data/human/kidney/LakeMenon2023/scRNAseq/lakemenon2023kidney.RDS")

Idents(lakemenonkidney2023) <- "all" #they're clustering wrong, so just shoving that in there.
QCplots <- VlnPlot(lakemenonkidney2023, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
QCplots

#normalizing didn't change anything about the below commands - therefore already done for me.
#lakemenonkidney2023_healthy <- NormalizeData(lakemenonkidney2023_healthy)

#checking gene expression distributions across cells:
Idents(lakemenonkidney2023_healthy) <- "all"
normcheck <- VlnPlot(lakemenonkidney2023_healthy, features = c("IGFBP7", "COL18A1", "VCAN"))

#feature-wise dispersion:
lakemenonkidney2023_healthy <- FindVariableFeatures(lakemenonkidney2023_healthy, selection.method="vst")
plot1 <- VariableFeaturePlot(lakemenonkidney2023_healthy)
plot2 <- LabelPoints(plot1, points=head(VariableFeatures(lakemenonkidney2023_healthy), 10), repel=T)
plot2

#object info & subsetting out the part we're interested in (healthy data)
#colnames(lakemenonkidney2023@meta.data)
#head(lakemenonkidney2023@meta.data,2)
unique(lakemenonkidney2023_healthy@meta.data$SampleType) #LD = healthy
lakemenonkidney2023_healthy <- subset(lakemenonkidney2023, subset = SampleType == "LD")
#saveRDS(object=lakemenonkidney2023_healthy, file="/local1/workdir/data/crosstissue_data/human/kidney/LakeMenon2023/scRNAseq/lakemenon2023kidney_healthy.RDS")
#rm(lakemenonkidney2023)
#DimPlot(lakemenonkidney2023_healthy, reduction="umap")

#unique(lakemenonkidney2023@meta.data$subclass.l2)
#unique(lakemenonkidney2023@meta.data$subclass.l1)
#unique(lakemenonkidney2023@meta.data$dataSource)
#unique(lakemenonkidney2023@meta.data$tissueType)

Idents(lakemenonkidney2023_healthy) <- lakemenonkidney2023_healthy@meta.data$subclass.l2
colnames(lakemenonkidney2023_healthy@meta.data)
length(unique(lakemenonkidney2023_healthy@meta.data$subclass.l2))

#### running against matrisome project categories: ####
# Loop through each unique category in the gene set
for (Matrisome.Category in unique(matrisomelist$Matrisome.Category)) {
  
  selected_genes <- matrisomelist$Gene.Symbol[matrisomelist$Matrisome.Category == Matrisome.Category]
  
  lakemenonkidney2023_healthy <- AddModuleScore(lakemenonkidney2023_healthy, features = list(selected_genes), name = Matrisome.Category)}

#colnames(lakemenonkidney2023_healthy@meta.data)

kidneyverif_Collagen <- VlnPlot(lakemenonkidney2023_healthy, features="Collagens1")
ggsave("/local1/workdir/agh227/lakemenon2023_kidneyverificationset/kidneyverif_Collagen.png", width=15, height=7)
kidneyverif_Glycopro <- VlnPlot(lakemenonkidney2023_healthy, features="ECM.Glycoproteins1")
ggsave("/local1/workdir/agh227/lakemenon2023_kidneyverificationset/kidneyverif_Glycopro.png", width=15, height=7)
kidneyverif_Proteogl <- VlnPlot(lakemenonkidney2023_healthy, features="Proteoglycans1")
ggsave("/local1/workdir/agh227/lakemenon2023_kidneyverificationset/kidneyverif_Proteogl.png", width=15, height=7)
kidneyverif_Affiliat <- VlnPlot(lakemenonkidney2023_healthy, features="ECM.affiliated.Proteins1")
ggsave("/local1/workdir/agh227/lakemenon2023_kidneyverificationset/kidneyverif_Affiliat.png", width=15, height=7)
kidneyverif_Secreted <- VlnPlot(lakemenonkidney2023_healthy, features="Secreted.Factors1")
ggsave("/local1/workdir/agh227/lakemenon2023_kidneyverificationset/kidneyverif_Secreted.png", width=15, height=7)
kidneyverif_Regulato <- VlnPlot(lakemenonkidney2023_healthy, features="ECM.Regulators1")
ggsave("/local1/workdir/agh227/lakemenon2023_kidneyverificationset/kidneyverif_Regulato.png", width=15, height=7)

#### following up with complexheatmaps ####
lakemenonkidney2023_healthy <- ScaleData(lakemenonkidney2023_healthy, features=rownames(lakemenonkidney2023_healthy))

categories = c("Collagens")
cluster_anno <- c
col_fun <- circlize::colorRamp2(c(-1,0,3),c("#FF00FF", "black", "#FFFF00"))

selected_genes <- matrisomelist$Gene.Symbol[matrisomelist$Matrisome.Category %in% categories]
selected_genes <- selected_genes[selected_genes %in% rownames(lakemenonkidney2023_healthy@assays$RNA@scale.data)]

mat <- lakemenonkidney2023_healthy[["RNA"]]@scale.data[selected_genes, ] %>% as.matrix()

#filtering heatmap entries by presence of expression:
nonzero_counts <- rowSums(mat !=0)
nonzero_proportion <- nonzero_counts/ncol(mat)
row_stdev <- apply(mat,1,sd)
row_median <- apply(mat, 1, median)

filtered_genes <- selected_genes[nonzero_proportion >= 0.1 & row_stdev > 0.1 & (row_median > 0.1 | row_median < -0.1)]
print(length(filtered_genes))

mat <- mat[filtered_genes, ]

#making color blocks per cell type, for annotation:
#category_annotation <- lakemenonkidney2023_healthy@meta.data$subclass.l2
#unique_cat <- unique(category_annotation)
#num_cat <- length(unique_cat)
#color_palette <- colorRampPalette(c("blue", "green", "pink"))(num_cat)
#col_annotation_colors <- list(Category=setNames(color_palette, unique_cat))

verifkidney_Collagens <- Heatmap(mat, 
                               name="Core Collagens", 
                               cluster_columns = T, 
                               show_column_names = F, 
                               show_column_dend = F,
                               cluster_column_slices = F,
                               column_split=factor(cluster_anno), 
                               column_gap = unit(0.2, "mm"), 
                               cluster_rows = T,
                               row_names_gp = gpar(fontsize = 6),
                               show_row_dend = F,
                               col = col_fun,
                               use_raster = T, 
                               raster_quality = 6,
                               column_title_rot=90)
unique(lakemenonkidney2023_healthy@meta.data$subclass.l2)
#ha=HeatmapAnnotation(df = data.frame(lakemenonkidney2023_healthy@assays$RNA@scale.data), col=col_annotation_colors)
#draw(ha, 1:10)
verifkidney_Collagens 
#ggsave("verifkidney_Collagens.png", width=15, height=6)

#### what differentiates the FIB celltypes? ####
fibroblasts <- c("FIB", "MyoF","aFIB")
verifkidney_Fibroblasts <- subset(lakemenonkidney2023_healthy, subset = subclass.l2 %in% fibroblasts)
table(verifkidney_Fibroblasts@meta.data$subclass.l2)

#head(verifkidney_Fibroblasts@meta.data)

fibroblast_markers <- FindAllMarkers(verifkidney_Fibroblasts)

top_markers <- fibroblast_markers %>%
  dplyr::filter(abs(avg_log2FC) >2) %>%
  arrange(desc(abs(avg_log2FC))) %>%
  slice_head(n=50) %>%
  ungroup()
top_genes <- top_markers$gene

DoHeatmap(verifkidney_Fibroblasts, features=top_genes)
ggsave("/local1/workdir/agh227/lakemenon2023_kidneyverificationset/heatmap_Fibroblasts.png", width=9, height=7)
dev.off()


##### DotPlot ALL Matrisome categories vs ALL cell types ####
plots_list <- list()
for (Matrisome.Category in unique(matrisomelist$Matrisome.Category)) {
  
  selected_genes <- matrisomelist$Gene.Symbol[matrisomelist$Matrisome.Category == Matrisome.Category]
  selected_genes <- selected_genes[selected_genes %in% rownames(lakemenonkidney2023_healthy@assays$RNA@scale.data)]
  
plot <-  DotPlot(lakemenonkidney2023_healthy,features=selected_genes) + theme(axis.text.x=element_text(angle=45, hjust=1, size=15), axis.text.y=element_text(size=15))

plots_list[[Matrisome.Category]] <- plot}

plots_list$`ECM Glycoproteins`
ggsave("/local1/workdir/agh227/lakemenon2023_kidneyverificationset/id_granular/DotPlotMatrisomeECMGlyco_alltypes.pdf", width=40, height=20)
dev.off()
plots_list$Collagens
ggsave("/local1/workdir/agh227/lakemenon2023_kidneyverificationset/id_granular/DotPlotMatrisomeCollagens_alltypes.pdf", width=40, height=20)
dev.off()
plots_list$Proteoglycans
ggsave("/local1/workdir/agh227/lakemenon2023_kidneyverificationset/id_granular/DotPlotMatrisomeProteogl_alltypes.pdf", width=20, height=20)
dev.off()
plots_list$`ECM-affiliated Proteins`
ggsave("/local1/workdir/agh227/lakemenon2023_kidneyverificationset/id_granular/DotPlotMatrisomeECMAffil_alltypes.pdf", width=40, height=20)
dev.off()
plots_list$`ECM Regulators`
ggsave("/local1/workdir/agh227/lakemenon2023_kidneyverificationset/id_granular/DotPlotMatrisomeECMReg_alltypes.pdf", width=40, height=20)
dev.off()
plots_list$`Secreted Factors`
ggsave("/local1/workdir/agh227/lakemenon2023_kidneyverificationset/id_granular/DotPlotMatrisomeSecrFac_alltypes.pdf", width=40, height=20)
dev.off()


##### WHICH COMPONENTS ARE ENRICHED IN THESE CELLS? ####
#MASS HEATMAP GENERATION GO

setwd("/local1/workdir/agh227/lakemenon2023_kidneyverificationset/id_granular/mass_heatmapping") #all files saved here
table(lakemenonkidney2023_healthy@meta.data$subclass.l2)
table(matrisomelist$Matrisome.Category)
Idents(lakemenonkidney2023_healthy) <- lakemenonkidney2023_healthy@meta.data$subclass.l2
#START
#saves to working directory!!!

ht_opt$TITLE_PADDING = unit(c(5,35),"points")
col_fun <- circlize::colorRamp2(c(-1,0,3),c("#FF00FF", "black", "#FFFF00"))

categories_l2 <- unique(lakemenonkidney2023_healthy@meta.data$subclass.l2)

for (category_l2 in categories_l2) {

  #Collagens, ECM Glycoproteins, Proteoglycans
  #ECM-affiliated Proteins, ECM Regulators, Secreted Factors
categories = ("Collagens")
selected_cells <- which(Idents(lakemenonkidney2023_healthy) == category_l2)

cluster_anno <- lakemenonkidney2023_healthy@meta.data$subclass.l2[selected_cells]

selected_genes <- matrisomelist$Gene.Symbol[matrisomelist$Matrisome.Category %in% categories]
selected_genes <- selected_genes[selected_genes %in% rownames(lakemenonkidney2023_healthy@assays$RNA@scale.data)]

mat <- lakemenonkidney2023_healthy[["RNA"]]@scale.data[selected_genes, selected_cells] %>% as.matrix()

nonzero_counts <- rowSums(mat !=0)
nonzero_proportion <- nonzero_counts/ncol(mat)
row_median <- apply(mat, 1, median)
row_stdev <- apply(mat,1,sd)

filtered_genes <- selected_genes[nonzero_proportion >= 0.1 & abs(row_median) > 0.2] # & row_stdev < 0.8

if(length(filtered_genes)<=1) {
  print(paste("filtered_genes is empty or has insufficient entries; skipping heatmap for", category_l2))
  } else {
    filtered_genes <- filtered_genes[!is.na(filtered_genes)]
    if(length(filtered_genes)==0) {
    print(paste("filtered_genes only contains NA entries; skipping heatmap for", category_l2))
} else {
  mat <- mat[filtered_genes[!is.na(filtered_genes)], ]  

category_name_clean <- gsub("/","_",category_l2)
heatmap_name <- paste("Heatmap_", categories, "_", category_name_clean, ".pdf",sep="")
pdf(heatmap_name)

heatmap <- Heatmap(mat, 
name=category_l2, 
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


####*pulling out - cell types for COLLAGENS ####
inquiry <- c("FIB", "MyoF","aFIB","POD", "NK1")
verifkidney_check <- subset(lakemenonkidney2023_healthy, subset = subclass.l2 %in% collageninquiry)
table(verifkidney_check@meta.data$subclass.l2)

categories = c("Collagens")
cluster_anno <- verifkidney_check@meta.data$subclass.l2
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

mat <- mat[filtered_genes, ]

ht_opt$TITLE_PADDING = unit(c(5,15),"points")
verifkidney_Collagens <- Heatmap(mat, 
                                 name="Fib. & Collagens", 
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
verifkidney_Collagens 
#ggsave("/local1/workdir/agh227/lakemenon2023_kidneyverificationset/id_granular/subsetinterest_MatrisomeCollagens.png")
dev.off()

#####*pulling out - cell types for PROTEOGLYCANS ####
inquiry <- c("FIB", "MyoF","aFIB","POD", "MDC", "ncMON","MAC-M2","MON")
verifkidney_check <- subset(lakemenonkidney2023_healthy, subset = subclass.l2 %in% inquiry)
table(verifkidney_check@meta.data$subclass.l2)

categories = c("Proteoglycans")
cluster_anno <- verifkidney_check@meta.data$subclass.l2
col_fun <- circlize::colorRamp2(c(-1,0,3),c("#FF00FF", "black", "#FFFF00"))

selected_genes <- matrisomelist$Gene.Symbol[matrisomelist$Matrisome.Category %in% categories]
selected_genes <- selected_genes[selected_genes %in% rownames(verifkidney_check@assays$RNA@scale.data)]

mat <- verifkidney_check[["RNA"]]@scale.data[selected_genes, ] %>% as.matrix()

#filtering heatmap entries by presence of expression:
#nonzero_counts <- rowSums(mat !=0)
#nonzero_proportion <- nonzero_counts/ncol(mat)
#row_stdev <- apply(mat,1,sd)
#row_median <- apply(mat, 1, median)

#filtered_genes <- selected_genes[nonzero_proportion >= 0.1 & row_stdev > 0.1 & (row_median > 0.1 | row_median < -0.1)]
#print(length(filtered_genes))

mat <- mat[selected_genes, ]

ht_opt$TITLE_PADDING = unit(c(5,15),"points")
verifkidney_Proteoglycans <- Heatmap(mat, 
                                 name="Proteoglycans", 
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
                                 column_title_rot=90)verifkidney_Proteoglycans 
#ggsave("/local1/workdir/agh227/lakemenon2023_kidneyverificationset/id_granular/subsetinterest_MatrisomeProteoglycans.png")
dev.off()

####*pulling out - cell types for ECM Glycoproteins ####
inquiry <- c("FIB", "MyoF","aFIB","POD", "PEC", "VSMC/P")
verifkidney_check <- subset(lakemenonkidney2023_healthy, subset = subclass.l2 %in% inquiry)
table(verifkidney_check@meta.data$subclass.l2)

categories = c("ECM Glycoproteins")
cluster_anno <- verifkidney_check@meta.data$subclass.l2
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

mat <- mat[filtered_genes, ]

verifkidney_Glycoprot <- Heatmap(mat, 
                                     name="Glycoproteins", 
                                     cluster_columns = T, 
                                     show_column_names = F, 
                                     show_column_dend = F,
                                     cluster_column_slices = F,
                                     column_title_gp = gpar(fontsize=8),
                                     column_split=factor(cluster_anno), 
                                     column_gap = unit(0.6, "mm"), 
                                     cluster_rows = T,
                                     row_names_gp = gpar(fontsize = 8),
                                     show_row_dend = F,
                                     col = col_fun,
                                     use_raster = T, 
                                     raster_quality = 2,
                                     column_title_rot=90)
verifkidney_Glycoprot
#ggsave("/local1/workdir/agh227/lakemenon2023_kidneyverificationset/id_granular/subsetinterest_MatrisomeECMGlycoprot.png", width=30, height=20)
dev.off()

####*pulling out - cell types for ECM Regulators ####
unique(lakemenonkidney2023_healthy@meta.data$subclass.l2) #immune cells had more going on 
inquiry <- c("MDC","NK1","NK2","NKT","MAC-M2","ncMON","B","T","T-CYT","cycT","IC-B")

verifkidney_check <- subset(lakemenonkidney2023_healthy, subset = subclass.l2 %in% inquiry)
table(verifkidney_check@meta.data$subclass.l2)

categories = c("ECM Regulators")
cluster_anno <- verifkidney_check@meta.data$subclass.l2
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

mat <- mat[filtered_genes, ]

ht_opt$TITLE_PADDING = unit(c(5,25),"points")
verifkidney_ECMReg <- Heatmap(mat, 
                                 name="ECM Regulators", 
                                 cluster_columns = T, 
                                 show_column_names = F, 
                                 show_column_dend = F,
                                 cluster_column_slices = F,
                                 column_title_gp = gpar(fontsize=8),
                                 column_split=factor(cluster_anno), 
                                 column_gap = unit(0.6, "mm"), 
                                 cluster_rows = T,
                                 row_names_gp = gpar(fontsize = 8),
                                 show_row_dend = F,
                                 col = col_fun,
                                 use_raster = T, 
                                 raster_quality = 2,
                                 column_title_rot=90)
verifkidney_ECMReg
#subsetinterest_MatrisomeECMReg.png", width=30, height=20)
dev.off()

####*pulling out - cell types for Secreted Factors ####
unique(lakemenonkidney2023_healthy@meta.data$subclass.l2)
inquiry <- c("MDC","NK1","NK2","NKT","MAC-M2","ncMON","B","T","T-CYT","cycT","IC-B") 

verifkidney_check <- subset(lakemenonkidney2023_healthy, subset = subclass.l2 %in% inquiry)
table(verifkidney_check@meta.data$subclass.l2)

cluster_anno <- verifkidney_check@meta.data$subclass.l2
col_fun <- circlize::colorRamp2(c(-1,0,3),c("#FF00FF", "black", "#FFFF00"))

categories = c("Secreted Factors")
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

mat <- mat[filtered_genes, ]

ht_opt$TITLE_PADDING = unit(c(5,25),"points")
verifkidney_SecFac <- Heatmap(mat, 
                              name="Secreted Factors", 
                              cluster_columns = T, 
                              show_column_names = F, 
                              show_column_dend = F,
                              cluster_column_slices = F,
                              column_title_gp = gpar(fontsize=8),
                              column_split=factor(cluster_anno), 
                              column_gap = unit(0.6, "mm"), 
                              cluster_rows = T,
                              row_names_gp = gpar(fontsize = 8),
                              show_row_dend = F,
                              col = col_fun,
                              use_raster = T, 
                              raster_quality = 2,
                              column_title_rot=90)
verifkidney_SecFac
#subsetinterest_MatrisomeSecFac.png", width=30, height=20)
dev.off()

####*pulling out - cell types for ECM-Affiliated ####
unique(lakemenonkidney2023_healthy@meta.data$subclass.l2)
inquiry <- c("PC","MAC-M2","cDC")

verifkidney_check <- subset(lakemenonkidney2023_healthy, subset = subclass.l2 %in% inquiry)
table(verifkidney_check@meta.data$subclass.l2)

cluster_anno <- verifkidney_check@meta.data$subclass.l2
col_fun <- circlize::colorRamp2(c(-1,0,3),c("#FF00FF", "black", "#FFFF00"))

categories = c("ECM-affiliated Proteins")
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

mat <- mat[filtered_genes, ]

ht_opt$TITLE_PADDING = unit(c(5,20),"points")
verifkidney_ECMAffil <- Heatmap(mat, 
                              name="ECM-Affiliated Proteins", 
                              cluster_columns = T, 
                              show_column_names = F, 
                              show_column_dend = F,
                              cluster_column_slices = F,
                              column_title_gp = gpar(fontsize=8),
                              column_split=factor(cluster_anno), 
                              column_gap = unit(0.6, "mm"), 
                              cluster_rows = T,
                              row_names_gp = gpar(fontsize = 8),
                              show_row_dend = F,
                              col = col_fun,
                              use_raster = T, 
                              raster_quality = 2,
                              column_title_rot=90)
verifkidney_ECMAffil
#subsetinterest_ECMAffil.png", width=30, height=20)
dev.off()