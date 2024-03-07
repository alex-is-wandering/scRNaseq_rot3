library(Seurat)
library(dplyr)
library(patchwork)
library(ComplexHeatmap)
library(ggplot2)

###### READING IN DATA ###### 
setwd("/workdir/data/crosstissue_data/kidney")
parenchymal<-readRDS("2020-06-07.human_baseline_paren_v3.RDS")
reactome_ECMOrg <- read.delim("/workdir/agh227/reference_ECM_signatures/reactome_ECMOrganization.tsv")
#genes at "reactome_ECMOrg$associations.gene.symbol"

#head(parenchymal@meta.data,5) 
#set ident to celltype, for granular vs high-level
Idents(object = parenchymal) <- parenchymal$celltype
DimPlot(parenchymal,reduction="umap")

###### ECM EXPLORATION ###### 

selected_genes <- reactome_ECMOrg$associations.gene.symbol
#*we only want to map those genes that are named in the matrisomelist
selected_genes <- selected_genes[selected_genes %in% rownames(parenchymal@assays$RNA@data)]

DotPlot(parenchymal,features=selected_genes) + theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave("/workdir/agh227/kidneydata_figures/parenchymal/Reactome_ECMOrg/paren_dotplot_alltypes_allgenes.png",width=25,height=10) #not sure why background is wonky, but it works
dev.off()


###### dotplot & vlnplot: looking at celltypes vs individual collagen subunits ###### 
#VIOLIN PLOT PER-GENE
VlnPlot(parenchymal,features="COL10A1",slot="data")
VlnPlot(parenchymal,features="COL1A1",slot="data")
ggsave("/workdir/agh227/kidneydata_figures/paren_COL1A1_violin.png",width=10,height=5)
dev.off()

#DOT PLOTS FOR Collagens in ECM Organization
paren_geneIDs <- rownames(parenchymal@assays$RNA@data) #not scale.data
COL_genes <- paren_geneIDs[grepl("COL",paren_geneIDs,ignore.case=FALSE)]

selected_genes <- reactome_ECMOrg$associations.gene.symbol
selected_genes <- selected_genes[selected_genes %in% COL_genes]

DotPlot(parenchymal,features=selected_genes) + theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave("/workdir/agh227/kidneydata_figures/parenchymal/Reactome_ECMOrg/paren_dotplot_alltypes_collagens.png",width=25,height=10)
dev.off()

###### Seurat AddModuleScore: looking at celltypes vs ECM component categories ######

# We only have 1 category here, so:
selected_genes <- reactome_ECMOrg$associations.gene.symbol

selected_genes <- selected_genes[selected_genes %in% rownames(parenchymal@assays$RNA@data)] 

parenchymal <- AddModuleScore(parenchymal, features = list(selected_genes), name = "ECMOrg")

#saveRDS(parenchymal, file="/workdir/agh227/processeddata/parenchymal_ECMOrgScored.RDS")

head(parenchymal@meta.data,5)

# next:
VlnPlot(parenchymal, features="ECMOrg1")
ggsave("/workdir/agh227/kidneydata_figures/parenchymal/Reactome_ECMOrg/paren_ViolinECMOrg_AddModuleScored.png",width=20,height=10)
dev.off()

###### Seurat DoHeatmap: gene expression vs ECM organization ###### 
unique_celltypes <- unique(parenchymal$celltype)
head(reactome_ECMOrg,5)

#*naming which categories to map
#*specifying geneID location
selected_genes <- reactome_ECMOrg$associations.gene.symbol
selected_genes <- selected_genes[selected_genes %in% rownames(parenchymal@assays$RNA@data)]

#*making heatmaps: ALL CORE MATRISOME COMPONENTS (3) categories, all cell types
DoHeatmap(parenchymal, features=selected_genes, slot="data")
ggsave("/workdir/agh227/kidneydata_figures/parenchymal/Reactome_ECMOrg/paren_heatmap_alltypes_ECMOrg.png",width=20,height=10)
dev.off()