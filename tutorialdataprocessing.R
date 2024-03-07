setwd("/local/workdir/crosstissue_data/kidney")

library(Seurat)
library(dplyr)
library(patchwork)

load("2019-10-02.human_sub.RObj")

head(human_sub@meta.data,10) #so yes, this does seem to be annotated
updated_humansub=UpdateSeuratObject(human_sub)
head(updated_humansub@meta.data,5)

percent_mitoDNA=VlnPlot(updated_humansub,features=c("nFeature_RNA","nCount_RNA","percent.mito"),ncol=3)
percent_mitoDNA=VlnPlot(updated_humansub,features=c("percent.mito"),ncol=1)

plot1<-FeatureScatter(updated_humansub,feature1="nCount_RNA",feature2="percent.mito")
plot2 <- FeatureScatter(updated_humansub, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

updated_humansub<-subset(updated_humansub,subset=nFeature_RNA>200 & nFeature_RNA<2500 & percent.mito<5) #this cutoff for mito% is arbitrary

updated_humansub<-NormalizeData(updated_humansub)

updated_humansub<-FindVariableFeatures(updated_humansub, selection.method="vst",nfeatures=2000)
top10 <- head(VariableFeatures(updated_humansub), 10)

plot1<-VariableFeaturePlot(updated_humansub)
plot2<-LabelPoints(plot=plot1,points=top10,repel=TRUE)
plotlayout<-plot1+plot2+plot_layout(guides='collect')


all.genes<-rownames(updated_humansub)
updated_humansub<-ScaleData(updated_humansub,features=all.genes) #"centering and scaling data matrix"
updated_humansub<-ScaleData(updated_humansub,vars.to.regress="percent.mito") #"regressing out percent.mito" then once more "centering and scaling data matrix"

updated_humansub<-RunPCA(updated_humansub,features=VariableFeatures(object=updated_humansub))

print(updated_humansub[["pca"]],dims=1:5,nfeatures=5)

VizDimLoadings(updated_humansub,dims=1:2,reduction="pca")

DimHeatmap(updated_humansub,dims=1,cells=500,balanced=TRUE)
DimHeatmap(updated_humansub,dims=1:15,cells=500,balanced=TRUE)

ElbowPlot(updated_humansub)

updated_humansub <- FindNeighbors(updated_humansub, dims=1:10)
updated_humansub <- FindClusters(updated_humansub, resolution=0.5)

head(Idents(updated_humansub),5)

updated_humansub<-RunUMAP(updated_humansub,dims=1:10)
DimPlot(updated_humansub,reduction="umap")