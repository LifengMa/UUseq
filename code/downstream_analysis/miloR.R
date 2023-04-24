setwd("/media/ggj/ggjlab2/UUseq/OB/")

library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
library(Seurat)
library(openxlsx)
library(SeuratDisk)

load("./tissue/tissue_all_type/heart/heart.RData")
# DimPlot(pbmc)
# anno <- read.xlsx("./tissue/heart/heart_marker_anno.xlsx",sheet=1)
# anno <- anno[,c(7,9)]
# anno <- na.omit(anno)
# 
# new.cluster.ids<-anno$cell.type
# names(new.cluster.ids) <- levels(pbmc)
# pbmc <- RenameIdents(pbmc, new.cluster.ids)
# pbmc@meta.data$celltype<-Idents(pbmc)
# save(pbmc,pbmc.markers,file="./tissue/heart/heart.RData")
dim(pbmc)
pbmc.sce <- as.SingleCellExperiment(pbmc) 
pbmc.sce

plotReducedDim(pbmc.sce,colour_by = "type",dimred = "UMAP")

######Differential abundance testing
#Construct KNN graph
pbmc.milo <- Milo(pbmc.sce)
#use the same value for k used for KNN graph building for clustering by Seurat
k <- pbmc@commands$FindNeighbors.RNA.pca$k.param
# use the same value for d used for KNN graph building
d <- length(pbmc@commands$FindNeighbors.RNA.pca$dims)
pbmc.milo <- buildGraph(pbmc.milo,k=k
                        ,d=d,reduced.dim = "PCA")

#Defining representative neighbourhoods on the KNN graph
{
if(length(colnames(pbmc))>60000){
  prop <- 0.05
}
else{
  prop <- 0.1
}
}
pbmc.milo <- makeNhoods(pbmc.milo,prop = prop,k=k,d=d,refined = TRUE,reduced_dims = "PCA")
plotNhoodSizeHist(pbmc.milo) ###50~100

#Counting cells in neighbourhoods
pbmc.milo <- countCells(pbmc.milo,meta.data=as.data.frame(colData(pbmc.milo)),sample="sample")
head(nhoodCounts(pbmc.milo))

#Defining experimental design
pbmc_design <- data.frame(colData(pbmc.milo))[,c("sample","type")]
pbmc_design <- distinct(pbmc_design)
rownames(pbmc_design) <- pbmc_design$sample

pbmc_design

#Computing neighbourhood connectivity
pbmc.milo <- calcNhoodDistance(pbmc.milo,d=d,reduced.dim = "PCA") 

#Testing
da_result <- testNhoods(pbmc.milo,design=~type,design.df = pbmc_design,reduced.dim = "PCA")
da_result%>%
  arrange(SpatialFDR) %>%
  head()

#####Inspecting DA testing results
ggplot(da_result,aes(PValue))+geom_histogram(bins=50)
ggplot(da_result,aes(logFC,-log10(SpatialFDR)))+
  geom_point()+
  geom_hline(yintercept=1)

pbmc.milo <- buildNhoodGraph(pbmc.milo)
umap_pl <- plotReducedDim(pbmc.milo,dimred = 'UMAP',
                          colour_by = "type",text_by="seurat_clusters",
                          text_size = 3,point_size=0.5)+
  guides(filee="none")
nh_graph_pl <- plotNhoodGraphDA(pbmc.milo,da_result,layout="UMAP",alpha=0.5)

pdf("./tissue/heart/miloR.pdf",width = 15,height = 8)
umap_pl+nh_graph_pl+
  plot_layout(guides="collect")
dev.off()

da_result <- annotateNhoods(pbmc.milo,da_result,coldata_col = "seurat_clusters")
da_result$seurat_clusters <- paste0("C",da_result$seurat_clusters)
head(da_result)
ggplot(da_result,aes(seurat_clusters_fraction))+geom_histogram(bins=50)
da_result$seurat_clusters <- ifelse(da_result$seurat_clusters_fraction<0.6,"Mixed",da_result$seurat_clusters)
plotDAbeeswarm(da_result,group.by = "seurat_clusters",alpha = 0.5)
ggsave("./tissue/heart/miloR_celltype_fraction.pdf",width = 8,height = 8)

####Finding markers of DA populations
#Automatic grouping of neighbourhoods
pbmc.milo <- buildNhoodGraph(pbmc.milo)
da_result <- groupNhoods(pbmc.milo,da_result,max.lfc.delta = 10,da.fdr = 0.1)
head(da_result)
plotNhoodGroups(pbmc.milo,da_result,layout="UMAP")
plotDAbeeswarm(da_result,"NhoodGroup",alpha = 0.5)

rm(list=c("pbmc","pbmc.markers","pbmc.sce"))
save.image("./tissue/heart/miloR.RData")
gc()
