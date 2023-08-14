setwd("/media/ggj/ggjlab2/UUseq/OB/Doubletfinder/save/")
library(data.table)
library(reshape2)
library(Seurat)
library(ggplot2)
library(Seurat)
rm(list=ls())
gc()
tissue <- 'spinal_cord'
#data <- readRDS("../../rawdata/OB_WT_alldataset_filter_UMI500more.rds")
#anno <- data@meta.data
#write.csv(anno,file="../../rawdata/OB_WT_alldataset_filter_UMI500more_anno.csv")
anno <- read.csv("../../rawdata/OB_WT_alldataset_filter_UMI500more_anno.csv",row.names = 1)
file <- list.files(pattern = tissue) 
combined <- readRDS(file[1])
for(i in 2:length(file)){
  temp <- readRDS(file[i])
  combined <- merge(combined,temp)
}

dim(combined)
anno <- anno[colnames(combined),]
combined$tissue <- anno$tissue
combined$sample <- anno$sample
combined$type <- anno$type
combined$reads <- anno$reads

table(combined$tissue)
table(combined$sample)


pbmc<-CreateSeuratObject(counts =GetAssayData(combined@assays$RNA,slot="counts"), min.cells = 3,  
                         project = tissue,meta.data = combined@meta.data)
dim(pbmc)
## Quality Control
dge<-as.data.frame(as.matrix(GetAssayData(pbmc,slot = "counts")))
summary(colSums(combined))#3729   mean UMI
hist(colSums(dge),breaks = 100)
abline(v=500)
hist(colSums(dge>0),breaks = 100)
abline(v=200)
summary(colSums(combined>0))#1437  mean gene
mito.genes <- grep(pattern = "^mt-", x = rownames(x = pbmc@assays$RNA), value = TRUE)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
par(mfrow = c(1, 2))
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))


pbmc <- subset(pbmc, subset = nFeature_RNA >=200 &nFeature_RNA < 5000 & percent.mt < 40 &
                 nCount_RNA < 20000 &nCount_RNA>=500)

dim(pbmc)#26327 92164

##### normalization
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", 
                      scale.factor = 10000)
dim(pbmc@assays$RNA)

## Find Variable Genes
pbmc <- FindVariableFeatures(object = pbmc )
hv.genes <- head(VariableFeatures(pbmc),2000)
#pbmc <- FindVariableFeatures(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, 
#                          x.low.cutoff = 0.015, y.cutoff = 0.4)
length(VariableFeatures(pbmc))#1869
a<-setdiff(VariableFeatures(pbmc),mito.genes)
hv.genes <-a
length(hv.genes)#1864

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

##### scale data
pbmc <- ScaleData(object =pbmc,features = hv.genes, vars.to.regress = c("percent.mt","nCount_RNA"))

##### runPCA
pbmc <- RunPCA(object =pbmc, features = hv.genes, pcs.compute = 50, do.print = TRUE, 
               pcs.print = 1:5, genes.print = 5)
ElbowPlot(object =pbmc, ndims = 40)#1:25

# run TSNE and FindCluster
pbmc <- FindNeighbors(pbmc,dims=1:40)
pbmc <- FindClusters(pbmc, resolution = 0.8)

# pbmc <- RunTSNE(object =pbmc, reduction = "pca", dims= 1:20,
#                 nthreads = 4, reduction.name = "FItSNE", reduction.key = "FItSNE_", fast_tsne_path = "/home/ggj/Documents/tools/fit-tsne/FIt-SNE-master/bin/fast_tsne", 
#                 max_iter = 2000)   #fast_tsne_path 
pbmc <- RunTSNE(object =pbmc, reduction = "pca", dims= 1:40,
                reduction.name = "tSNE") 
pbmc <- RunUMAP(object =pbmc, reduction = "pca", dims= 1:40,
                reduction.name = "umap") 

library(RColorBrewer)
col_flg<-colorRampPalette(brewer.pal(8,"Paired"))(length(levels(as.factor(pbmc$seurat_clusters))))
p1 <- DimPlot(object = pbmc, raster=FALSE,reduction= "umap",label  = T , pt.size = 0.1,
              label.size = 4,cols = col_flg) 
p1
DimPlot(object = pbmc, reduction= "umap",label  = T , pt.size = 0.1,label.size = 5,cols = col_flg) 
#DimPlot(object = pbmc, reduction= "tSNE",label  = T , pt.size = 0.1,label.size = 5,cols = col_flg) 

dir.create(paste0("../../tissue/tissue_all_type/",tissue))
ggsave(paste0("../../tissue/tissue_all_type/",tissue,"/",tissue,"_cluster.pdf"),width = 10,height = 8)
# Batch
# library(reshape2)
# cellname<-colnames(pbmc@assays$RNA)
# cellname1<-colsplit(cellname,pattern = "\\.",names = c("Batch","n2"))
# pbmc$Batch<-cellname1$Batch
# table(cellname1$Batch)
# p2 <- DimPlot(object = pbmc, reduction = "tSNE", label  = F ,
#               pt.size = 0.1,group.by = "Batch")
# p2


p3 <- DimPlot(object = pbmc, reduction = "umap", label  = F ,
              pt.size = 0.1,group.by = "sample")
p3
ggsave(paste0("../../tissue/tissue_all_type/",tissue,"/",tissue,"_sample.pdf"),width = 10,height = 8)

pbmc$type <- pbmc$sample
pbmc$type <- gsub("[0-9]","",pbmc$type)
p4 <- DimPlot(object = pbmc, reduction = "umap", label  = F ,
              pt.size = 0.1,split.by = "type",cols = col_flg)
p4
ggsave(paste0("../../tissue/tissue_all_type/",tissue,"/",tissue,"_type.pdf"),width = 15,height = 8)

library(future)
plan(multisession, workers=10)
options(future.globals.maxSize= 1500000000)
pbmc.markers <- FindAllMarkers(object =pbmc, only.pos = TRUE, min.pct = 0.25, 
                               logfc.threshold  = 0.25)
pbmc.markers<-pbmc.markers[order(pbmc.markers$cluster,-pbmc.markers$avg_log2FC,pbmc.markers$p_val  ),]
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = pbmc, features  = top10$gene,  label = T, size=5)

#check some genes 
FeaturePlot(pbmc,features = "Runx2")
DotPlot(pbmc,features = "FABP4")
VlnPlot(pbmc,features = "FABP4")
######
load("./brainstem/brainstem.RData")
anno <- openxlsx::read.xlsx("./brainstem/brainstem_total_marker_version1.xlsx")
anno <- anno[,c(7,9)]
anno <- na.omit(anno)
anno$cell.type <- paste0("C",anno$cluster,"_",anno$cell.type)
#anno$cell.type <- gsub("adipo","Adipo",anno$cell.type)

new.cluster.ids<-anno$cell.type
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
pbmc@meta.data$celltype<-Idents(pbmc)
library(RColorBrewer)
col_flg<-colorRampPalette(brewer.pal(8,"Paired"))(length(levels(as.factor(pbmc$celltype))))
p5 <- DimPlot(object = pbmc, reduction = "umap",label  = T ,cols = col_flg  ,pt.size = 0.5)
p5
ggsave(paste0("../../tissue/tissue_all_type/",tissue,"/",tissue,"_anno.pdf"),width = 15,height = 8)

Idents(pbmc) <- pbmc$celltype
p6 <- DimPlot(object = pbmc, reduction = "umap", label  = F ,
              pt.size = 0.5,split.by = "type",cols = col_flg)
p6
ggsave(paste0("../../tissue/tissue_all_type/",tissue,"/",tissue,"_anno_split.pdf"),width = 15,height = 8)

#save
save(pbmc,pbmc.markers, file = paste0("../../tissue/tissue_all_type/",tissue,"/",tissue,".RData"))
WriteXLS::WriteXLS(pbmc.markers, paste0("../../tissue/tissue_all_type/",tissue,"/",tissue,"_marker.xlsx"),row.names=T)

######harmony
library(harmony)
pbmc <- harmony::RunHarmony(pbmc,"sample", plot_convergence = TRUE)
ElbowPlot(pbmc, ndims = 50,reduction = "harmony")

pbmc <- RunUMAP(pbmc,reduction = "harmony", dims = 1:20) 
pbmc <- FindNeighbors(pbmc,reduction = "harmony", dims = 1:20) 

pbmc <- FindClusters(pbmc,resolution =0.8)
col_flg<-colorRampPalette(brewer.pal(8,"Paired"))(length(levels(as.factor(pbmc$seurat_clusters))))
p1 <- DimPlot(object = pbmc, raster=FALSE,reduction= "umap",label  = T , pt.size = 0.5,label.size = 4,cols = col_flg) 
p1
DimPlot(object = pbmc, reduction= "umap",label  = T , pt.size = 0.1,label.size = 5,cols = col_flg) 

p3 <- DimPlot(object = pbmc, reduction = "umap", label  = F ,
              pt.size = 0.1,group.by = "sample")
p3

p4 <- DimPlot(object = pbmc, reduction = "umap", label  = F ,
              pt.size = 0.1,split.by = "type",cols = col_flg)
p4
