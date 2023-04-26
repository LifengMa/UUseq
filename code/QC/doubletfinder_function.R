library(DoubletFinder)
library(Seurat)
library(dplyr)
library(Seurat)
library(ggplot2)
rm_doublet <- function(seuratobject=data,dim.usage=20,res=0.8,DoubletRate=0.05,name='seurat',
                       outDir='/media/ggj/figures',
                       saveDir='/media/ggj/ggjlab2/UUseq/OB/Doubletfinder/') {
  
  if(!file.exists(saveDir)){
    dir.create(saveDir)
  }
  if(!file.exists(outDir)){
    dir.create(outDir)
  }
  
  seuratobject <- CreateSeuratObject(seuratobject,min.cells=3,project = name) %>%
        NormalizeData() %>%
        FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
        ScaleData() %>%
        RunPCA() %>%
        FindNeighbors() %>%
        FindClusters(resolution = res) %>%
        RunUMAP(dims = 1:dim.usage)
  
  ##DoubletFinder去双胞的标准流程，封装成一个函数
  Find_doublet <- function(data){
    sweep.res.list <- paramSweep_v3(data, PCs = 1:dim.usage, sct = FALSE)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    homotypic.prop <- modelHomotypic(data$seurat_clusters)
    nExp_poi <- round(DoubletRate*ncol(data))
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    p<-as.numeric(bcmvn$pK[which.max(bcmvn$BCmetric)])
    # data <- doubletFinder_v3(data, PCs = 1:dim.usage, pN = 0.25, pK = p, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
    # pANNuse<-colnames(data@meta.data)[grep(paste0("pANN_0.25_",as.numeric(pKuse)),colnames(data@meta.data))]
    # pANNuse
    data <- doubletFinder_v3(data, PCs = 1:dim.usage, pN = 0.25, pK = p, nExp = nExp_poi.adj, reuse.pANN =FALSE, sct = FALSE)
    colnames(data@meta.data)[ncol(data@meta.data)] = "doublet_info"
    return(data)
  }
  
  ##调用上面写好的函数，返回的是一个Seurat对象，meta.data信息里会有双胞信息，需要自己手动删除
  seuratobject<-Find_doublet(seuratobject)
  p1 <- DimPlot(seuratobject,group.by = 'doublet_info')+ggtitle(paste0(name,"_doublet_info"))
  p2 <- DimPlot(seuratobject,group.by = 'seurat_clusters',label=T)+ggtitle(paste0(name,"_seurat_clusters"))
  p1|p2
  ggsave(paste0(outDir,"/",name,"_Doubletfinder_result.png"),width = 15,height = 6)
  seuratobject<-subset(seuratobject,subset=doublet_info=="Singlet")
  #seuratobject@meta.data$library = name #顺便打上这个样本的label
  
  ##生成的Seurat对象有个问题，会在meta.data里多了很多pANN_开头的列，需要手动删除
  c <- grep("pANN_",colnames(seuratobject@meta.data))
  seuratobject@meta.data <- seuratobject@meta.data[,-c]
  saveRDS(seuratobject,file=paste0(saveDir,"/",name,"_Doubletfinder.rds"))
  ##输出此样本的细胞数
  print(paste0(" cells: ", length(seuratobject@meta.data$orig.ident)," is remained!"," Time:",format(Sys.time(), "%Y%m%d %X"),sep = " "))
  return(seuratobject)
}