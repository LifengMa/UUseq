setwd("/media/ggj/ggjlab2/UUseq/brain/Doubletfinder/")

library(Seurat)
library(reshape2)
library(scDblFinder)
library(SingleCellExperiment)
library(BiocParallel)

seob <- readRDS("../brain_rawdata_1000gene.rds")
#seob$index <- colsplit(colnames(seob),"_",names = c("n1","n2"))$n1
seob=as.SingleCellExperiment(seob)

seob=scDblFinder(seob,samples = 'index',BPPARAM = MulticoreParam(8)) 
seob=as.Seurat(seob)

table(seob$scDblFinder.class)
dir.create('../scDblFinder')
saveRDS(seob,"../scDblFinder/rmdoublet_result_0403.rds")
