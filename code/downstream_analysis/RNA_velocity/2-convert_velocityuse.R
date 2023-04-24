setwd("/media/ggj/ggjlab2/UUseq/brain/benchmark/brain/10X/")
library(SeuratDisk)
library(reshape2)
library(Seurat)
load("./brain.RData")

anno <- pbmc@meta.data
use <- pbmc
head(anno)
anno <- anno[colnames(use),]
#use$celltype <- anno$celltype_l7

# cellname <- colsplit(colnames(use),"_",names=c("n1","n2","n3","n4","n5"))
# length(unique(cellname$n4))
# cellname$n4 <- gsub("-1","",cellname$n4)

# RNAvelocity_index <- read.csv("./发表文章整理/alldata/RNAvelocity/velocity_index.csv")
# length(unique(RNAvelocity_index$X))
# 
# length(intersect(cellname$n4,RNAvelocity_index$X))

# use$celltype <- gsub("B_12_LMO2_LZGC","B_c12_LMO2_LZGC",use$celltype)
# use$celltype <- gsub("B_11_CXCR4_DZGC","B_c11_CXCR4_DZGC",use$celltype)
# use$celltype <- gsub("B_10_PSME2_PreGC","B_c10_PSME2_PreGC",use$celltype)
# use$celltype <- gsub("B_02_IFIT3_B","B_c02_IFIT3_B",use$celltype)
# use$celltype <- gsub("B_01_TCL1A_naïveB","B_c01_TCL1A_naïveB",use$celltype)
SaveH5Seurat(use,filename="/media/ggj/Guo-4T-AB3/UU/OB/brain/RNAvelocity/10X_brain.h5Seurat",overwrite=T)
Convert("/media/ggj/Guo-4T-AB3/UU/OB/brain/RNAvelocity/10X_brain.h5Seurat",dest="h5ad",overwrite = TRUE)

write.table(Embeddings,file = "/media/ggj/ggjlab2/hezuo/mjq/RNAvelocity/eb_strict.csv",quote = F,col.names = T,sep = "\t")
#write.table(col,file = "/media/ggj/ggjlab2/hezuo/mjq/RNAvelocity/col.csv",quote = F,sep = "\t")
write.table(anno,file = "/media/ggj/ggjlab2/hezuo/mjq/RNAvelocity/anno_strict.csv",quote = F,col.names = T,sep = "\t")
