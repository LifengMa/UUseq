setwd("/media/ggj/ggj/CJY/ob/doubletfinder/")
combined <- readRDS("../OB_WT_alldataset_filter_UMI500more.rds")
#combined=data
library(Seurat)
library(reshape2)
source("./doubletfinder_function.R")

#combined$type <- gsub("_8w|_11w","",combined$type)
combined$index=barcode$index
barcode=colsplit(combined$barcode,"\\.",names = c("index","cb"))
combined$index=gsub('col','',combined$index)
#combined$index.sample=paste(combined$sample,combined$index,sep='_')


sample=as.character(unique(Idents(combined)))
sample=grep('WT1|WT2|OB1|OB2',sample,value = T)

type <- as.character(unique(Idents(combined)))
type <- grep('OB9',type,value = T)
table(combined$tissue,combined$sample)
for(i in 1:length(type)){
  temp <- subset(batch1,idents=type[i])
  temp <- GetAssayData(temp@assays$RNA,slot="counts")
  temp <- rm_doublet(seuratobject =temp ,dim.usage = 20,res=0.8,name=paste0('batch1_',type[i]),
                     saveDir='/media/ggj/ggj/CJY/ob/doubletfinder/save',
                     outDir='/media/ggj/ggj/CJY/ob/doubletfinder/figures'
                     )
}

#list < list.files(pattern = '*.rds')
path = '/media/ggj/ggj/CJY/ob/doubletfinder/save'
fileName = dir(path)
for(i in 1:length(fileName)){
  rds <- readRDS(file = paste(path, fileName[i],sep='/'))
  message(paste(fileName[i],' left ', length(colnames(rds))))
}
