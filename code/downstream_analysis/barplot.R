setwd("/media/ggj/ggj/CJY/ob/tissue_all")
library(openxlsx)
library(Seurat)
library(ggplot2)
###barplot ****
dir=list.dirs('.',recursive = F)
dir=dir[c(2,4,6:9)]
dir
anno <- read.xlsx("./pituitary/pituitary_marker_anno.xlsx",sheet=1)
anno <- anno[,c(7,9)]
anno <- na.omit(anno)
new.cluster.ids<-anno$cell.type
Idents(pbmc)=pbmc$seurat_clusters
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
pbmc@meta.data$celltype<-Idents(pbmc)
unique(pbmc$celltype)

pbmc.wt=subset(pbmc,subset=type1=='WT')
pbmc.wt=colnames(pbmc.wt)
pbmc.ob=subset(pbmc,subset=type1=='OB')
pbmc.ob=colnames(pbmc.ob)
pbmc.ob=sample(pbmc.ob,size = 10000)
pbmc.wt=sample(pbmc.wt,size = 10000)
cellname=c(pbmc.ob,pbmc.wt)
Idents(pbmc)=colnames(pbmc)
pbmc=subset(pbmc,idents=cellname)
for (i in 1:length(dir)){
  file.path=list.files(dir[i],pattern = '\\.RData',full.names = T)
  load(file.path)
  Idents(pbmc)=pbmc$celltype
  final<-data.frame(Idents(pbmc),pbmc$type1)
  colnames(final)<-c("celltype","state")
  
  path=paste0(dir[i],'/barplot.pdf')
  ggplot(final,aes(x=celltype))+geom_bar(aes(fill=factor(state)),position = "fill")+xlab("celltype")+ylab("percentage")+
  theme(axis.text.x = element_text(angle=45,vjust = 0.6,face="bold"))+scale_fill_manual(values = c('#f3994c','#6e97c4'))+
    scale_x_discrete(limits=rev(new.cluster.ids))+
  coord_flip()
  
 ggsave(path, device = 'pdf',width = 6,height = 8)
  message(dir[i])
}

