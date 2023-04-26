library(openxlsx)
library(Seurat)
library(ggplot2)
library(ggrepel)

anno <- read.xlsx("./ovary/ovary_marker_anno.xlsx",sheet=1)
anno <- anno[,c(7,10)]
anno <- na.omit(anno)
new.cluster.ids<-anno$celltype
Idents(pbmc)=pbmc$seurat_clusters
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
pbmc@meta.data$celltype<-Idents(pbmc)
unique(pbmc$celltype)

Idents(pbmc) <- pbmc$celltype
celltype=unique(pbmc$celltype)
celltype=celltype[c(1:6,8,11:12)]
celltype
library(future)
plan('multisession',workers=12)
options(future.globals.maxSize=50000*1024^2)
##create deg dataframe
deg <- NULL
for(i in 1:length(celltype)){
  temp <- subset(pbmc,idents=celltype[i])
  sum <- as.data.frame(table(temp$type1))
  Idents(temp) <- temp$type1
  if(length(unique(temp$type1)) == 2 & length(which(temp$type1 %in% 'OB'))>3 & length(which(temp$type1 %in% 'WT'))>3){
    deg.temp <- FindMarkers(temp,ident.1 = 'OB',ident.2 = 'WT',logfc.threshold = 0)
    #deg.temp <- deg.temp[deg.temp$p_val_adj<0.05,]
    deg.temp$type <- ifelse(deg.temp$avg_log2FC>0,"OB","WT")
    deg.temp$celltype <- celltype[i]
    deg.temp$gene <- rownames(deg.temp)
    deg <- rbind(deg,deg.temp)
  }
  else{
    next
  }
}
###绘制火山图
for (i in 1:length(celltype)) {
temp=deg[deg$celltype==celltype[i],]
temp=temp[order(-abs(temp$avg_log2FC)),]
temp=temp[!duplicated(temp$gene),]
temp$rank=paste0(1:length(rownames(temp)))
temp$rank=as.numeric(temp$rank)
temp$change=factor(ifelse(abs(temp$avg_log2FC)>1,ifelse(temp$avg_log2FC>1,'OB','WT'),'mid'))
temp$log10p_val=-log10(temp$p_val)
#temp$text=ifelse(abs(temp$avg_log2FC)>1 & temp$log10p_val>35,as.character(temp$gene),'')
temp$text=ifelse(temp$rank<21,as.character(temp$gene),'')

ggplot(temp,aes(avg_log2FC,log10p_val,color=change))+geom_point()+
  scale_y_continuous(limits = c(0,375))+
  scale_color_manual(values=c('grey','#546de5','#ff4757'))+
  geom_vline(xintercept = c(-1,1),lty=4,col='black',lwd=0.8)+
  labs(x='log2FC',y='-log10 p_value')+
  geom_text_repel(aes(label=temp$text),color='black',size=3, box.padding=unit(0.5,'lines'),max.overlaps = 50 ,point.padding = NA)

ggsave(paste0("../basic_analysis/diffgene/ovary/",celltype[i],'.png'),width = 6,height = 6)

}


#####GO####
library(org.Mm.eg.db)    
library(ggplot2)   
library(clusterProfiler)
library(DOSE)
# dir.create("../basic_analysis/GO_KEGG")
# deg$celltype <- gsub("Lac/Som","Lac|Som",deg$celltype)
# temp <- OB.deg[1:30,]
# temp.WT <- WT.deg[1:19,]

#####OB
for (i in 1:length(celltype)) {
  ob.temp <- deg[deg$celltype==celltype[i],]
  ob.temp=ob.temp[order(-ob.temp$avg_log2FC),]
  gene.OB<- ob.temp$gene[1:100]
  eg <- bitr(gene.OB, fromType="SYMBOL", 
             toType="ENTREZID", OrgDb="org.Mm.eg.db")
  go <- enrichGO(eg$ENTREZID, OrgDb = org.Mm.eg.db, ont='BP',pAdjustMethod = 'BH',
                 pvalueCutoff = 0.05, qvalueCutoff = 0.05) 
  dotplot(go,showCategory=20)+
    theme(axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = 7))+ggtitle(paste0("OB"))
  ggsave(paste0("../basic_analysis/GO_KEGG/ovary/GO_OB_",celltype[i],".png"),width = 7,height = 7)
}
  
  


  #barplot(go,showCategory=20,font.size = 12)
 
  
#######WT
for (i in 1:length(celltype)) {
  wt.temp <- deg[deg$celltype==celltype[i]&deg$type=='WT',]
  wt.temp=wt.temp[order(wt.temp$avg_log2FC),]
  gene.WT<- wt.temp$gene[1:100]
  eg <- bitr(gene.WT, fromType="SYMBOL", 
             toType="ENTREZID", OrgDb="org.Mm.eg.db")
  go <- enrichGO(eg$ENTREZID, OrgDb = org.Mm.eg.db, ont='BP',pAdjustMethod = 'BH',
                 pvalueCutoff = 0.05, qvalueCutoff = 0.05) 
  dotplot(go,showCategory=20)+
    theme(axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = 7))+ggtitle(paste0("WT"))
  ggsave(paste0("../basic_analysis/GO_KEGG/ovary/GO_WT_",celltype[i],".png"),width = 7,height = 7)
}
  
  ###KEGG
  #OB
  kk.OB <- enrichKEGG(gene = eg$ENTREZID,
                   organism ="mmu",
                   pvalueCutoff = 0.1,
                   qvalueCutoff = 0.2,
                   minGSSize = 1,
                   use_internal_data =FALSE)
  
  #barplot(kk,showCategory=20,font.size = 12)
  dotplot(kk.OB,showCategory=20)+
    theme(axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = 7))+ggtitle(paste0("WT"))
  ggsave("../basic_analysis/GO_KEGG/GO_WT.pdf",width = 10,height = 10)

  #WT
  kk.WT <- enrichKEGG(gene = eg.WT$ENTREZID,
                      organism ="mmu",
                      pvalueCutoff = 0.1,
                      qvalueCutoff = 0.2,
                      minGSSize = 1,
                      use_internal_data =FALSE)
  
  #barplot(kk,showCategory=20,font.size = 12)
  dotplot(kk.WT,showCategory=20)+
    theme(axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = 7))+ggtitle(paste0("WT"))
  ggsave("../basic_analysis/GO_KEGG/GO_WT.pdf",width = 10,height = 10)
  
  
  save.image('/media/ggj/ggj/CJY/ob/basic_analysis/diffgene/diffgene.RData')
  
