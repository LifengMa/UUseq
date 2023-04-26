install.packages("sva")
devtools::install_local("~/Downloads/CytoTRACE_0.3.3.tar.gz")
setwd('/media/ggj/ggj/CJY/ob/basic_analysis/')
###################################################### cytotrace
library(CytoTRACE)
library(ggplot2)
library(patchwork)
library(Seurat)
marrow_10x_expr[1:4,1:4]

#load("/media/ggj/ggj/CJY/ob/tissue_all/pituitary/pituitary_anno.RData")
### cytotrace
setwd('/media/ggj/ggj/CJY/ob/tissue_all/')
dir=list.dirs('.',recursive = F)
dir=dir[c(2,4,6:9)]
dir

p=list()
for(i in 1:length(dir)){
  file.path=list.files(dir[i],pattern = '\\.RData',full.names = T)
  load(file.path)
  
  Idents(pbmc)=colnames(pbmc)
idents=Idents(pbmc)
idents_small=sample(idents, size = 10000)
pbmc_small=subset(pbmc,idents=idents_small)
dataset <- as.data.frame(as.matrix(pbmc_small@assays$RNA@counts))
rm(pbmc)
gc()

result <- CytoTRACE(dataset,ncores = 10)
gcs=as.data.frame(result$GCS)
gcs$type=pbmc_small$type1
colnames(gcs)[1]='GCS'
p[[i]]=ggplot(gcs,aes(x=type,y=GCS,color=type))+geom_boxplot()+
  scale_color_manual(values = c('#f3994c','#6e97c4'))
rm(dataset)
gc()
}

p[[1]]+ theme( panel.background = element_blank() )+p[[2]]+ theme( panel.background = element_blank() )+
p[[3]]+ theme( panel.background = element_blank() )+p[[4]]+ theme( panel.background = element_blank() )+
p[[5]]+ theme( panel.background = element_blank() )+p[[6]]+ theme( panel.background = element_blank() )+plot_layout(guides="collect")
ggsave('/media/ggj/ggj/CJY/ob/basic_analysis/cytotrace_SLICE/cytotrace.pdf',width = 12,height = 8)

p2=list()
for(i in 1:length(dir)){
  file.path=list.files(dir[i],pattern = '\\.RData',full.names = T)
  load(file.path)
  
  Idents(pbmc)=colnames(pbmc)
  idents=Idents(pbmc)
  idents_small=sample(idents, size = 10000)
  pbmc_small=subset(pbmc,idents=idents_small)
  dataset <- as.data.frame(as.matrix(pbmc_small@assays$RNA@counts))
  rm(pbmc)
  gc()
  
  result <- CytoTRACE(dataset,ncores = 10)
  cytotrace=as.data.frame(result$CytoTRACE)
  cytotrace$type=pbmc_small$type1
  colnames(cytotrace)[1]='CytoTRACE'
  p2[[i]]=ggplot(cytotrace,aes(x=type,y=CytoTRACE,color=type))+geom_boxplot()+
    xlab(paste0(gsub('./','',dir[i])))+theme( panel.background = element_blank() )+
    scale_color_manual(values = c('#f3994c','#6e97c4'))
  rm(dataset)
  gc()
}
p2[[1]]+p2[[2]]+p2[[3]]+p2[[4]]+p2[[5]]+p2[[6]]+plot_layout(guides="collect")
ggsave('/media/ggj/ggj/CJY/ob/basic_analysis/cytotrace_SLICE/cytotrace_2.pdf',width = 12,height = 8)

# cytotrace_pheno <- read.csv('./intestine_cytotrace.csv',stringsAsFactors = F,row.names = 1)
# cytotrace_pheno <- cbind(pheno[,1:3],cytotrace_pheno[rownames(pheno),])
# 
# cytotrace_pheno$pheno <- factor(cytotrace_pheno$pheno,levels = as.character(zz$pheno))
p1 <- ggplot(cytotrace_pheno, aes(x=pheno, y=gcs,fill=pheno)) +
  geom_boxplot(color="black", notch = F) +theme(axis.text.x = element_text(angle = 90))
p1

p1 <- ggplot(cytotrace_pheno, aes(x=pheno, y=cytotrace,fill=pheno)) +
  geom_boxplot(color="black", notch = F) +theme(axis.text.x = element_text(angle = 90))
p1

save(p,p2,file = '../basic_analysis/cytotrace_SLICE/plot.RData')
###################################################### SLICE
#cpm <- apply(datause1,2, function(x) (x/sum(x))*10000) 
cpm<-as.data.frame(logcpm)
library(SLICE)
sc <- construct(exprmatrix=cpm, 
                cellidentity=pheno$pheno,
                projname="125")
## load the kappa similarity matrix
load("/media/ggj/home/ggj/Documents/data/fei/Rpackage/SLICE-master/SLICE.a12242016/data/mm_km.Rda")
#data(mm_kappasim)

rm(pbmc)
rm(datause1)
rm(cpm)
gc()

# bootstrap calculation of scEntropy
sc <- getEntropy(sc, km=km,                             # use the pre-computed kappa similarity matrix of mouse genes
                 calculation="bootstrap",               # choose the bootstrap calculation
                 B.num=25,                             # 100 iterations
                 exp.cutoff=1,                          # the threshold for expressed genes
                 B.size=1000,                           # the size of bootstrap sample
                 clustering.k=floor(sqrt(1000/2)),      # the number of functional clusters  
                 random.seed=123)                    # set the random seed to reproduce the results in the paper


write.csv(sc@entropies,file = "./intestine_entropies.csv",quote = F)

entropy <- sc@entropies
head(entropy)

dim(entropy)
dim(annouse)

slice_pheno <- cbind(pheno[rownames(entropy),],entropy)
head(slice_pheno)

slice_pheno$pheno <- factor(slice_pheno$pheno,levels = as.character(zz$pheno))
p1 <- ggplot(slice_pheno, aes(x=pheno, y=scEntropy.bootstrap,fill=pheno)) +
  geom_boxplot(color="black", notch = F) +theme(axis.text.x = element_text(angle = 90))
p1
write.csv(slice_pheno,file = 'slice_pheno.csv',quote = F)
magiccpm_pheno <- pheno
cor(magiccpm_pheno$order,magiccpm_pheno$PSI,method = "spearman")
cor(slice_pheno$order,slice_pheno$scEntropy.bootstrap,method = "spearman")
cor(cytotrace_pheno$order,cytotrace_pheno$gcs,method = "spearman")
cor(cytotrace_pheno$order,cytotrace_pheno$cytotrace,method = "spearman")
cor(logcpm_pheno$order,logcpm_pheno$PSI,method = "spearman")
save(magiccpm_pheno,slice_pheno,cytotrace_pheno,logcpm_pheno,file = "intestine_pheon.RData")



library(CytoTRACE)
setwd("/media/ggj/YF2/Qile/frog/")
library(Seurat)
library(reshape2)
pbmc <- readRDS("./allstage.pbmc.RDS")

### cytotrace

Idents(pbmc)<-pbmc$orig.ident
st48<-subset(pbmc,idents = "St48")
datause_st48<-st48@assays$RNA@counts
result <- CytoTRACE(datause_st48,ncores = 10,subsamplesize=10000)
save(result,file="cytotrace_entropy_St48.RData")

st54<-subset(pbmc,idents = "St54")
datause_st54<-st54@assays$RNA@counts
result <- CytoTRACE(datause_st54,ncores = 10)
save(result,file="cytotrace_entropy_St54.RData")


st59<-subset(pbmc,idents = "St59")
datause_st59<-st59@assays$RNA@counts
result <- CytoTRACE(datause_st59,ncores = 10)
save(result,file="cytotrace_entropy_St59.RData")


st66<-subset(pbmc,idents = "St66")
datause_st66<-st66@assays$RNA@counts
result <- CytoTRACE(datause_st66,ncores = 10)
save(result,file="cytotrace_entropy_St66.RData")
