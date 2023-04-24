setwd("/media/ggj/ggjlab2/UUseq/brain/benchmark/293T/genecoverage/")

library(reshape2)
library(ggplot2)
library(RColorBrewer)

file <- list.files()
file
file <- file[c(1:4,6,8,9)]

use <- NULL
for(i in 1:length(file)){
  data <- read.table(paste0(file[i],"/QC.geneBodyCoverage.byExpr.avgPct.txt.gz"),header =T)
  data <- data[,c(1,2)] 
  data$tech <- colsplit(file[i],"_",names = c("n1","n2"))$n2
  use <- rbind(use,data)
}

unique(use$tech)
use$tech <- gsub("VASA","VASA-plate",use$tech)
use$tech <- gsub("10X","10X Chromium",use$tech)
use$tech <- gsub("Smartseq3","Smart-seq3",use$tech)
use$tech <- gsub("split_3000cell","Split-seq",use$tech)
use$tech <- gsub("sci","Sci-seq",use$tech)
use <- use[-which(use$tech=="DroNc"),]
col_flg<-colorRampPalette(brewer.pal(8,"Set1"))(6)

ggplot(use,aes(x=QUANTILE,y=TOTAL,color=tech))+
  theme_classic()+geom_line()+
  ylab('Coverage')+xlab("Percentile of Gene Body (5'->3')")+
  scale_color_manual(values=col_flg)

ggsave("./result_0416.pdf",width = 8,height = 5)
