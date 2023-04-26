###########load OB12 WT12
setwd("/media/ggj/ggj/CJY/ob/doubletfinder/save3/")
library(reshape2)
file <- list.files(pattern = "^OB1_|^OB2|^WT1|^WT2")
file
cellname <- NULL
for(j in 1:length(file)){
  temp <- readRDS(file[j])
  cellname.temp <- colnames(temp)
  cellname.temp <- colsplit(cellname.temp,"\\.",names = c("index","cb"))
  cellname.temp$barcode <- colnames(temp)
  cellname.temp$cb <- gsub("_[0-9]","",cellname.temp$cb)
  cellname <- rbind(cellname,cellname.temp)
  
}
dim(cellname)
#i=1
#head(cellname)

#############check cell barcode
setwd("/media/ggj/UU-4T-5/E150000430/L01/process/OB_WT_1_2/")
dir <- list.files(".")
dir <- dir[c(13:28)]
dir
#i=1
for(i in 1:length(dir)){
  setwd(dir[i])
  summary.temp <- read.table("_dge.summary.txt",header = T)
  head(summary.temp)
  
  aa <- as.data.frame(table(cellname$index==paste0("col",dir[i])))
  head(aa)
  if(length(rownames(aa))==1){
    setwd("../")
    next
  }
  else{
    cellname1 <- cellname[cellname$index==paste0("col",dir[i]),]
    if(length(intersect(cellname1$cb,summary.temp$CELL_BARCODE))==length(cellname1$cb)){
      message("col ",dir[i], " is OK!")
    }
    else{
      message("col ",dir[i], " is not OK!")
    }
    
    write.table(cellname1$cb,file = paste0('/media/ggj/ggj/CJY/ob/velocity/col',dir[i],'_barcode.tsv'),quote = FALSE,row.names = F,col.names = F)
    
    setwd("../")
    
  }
  
  
  }
