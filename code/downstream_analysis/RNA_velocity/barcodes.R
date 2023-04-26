ob2.pitui=colnames(OB2_pituitary_Doubletfinder)
wt2.pitui=colnames(WT2_pituitary_Doubletfinder)
barcode.pitui=c(ob2.pitui,wt2.pitui)
barcode.pituitary=gsub( '^[^.]+','',barcode.pitui)
# barcode.pitui=gsub('[.]','',barcode.pituitary)
# barcode.pitui=gsub('_*','',barcode.pitui)
# barcode.pitui1=gsub('1','',barcode.pitui1)
write.table(barcode.pitui1,file = '/media/ggj/ggj/CJY/ob/velocity/barcode.pitui.tsv',quote = FALSE,row.names = F,col.names =F )

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


intersect(barcode.pitui1)
write.table(ob2.pitui,file = '/media/ggj/ggj/CJY/ob/velocity/barcode.tsv',quote = FALSE,row.names = F)

###########load OB3 WT3
setwd("/media/ggj/ggj/CJY/ob/doubletfinder/save3/")
library(reshape2)
file <- list.files(pattern = "^OB3|^WT3")
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

#############check cell barcode 491
setwd("/media/ggj/UU-4T-5/E150000491/L01/process")
dir <- list.files(".")
dir <- dir[c(49:60)]
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
    
    write.table(cellname1$cb,file = paste0('/media/ggj/ggj/CJY/ob/velocity/barcode/E150000491/col',dir[i],'_barcode.tsv'),quote = FALSE,row.names = F,col.names = F)
    
    setwd("../")
    
  }
  
  
}
#############check cell barcode 410
setwd("/media/ggj/UU-4T-3/E150000410/L01/process")
dir <- list.files(".")
dir <- dir[c(1:6,19:24)]
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
    
    write.table(cellname1$cb,file = paste0('/media/ggj/ggj/CJY/ob/velocity/barcode/E150000410/col',dir[i],'_barcode.tsv'),quote = FALSE,row.names = F,col.names = F)
    
    setwd("../")
    
  }
  
  
}
###########load WT6
setwd("/media/ggj/ggj/CJY/ob/doubletfinder/save3/")
library(reshape2)
file <- list.files(pattern = "^WT6")
file
cellname <- NULL
for(j in 1:length(file)){
  temp <- readRDS(file[j])
  cellname.temp <- colnames(temp)
  cellname.temp <- colsplit(cellname.temp,"\\.",names = c("index","cb"))
  cellname.temp$barcode <- colnames(temp)
  cellname.temp$cb <- gsub("_[0-9]*","",cellname.temp$cb)
  cellname <- rbind(cellname,cellname.temp)
  
}
dim(cellname)

cellname$index=gsub('COL','col',cellname$index)
#i=1
#head(cellname)

#############check cell barcode
setwd("/media/ggj/UU-4T-3/E150002225/L01/OUTPUT")
dir <- list.files(".")
dir <- dir[c(1,11:16,37:41)]
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
    
    write.table(cellname1$cb,file = paste0('/media/ggj/ggj/CJY/ob/velocity/barcode/E150002225/col',dir[i],'_barcode.tsv'),quote = FALSE,row.names = F,col.names = F)
    
    setwd("../")
    
  }
  
  
}
#############check cell barcode
setwd("/media/ggj/UU-4T-3/E150002200/L01/process/")
dir <- list.files(".")
dir <- dir[c(15:26)]
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
    
    write.table(cellname1$cb,file = paste0('/media/ggj/ggj/CJY/ob/velocity/barcode/E150002200/col',dir[i],'_barcode.tsv'),quote = FALSE,row.names = F,col.names = F)
    
    setwd("../")
    
  }
  
  
}
