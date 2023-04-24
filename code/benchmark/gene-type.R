setwd("/media/ggj/ggjlab2/UUseq/brain/benchmark/293T/")
library(data.table)
#gene-type
rm(list = ls())
gc()
##get gtf

dat<-fread('./gene_anno/Homo_sapiens.GRCh38.gtf',header = F)
head(dat)
dat <- dat[dat$V3=='gene',]
dat <- as.data.frame(dat)
dat$biotype <- colsplit(dat$V9,";",names = c("n1","n2","n3","n4","n5","n6","n7"))$n5
dat$biotype <- gsub("gene_biotype","",dat$biotype)
dat$biotype <- gsub("\"","",dat$biotype)
dat$biotype <- gsub(" ","",dat$biotype)

dat$gene <- colsplit(dat$V9,";",names = c("n1","n2","n3","n4","n5","n6","n7"))$n3
dat$gene <- gsub("gene_name|\"| ","",dat$gene)

use <- dat[,c(11,10)]
saveRDS(use,"./gene_anno/GRCh38_gene_biotype.rds")



rm(list = ls())

############
setwd('/media/ggj/ggjlab2/UUseq/brain/benchmark/293T/gene_anno/data/')
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(dplyr)

use <- readRDS("../GRCh38_gene_biotype.rds")


###modify biotype
biotype <-  as.data.frame(unique(use$biotype))
othertype <- setdiff(unique(use$biotype),c("protein_coding","lincRNA",
                                           "snRNA","misc_RNA" ,"snoRNA","scaRNA","miRNA","Mt_tRNA", "ribozyme",####small non coding RNA
                                           "pseudogene","unprocessed_pseudogene","processed_pseudogene",
                                           "transcribed_processed_pseudogene","unitary_pseudogene","polymorphic_pseudogene",
                                           "transcribed_unitary_pseudogene","IG_V_pseudogene" ,"TR_V_pseudogene","IG_C_pseudogene",
                                           "TR_J_pseudogene","IG_J_pseudogene","IG_pseudogene","transcribed_unprocessed_pseudogene",
                                           "rRNA","Mt_rRNA"))
use$biotype2 <- ifelse(use$biotype%in%othertype,"Other",use$biotype)
use[which(use$gene %like% "^MT-" & use$biotype2=="protein_coding") ,]$biotype2 <- 'MT genes'
use$biotype2 <- gsub("snRNA|misc_RNA|snoRNA|miRNA|ribozyme|Mt_tRNA|scaRNA","sncRNA",use$biotype2)
use$biotype2 <- gsub("^pseudogene|unprocessed_pseudogene|processed_pseudogene|transcribed_processed_pseudogene|unitary_pseudogene|polymorphic_pseudogene|transcribed_unitary_pseudogene|IG_V_pseudogene|TR_V_pseudogene|IG_C_pseudogene|TR_J_pseudogene|IG_J_pseudogene|IG_pseudogene|transcribed_unprocessed_pseudogene" ,"Pseudogenes",use$biotype2)
table(use$biotype2)
use$biotype3 <- ifelse(use$biotype%in%othertype,"Other",use$biotype)

#####
UUseq = read.table('./type_UUseq.txt.gz', header = F)
#x10 = read.table('./type_10X.txt.gz', header = F)
#vasa = read.table('./type_VASA.txt.gz', header = F)
#smart = read.table('./type_smartseq3.txt.gz', header = F)
MATQ = read.table('./type_MATQ.txt.gz', header = F)
sci= read.table('./type_sci.txt.gz', header = F)
split= read.table('./type_split_3000cell.txt.gz', header = F)

UUseq$V1 = gsub('gn:Z:', '', UUseq$V1) ;
#x10$V1 = gsub('gn:Z:', '', x10$V1); 
#vasa$V1 = gsub('gn:Z:', '', vasa$V1); 
#smart$V1 = gsub('gn:Z:', '', smart$V1);
MATQ$V1 = gsub('gn:Z:', '', MATQ$V1) 
sci$V1 = gsub('gn:Z:', '', sci$V1) 
split$V1 = gsub('gn:Z:', '', split$V1) 

# length(grep(',', UUseq$V1)) #35547
# length(grep(',', x10$V1)) #113530
# length(grep(',', vasa$V1)) #33569

t = as.data.frame(table(use$biotype)); colnames(t) = c('biotype', 'num')
#write.csv(t, file = './output/mix_gtf_biotype.csv')

#colnames(UUseq) = colnames(x10) = colnames(vasa) = colnames(smart) =colnames(MATQ)= 'gene'
colnames(UUseq) = colnames(sci) = colnames(split) = colnames(MATQ)= 'gene'

UUseq.full = left_join(UUseq, use, by = 'gene')
UUseq.full[grep(',', UUseq.full$gene),]$biotype2 = 'Multi'
UUseq.full[grep(',', UUseq.full$gene),]$biotype = 'Multi'

# x10.full = left_join(x10, use, by = 'gene')
# x10.full[grep(',', x10.full$gene),]$biotype2 = 'Multi'
# x10.full[grep(',', x10.full$gene),]$biotype = 'Multi'


# vasa.full = left_join(vasa, use, by = 'gene')
# vasa.full[grep(',', vasa.full$gene),]$biotype2 = 'Multi'
# vasa.full[grep(',', vasa.full$gene),]$biotype = 'Multi'
# 
# smart.full = left_join(smart, use, by = 'gene')
# smart.full[grep(',', smart.full$gene),]$biotype2 = 'Multi'
# smart.full[grep(',', smart.full$gene),]$biotype = 'Multi'

MATQ.full = left_join(MATQ, use, by = 'gene')
MATQ.full[grep(',', MATQ.full$gene),]$biotype2 = 'Multi'
MATQ.full[grep(',', MATQ.full$gene),]$biotype = 'Multi'

sci.full = left_join(sci, use, by = 'gene')
sci.full[grep(',', sci.full$gene),]$biotype2 = 'Multi'
sci.full[grep(',', sci.full$gene),]$biotype = 'Multi'

split.full = left_join(split, use, by = 'gene')
split.full[grep(',', split.full$gene),]$biotype2 = 'Multi'
split.full[grep(',', split.full$gene),]$biotype = 'Multi'

UUseq.final = UUseq.full[-grep('ERCC', UUseq.full$gene),]
# x10.final = x10.full[-grep('ERCC', x10.full$gene),]
# vasa.final = vasa.full[-grep('ERCC', vasa.full$gene),]
# smart.final = smart.full[-grep('ERCC', smart.full$gene),]
MATQ.final = MATQ.full[-grep('ERCC', MATQ.full$gene),]
sci.final = sci.full[-grep('ERCC', sci.full$gene),]
split.final = split.full[-grep('ERCC', split.full$gene),]

## p1
## Protein coding, lincRNA, miscRNA, Pseudogenes, rRNA & MT-rRNA, Antisense, Multi Annotated, Other
## protein_coding, lincRNA, misc_RNA, pseudogene, rRNA & Mt_rRNA, antisense, Multi, 

data.matrix = data.frame()
# target = c('protein_coding', 'lincRNA', 'misc_RNA', 'pseudogene', 'rRNA', 'Mt_rRNA', 'antisense')
# rest = setdiff(unique(use$biotype), target)
# rest = na.omit(rest)

# data.matrix = rbind(c(sum(UUseq.final$biotype2 == 'protein_coding'), sum(UUseq.final$biotype2 == 'lincRNA'), sum(UUseq.final$biotype2 == 'sncRNA'), sum(UUseq.final$biotype2 == 'Pseudogenes'), sum(UUseq.final$biotype2 %in% c('rRNA', 'Mt_rRNA')), sum(UUseq.final$biotype2 == 'MT genes'), sum(UUseq.final$biotype2 == 'Multi'), sum(UUseq.final$biotype2 %in% 'Other')),
#                     c(sum(vasa.final$biotype2 == 'protein_coding'), sum(vasa.final$biotype2 == 'lincRNA'), sum(vasa.final$biotype2 == 'sncRNA'), sum(vasa.final$biotype2 == 'Pseudogenes'), sum(vasa.final$biotype2 %in% c('rRNA', 'Mt_rRNA')), sum(vasa.final$biotype2 == 'MT genes'), sum(vasa.final$biotype2 == 'Multi'), sum(vasa.final$biotype2 %in% 'Other')),
#                     c(sum(x10.final$biotype2 == 'protein_coding'), sum(x10.final$biotype2 == 'lincRNA'), sum(x10.final$biotype2 == 'sncRNA'), sum(x10.final$biotype2 == 'Pseudogenes'), sum(x10.final$biotype2 %in% c('rRNA', 'Mt_rRNA')), sum(x10.final$biotype2 == 'MT genes'), sum(x10.final$biotype2 == 'Multi'), sum(x10.final$biotype2 %in% 'Other')),
#                     c(sum(smart.final$biotype2 == 'protein_coding'), sum(smart.final$biotype2 == 'lincRNA'), sum(smart.final$biotype2 == 'sncRNA'), sum(smart.final$biotype2 == 'Pseudogenes'), sum(smart.final$biotype2 %in% c('rRNA', 'Mt_rRNA')), sum(smart.final$biotype2 == 'MT genes'), sum(smart.final$biotype2 == 'Multi'), sum(smart.final$biotype2 %in% 'Other')),
#                     c(sum(MATQ.final$biotype2 == 'protein_coding'), sum(MATQ.final$biotype2 == 'lincRNA'), sum(MATQ.final$biotype2 == 'sncRNA'), sum(MATQ.final$biotype2 == 'Pseudogenes'), sum(MATQ.final$biotype2 %in% c('rRNA', 'Mt_rRNA')), sum(MATQ.final$biotype2 == 'MT genes'), sum(MATQ.final$biotype2 == 'Multi'), sum(MATQ.final$biotype2 %in% 'Other'))
# )

data.matrix = rbind(c(sum(UUseq.final$biotype2 == 'protein_coding'), sum(UUseq.final$biotype2 == 'lincRNA'), sum(UUseq.final$biotype2 == 'sncRNA'), sum(UUseq.final$biotype2 == 'Pseudogenes'), sum(UUseq.final$biotype2 %in% c('rRNA', 'Mt_rRNA')), sum(UUseq.final$biotype2 == 'MT genes'), sum(UUseq.final$biotype2 == 'Multi'), sum(UUseq.final$biotype2 %in% 'Other')),
                    #c(sum(vasa.final$biotype2 == 'protein_coding'), sum(vasa.final$biotype2 == 'lincRNA'), sum(vasa.final$biotype2 == 'sncRNA'), sum(vasa.final$biotype2 == 'Pseudogenes'), sum(vasa.final$biotype2 %in% c('rRNA', 'Mt_rRNA')), sum(vasa.final$biotype2 == 'MT genes'), sum(vasa.final$biotype2 == 'Multi'), sum(vasa.final$biotype2 %in% 'Other')),
                    #c(sum(x10.final$biotype2 == 'protein_coding'), sum(x10.final$biotype2 == 'lincRNA'), sum(x10.final$biotype2 == 'sncRNA'), sum(x10.final$biotype2 == 'Pseudogenes'), sum(x10.final$biotype2 %in% c('rRNA', 'Mt_rRNA')), sum(x10.final$biotype2 == 'MT genes'), sum(x10.final$biotype2 == 'Multi'), sum(x10.final$biotype2 %in% 'Other')),
                    #c(sum(smart.final$biotype2 == 'protein_coding'), sum(smart.final$biotype2 == 'lincRNA'), sum(smart.final$biotype2 == 'sncRNA'), sum(smart.final$biotype2 == 'Pseudogenes'), sum(smart.final$biotype2 %in% c('rRNA', 'Mt_rRNA')), sum(smart.final$biotype2 == 'MT genes'), sum(smart.final$biotype2 == 'Multi'), sum(smart.final$biotype2 %in% 'Other')),
                    c(sum(MATQ.final$biotype2 == 'protein_coding'), sum(MATQ.final$biotype2 == 'lincRNA'), sum(MATQ.final$biotype2 == 'sncRNA'), sum(MATQ.final$biotype2 == 'Pseudogenes'), sum(MATQ.final$biotype2 %in% c('rRNA', 'Mt_rRNA')), sum(MATQ.final$biotype2 == 'MT genes'), sum(MATQ.final$biotype2 == 'Multi'), sum(MATQ.final$biotype2 %in% 'Other')),
                    c(sum(split.final$biotype2 == 'protein_coding'), sum(split.final$biotype2 == 'lincRNA'), sum(split.final$biotype2 == 'sncRNA'), sum(split.final$biotype2 == 'Pseudogenes'), sum(split.final$biotype2 %in% c('rRNA', 'Mt_rRNA')), sum(split.final$biotype2 == 'MT genes'), sum(split.final$biotype2 == 'Multi'), sum(split.final$biotype2 %in% 'Other')),
                    c(sum(sci.final$biotype2 == 'protein_coding'), sum(sci.final$biotype2 == 'lincRNA'), sum(sci.final$biotype2 == 'sncRNA'), sum(sci.final$biotype2 == 'Pseudogenes'), sum(sci.final$biotype2 %in% c('rRNA', 'Mt_rRNA')), sum(sci.final$biotype2 == 'MT genes'), sum(sci.final$biotype2 == 'Multi'), sum(sci.final$biotype2 %in% 'Other'))
)


data.matrix = as.data.frame(data.matrix)

colnames(data.matrix) = c('Protein coding', 'lincRNA', 'sncRNA', 'Pseudogenes', 'rRNA & MT-rRNA', 'MT genes','Multi Annotated', 'Other')
data.matrix = data.matrix / rowSums(data.matrix)
#data.matrix$expr = c('UU-seq', 'VASA-plate', '10X Chromium', 'Smart-seq3', 'MATQ')
data.matrix$expr = c('UU-seq', 'MATQ', 'Split-seq', 'Sci-seq')
data.melt = reshape2::melt(data.matrix)


library(RColorBrewer)
colors <-colorRampPalette(brewer.pal(8,"Dark2"))(8)

#data.melt$expr = factor(data.melt$expr, levels = rev(c('UU-seq', 'VASA-plate', 'Smart-seq3', 'MATQ','10X Chromium')))
data.melt$expr = factor(data.melt$expr, levels = rev(c('UU-seq', 'MATQ', 'Split-seq', 'Sci-seq')))

p1 = ggplot(data.melt, aes(x = expr, y = value * 100, fill = variable)) + geom_bar(stat = 'identity') +
  scale_fill_manual(values = colors) + theme_classic() + xlab('') + 
  ylab('% of annotated reads') + coord_flip() + 
  theme(legend.position = 'top', 
        legend.title = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x = element_text(colour = 'black', size = 12),
        axis.text.y = element_text(color = 'black', size = 12))
p1
ggsave(p1, filename = '../result_gene_anno_all_0416.pdf', width = 9, height = 6)


## p2
# snoRNA snRNA rRNA MiRNA MiscRNA  MT-tRNA这几个
# snoRNA, snRNA, miRNA, misc_RNA, Mt_tRNA
#"snRNA","misc_RNA" ,"snoRNA","scaRNA","miRNA","Mt_tRNA", "ribozyme"
p2.data = data.frame()
# p2.data = rbind(c(sum(UUseq.final$biotype == 'snoRNA'),sum(UUseq.final$biotype == 'scaRNA'), sum(UUseq.final$biotype == 'snRNA'), sum(UUseq.final$biotype == 'miRNA'), sum(UUseq.final$biotype %in% c('misc_RNA')), sum(UUseq.final$biotype == 'Mt_tRNA'),sum(UUseq.final$biotype == "ribozyme")),
#                 c(sum(vasa.final$biotype == 'snoRNA'),sum(vasa.final$biotype == 'scaRNA'), sum(vasa.final$biotype == 'snRNA'), sum(vasa.final$biotype == 'miRNA'), sum(vasa.final$biotype %in% c('misc_RNA')), sum(vasa.final$biotype == 'Mt_tRNA'),sum(vasa.final$biotype == "ribozyme")),
#                 c(sum(x10.final$biotype == 'snoRNA'),sum(x10.final$biotype == 'scaRNA'), sum(x10.final$biotype == 'snRNA'), sum(x10.final$biotype == 'miRNA'), sum(x10.final$biotype %in% c('misc_RNA')), sum(x10.final$biotype == 'Mt_tRNA'),sum(x10.final$biotype == "ribozyme")),
#                 c(sum(smart.final$biotype == 'snoRNA'),sum(smart.final$biotype == 'scaRNA'), sum(smart.final$biotype == 'snRNA'), sum(smart.final$biotype == 'miRNA'), sum(smart.final$biotype %in% c('misc_RNA')), sum(smart.final$biotype == 'Mt_tRNA'),sum(smart.final$biotype == "ribozyme")),
#                 c(sum(MATQ.final$biotype == 'snoRNA'),sum(MATQ.final$biotype == 'scaRNA'), sum(MATQ.final$biotype == 'snRNA'), sum(MATQ.final$biotype == 'miRNA'), sum(MATQ.final$biotype %in% c('misc_RNA')), sum(MATQ.final$biotype == 'Mt_tRNA'),sum(MATQ.final$biotype == "ribozyme"))
# )
p2.data = rbind(c(sum(UUseq.final$biotype == 'snoRNA'),sum(UUseq.final$biotype == 'scaRNA'), sum(UUseq.final$biotype == 'snRNA'), sum(UUseq.final$biotype == 'miRNA'), sum(UUseq.final$biotype %in% c('misc_RNA')), sum(UUseq.final$biotype == 'Mt_tRNA'),sum(UUseq.final$biotype == "ribozyme")),
                #c(sum(vasa.final$biotype == 'snoRNA'),sum(vasa.final$biotype == 'scaRNA'), sum(vasa.final$biotype == 'snRNA'), sum(vasa.final$biotype == 'miRNA'), sum(vasa.final$biotype %in% c('misc_RNA')), sum(vasa.final$biotype == 'Mt_tRNA'),sum(vasa.final$biotype == "ribozyme")),
                #c(sum(x10.final$biotype == 'snoRNA'),sum(x10.final$biotype == 'scaRNA'), sum(x10.final$biotype == 'snRNA'), sum(x10.final$biotype == 'miRNA'), sum(x10.final$biotype %in% c('misc_RNA')), sum(x10.final$biotype == 'Mt_tRNA'),sum(x10.final$biotype == "ribozyme")),
                #c(sum(smart.final$biotype == 'snoRNA'),sum(smart.final$biotype == 'scaRNA'), sum(smart.final$biotype == 'snRNA'), sum(smart.final$biotype == 'miRNA'), sum(smart.final$biotype %in% c('misc_RNA')), sum(smart.final$biotype == 'Mt_tRNA'),sum(smart.final$biotype == "ribozyme")),
                c(sum(MATQ.final$biotype == 'snoRNA'),sum(MATQ.final$biotype == 'scaRNA'), sum(MATQ.final$biotype == 'snRNA'), sum(MATQ.final$biotype == 'miRNA'), sum(MATQ.final$biotype %in% c('misc_RNA')), sum(MATQ.final$biotype == 'Mt_tRNA'),sum(MATQ.final$biotype == "ribozyme")),
                c(sum(split.final$biotype == 'snoRNA'),sum(split.final$biotype == 'scaRNA'), sum(split.final$biotype == 'snRNA'), sum(split.final$biotype == 'miRNA'), sum(split.final$biotype %in% c('misc_RNA')), sum(split.final$biotype == 'Mt_tRNA'),sum(split.final$biotype == "ribozyme")),
                c(sum(sci.final$biotype == 'snoRNA'),sum(sci.final$biotype == 'scaRNA'), sum(sci.final$biotype == 'snRNA'), sum(sci.final$biotype == 'miRNA'), sum(sci.final$biotype %in% c('misc_RNA')), sum(sci.final$biotype == 'Mt_tRNA'),sum(sci.final$biotype == "ribozyme"))
)

p2.data = as.data.frame(p2.data)

colnames(p2.data) = c('snoRNA', 'scaRNA','snRNA', 'miRNA', 'miscRNA', 'MT-tRNA','ribozyme')


#p2.data = p2.data / rowSums(p2.data)
p2.data[1,] = p2.data[1,] / dim(UUseq.final)[1]
#p2.data[2,] = p2.data[2,] / dim(vasa.final)[1]
#p2.data[3,] = p2.data[3,] / dim(x10.final)[1]
#p2.data[4,] = p2.data[4,] / dim(smart.final)[1]
p2.data[2,] = p2.data[2,] / dim(MATQ.final)[1]
p2.data[3,] = p2.data[3,] / dim(split.final)[1]
p2.data[4,] = p2.data[4,] / dim(sci.final)[1]

#p2.data$expr = c('UU-seq', 'VASA-plate', '10X Chromium', 'Smart-seq3', 'MATQ')
p2.data$expr = c('UU-seq',  'MATQ','Split-seq','Sci-seq')

p2.melt = reshape2::melt(p2.data)

#p2.melt$expr = factor(p2.melt$expr, levels = rev(c('UU-seq', 'VASA-plate', 'Smart-seq3', 'MATQ','10X Chromium')))
p2.melt$expr = factor(p2.melt$expr, levels = rev(c('UU-seq', 'MATQ','Split-seq','Sci-seq')))

colors2 <- colorRampPalette(brewer.pal(7,"Accent"))(7)
p2 = ggplot(p2.melt, aes(x = expr, y = value * 100, fill = variable)) + geom_bar(stat = 'identity') + 
  scale_fill_manual(values = colors2) + theme_classic() + xlab('') + 
  ylab('% of annotated reads') + coord_flip() + 
  theme(legend.position = 'top', 
        legend.title = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x = element_text(colour = 'black', size = 12),
        axis.text.y = element_text(color = 'black', size = 12))
p2
ggsave(p2, filename = '../result_gene_anno_sncRNA_0416.pdf', width = 9, height = 6)


