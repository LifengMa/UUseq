rm(list = ls())

setwd('~/Desktop/MW3.0')
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(dplyr)

load('./biotype_anno/mix_gene_biotype.rdata')
mw3 = read.table('./biotype_anno/mw3-all_type.txt', header = F)
x10 = read.table('./biotype_anno/10x-all_type.txt', header = F)
# one sample.
#vasa = read.table('./biotype_anno/vasa_type.txt', header = F)
# all sample combined.
vasa = read.table('./biotype_anno/vasa_type_new.txt', header = F)
smart = read.table('./biotype_anno/smart-single_type.txt', header = F)

mw3$V1 = gsub('gn:Z:', '', mw3$V1) ; x10$V1 = gsub('gn:Z:', '', x10$V1); vasa$V1 = gsub('gn:Z:', '', vasa$V1); smart$V1 = gsub('gn:Z:', '', smart$V1)

# length(grep(',', mw3$V1)) #35547
# length(grep(',', x10$V1)) #113530
# length(grep(',', vasa$V1)) #33569

t = as.data.frame(table(use$biotype)); colnames(t) = c('biotype', 'num')
write.csv(t, file = './output/mix_gtf_biotype.csv')

colnames(mw3) = colnames(x10) = colnames(vasa) = colnames(smart) = 'gene'

mw3.full = left_join(mw3, use, by = 'gene')
mw3.full[grep(',', mw3.full$gene),]$biotype = 'Multi'

x10.full = left_join(x10, use, by = 'gene')
x10.full[grep(',', x10.full$gene),]$biotype = 'Multi'

vasa.full = left_join(vasa, use, by = 'gene')
vasa.full[grep(',', vasa.full$gene),]$biotype = 'Multi'

smart.full = left_join(smart, use, by = 'gene')
smart.full[grep(',', smart.full$gene),]$biotype = 'Multi'


mw3.final = mw3.full[-grep('ERCC', mw3.full$gene),]
x10.final = x10.full[-grep('ERCC', x10.full$gene),]
vasa.final = vasa.full[-grep('ERCC', vasa.full$gene),]
smart.final = smart.full[-grep('ERCC', smart.full$gene),]

## p1
## Protein coding, lincRNA, miscRNA, Pseudogenes, rRNA & MT-rRNA, Antisense, Multi Annotated, Other
## protein_coding, lincRNA, misc_RNA, pseudogene, rRNA & Mt_rRNA, antisense, Multi, 

data.matrix = data.frame()
target = c('protein_coding', 'lincRNA', 'misc_RNA', 'pseudogene', 'rRNA', 'Mt_rRNA', 'antisense')
rest = setdiff(unique(use$biotype), target)
rest = na.omit(rest)

data.matrix = rbind(c(sum(mw3.final$biotype == 'protein_coding'), sum(mw3.final$biotype == 'lincRNA'), sum(mw3.final$biotype == 'misc_RNA'), sum(mw3.final$biotype == 'pseudogene'), sum(mw3.final$biotype %in% c('rRNA', 'Mt_rRNA')), sum(mw3.final$biotype == 'antisense'), sum(mw3.final$biotype == 'Multi'), sum(mw3.final$biotype %in% rest)),
                    c(sum(vasa.final$biotype == 'protein_coding'), sum(vasa.final$biotype == 'lincRNA'), sum(vasa.final$biotype == 'misc_RNA'), sum(vasa.final$biotype == 'pseudogene'), sum(vasa.final$biotype %in% c('rRNA', 'Mt_rRNA')), sum(vasa.final$biotype == 'antisense'), sum(vasa.final$biotype == 'Multi'), sum(vasa.final$biotype %in% rest)),
                    c(sum(x10.final$biotype == 'protein_coding'), sum(x10.final$biotype == 'lincRNA'), sum(x10.final$biotype == 'misc_RNA'), sum(x10.final$biotype == 'pseudogene'), sum(x10.final$biotype %in% c('rRNA', 'Mt_rRNA')), sum(x10.final$biotype == 'antisense'), sum(x10.final$biotype == 'Multi'), sum(x10.final$biotype %in% rest)),
                    c(sum(smart.final$biotype == 'protein_coding'), sum(smart.final$biotype == 'lincRNA'), sum(smart.final$biotype == 'misc_RNA'), sum(smart.final$biotype == 'pseudogene'), sum(smart.final$biotype %in% c('rRNA', 'Mt_rRNA')), sum(smart.final$biotype == 'antisense'), sum(smart.final$biotype == 'Multi'), sum(smart.final$biotype %in% rest)))


data.matrix = as.data.frame(data.matrix)

colnames(data.matrix) = c('Protein coding', 'lincRNA', 'miscRNA', 'Pseudogenes', 'rRNA & MT-rRNA', 'Antisense','Multi Annotated', 'Other')
data.matrix = data.matrix / rowSums(data.matrix)
data.matrix$expr = c('Microwell-seq 3.0', 'VASA-plate', '10X Genomics', 'SMART-seq Total')

data.melt = reshape2::melt(data.matrix)


colors = brewer.pal(8, 'Spectral')

data.melt$expr = factor(data.melt$expr, levels = c('SMART-seq Total','10X Genomics', 'VASA-plate', 'Microwell-seq 3.0'))

p1 = ggplot(data.melt, aes(x = expr, y = value * 100, fill = variable)) + geom_bar(stat = 'identity') + scale_fill_manual(values = colors) + theme_classic() + xlab('') + ylab('% of annotated reads') + coord_flip() + 
  theme(legend.position = 'top', 
        legend.title = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x = element_text(colour = 'black', size = 9),
        axis.text.y = element_text(color = 'black', size = 9))
ggsave(p1, filename = '~/Desktop/MW3.0/output/RNA/reads_anno_p1_full.pdf', width = 9, height = 6)


## p2
# snoRNA snRNA rRNA MiRNA MiscRNA  MT-tRNA这几个
# snoRNA, snRNA, miRNA, misc_RNA, Mt_tRNA

p2.data = data.frame()
p2.data = rbind(c(sum(mw3.final$biotype == 'snoRNA'), sum(mw3.final$biotype == 'snRNA'), sum(mw3.final$biotype == 'miRNA'), sum(mw3.final$biotype %in% c('misc_RNA')), sum(mw3.final$biotype == 'Mt_tRNA')),
                    c(sum(vasa.final$biotype == 'snoRNA'), sum(vasa.final$biotype == 'snRNA'), sum(vasa.final$biotype == 'miRNA'), sum(vasa.final$biotype %in% c('misc_RNA')), sum(vasa.final$biotype == 'Mt_tRNA')),
                    c(sum(x10.final$biotype == 'snoRNA'), sum(x10.final$biotype == 'snRNA'), sum(x10.final$biotype == 'miRNA'), sum(x10.final$biotype %in% c('misc_RNA')), sum(x10.final$biotype == 'Mt_tRNA')),
                c(sum(smart.final$biotype == 'snoRNA'), sum(smart.final$biotype == 'snRNA'), sum(smart.final$biotype == 'miRNA'), sum(smart.final$biotype %in% c('misc_RNA')), sum(smart.final$biotype == 'Mt_tRNA')))

p2.data = as.data.frame(p2.data)

colnames(p2.data) = c('snoRNA', 'snRNA', 'miRNA', 'miscRNA', 'MT-tRNA')


#p2.data = p2.data / rowSums(p2.data)
p2.data[1,] = p2.data[1,] / dim(mw3.final)[1]
p2.data[2,] = p2.data[2,] / dim(vasa.final)[1]
p2.data[3,] = p2.data[3,] / dim(x10.final)[1]
p2.data[4,] = p2.data[4,] / dim(smart.final)[1]

p2.data$expr = c('Microwell-seq 3.0', 'VASA-plate', '10X Genomics', 'SMART-seq Total')

p2.melt = reshape2::melt(p2.data)

p2.melt$expr = factor(p2.melt$expr, levels = c('SMART-seq Total', '10X Genomics', 'VASA-plate', 'Microwell-seq 3.0'))

p2 = ggplot(p2.melt, aes(x = expr, y = value * 100, fill = variable)) + geom_bar(stat = 'identity') + scale_fill_manual(values = colors) + theme_classic() + xlab('') + ylab('% of annotated reads') + coord_flip() + 
  theme(legend.position = 'top', 
        legend.title = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x = element_text(colour = 'black', size = 9),
        axis.text.y = element_text(color = 'black', size = 9))
ggsave(p2, filename = '~/Desktop/MW3.0/output/RNA/reads_anno_p2_full.pdf', width = 9, height = 6)

