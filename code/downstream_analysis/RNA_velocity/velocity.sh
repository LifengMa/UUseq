#!/bin/sh
cd /media/ggj/UU-4T-5/E150000430/L01/process/OB_WT_1_2/
for i in 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80; do
   cd $i
   str1=/media/ggj/ggj/CJY/ob/velocity/barcode/E150000430/
   str2=_barcode.tsv
   str3=/media/ggj/ggj/CJY/ob/velocity/output/E150000430/
   str4=col
    velocyto run -b $str1$str4$i$str2 -o output star_gene_exon_tagged.bam /media/ggj/ggj/CJY/ob/velocity/Mus_musculus.GRCm38.88.gtf
    mkdir $str3$i
    cd output
    mv *.loom $str3$i
    cd ../../
done


   
