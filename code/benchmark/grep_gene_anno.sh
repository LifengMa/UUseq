#BSUB -q normal
#BSUB -J gene_anno
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -R "span[ptile=1]"
#BSUB -n 1

sample_name=$(basename `pwd`)
dropseq_root=/share/home/hanxiaoping/tools/Drop-seq_tools-2.5.1
other_root=/share/home/hanxiaoping/RAWDATA/OB/benchmark/genebody_coverage
MATQ_root=/share/home/hanxiaoping/RAWDATA/OB/benchmark/MATQ/GRCh38
Smartseq3_root=/share/home/hanxiaoping/RAWDATA/OB/benchmark/smartseq3/correct/


#samtools view $other_root/VASA_downsample.bam |awk '{for (i=1; i<=NF; ++i) {if($i ~ "^gn:Z:"){print $i}}}' > type_VASA.txt
#samtools view $other_root/UUseq_downsample.bam |awk '{for (i=1; i<=NF; ++i) {if($i ~ "^gn:Z:"){print $i}}}' > type_UUseq.txt
#samtools view $other_root/10X_downsample.bam |awk '{for (i=1; i<=NF; ++i) {if($i ~ "^gn:Z:"){print $i}}}' > type_10X.txt
#samtools view $MATQ_root/star_gene_exon_tagged.bam |awk '{for (i=1; i<=NF; ++i) {if($i ~ "^gn:Z:"){print $i}}}' > type_MATQ.txt

#samtools view -@ 16 -H $Smartseq3_root/star_gene_exon_tagged_corrected.bam > SAM_header
#samtools view $Smartseq3_root/star_gene_exon_tagged_corrected.bam | shuf -n 9000000 > filtered_SAM_body
#cat SAM_header filtered_SAM_body > Smartseq3_downsample.sam
#samtools view -@ 16 -b Smartseq3_downsample.sam > Smartseq3_downsample.bam &&rm filtered_SAM_body Smartseq3_downsample.sam SAM_header
#samtools view Smartseq3_downsample.bam |awk '{for (i=1; i<=NF; ++i) {if($i ~ "^gn:Z:"){print $i}}}' > type_smartseq3.txt

samtools view $other_root/split_300cell_downsample.bam |awk '{for (i=1; i<=NF; ++i) {if($i ~ "^gn:Z:"){print $i}}}' > type_split_300cell.txt
samtools view $other_root/split_3000cell_downsample.bam |awk '{for (i=1; i<=NF; ++i) {if($i ~ "^gn:Z:"){print $i}}}' > type_split_3000cell.txt
samtools view $other_root/DroNc_downsample.bam |awk '{for (i=1; i<=NF; ++i) {if($i ~ "^gn:Z:"){print $i}}}' > type_DroNc.txt
samtools view $other_root/sci_downsample.bam |awk '{for (i=1; i<=NF; ++i) {if($i ~ "^gn:Z:"){print $i}}}' > type_sci.txt

gzip *.txt