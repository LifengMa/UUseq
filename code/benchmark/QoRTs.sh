#BSUB -q normal
#BSUB -J QoRT
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -R "span[ptile=18]"
#BSUB -n 18

cd pwd

QoRT_root=/share/home/hanxiaoping/tools/QoRTs-master/
dropseq_root=/share/home/hanxiaoping/tools/Drop-seq_tools-2.5.1
VASA_root=/share/home/hanxiaoping/RAWDATA/OB/benchmark/VASA/
UUseq_root=/share/home/hanxiaoping/RAWDATA/OB/benchmark/UUseq/
TenX_root=/share/home/hanxiaoping/RAWDATA/OB/benchmark/10X/GRch38/
Smartseq3_root=/share/home/hanxiaoping/RAWDATA/OB/benchmark/smartseq3/all/
MATQ_root=/share/home/hanxiaoping/RAWDATA/OB/benchmark/MATQ/GRCh38
sci_root=/share/home/hanxiaoping/RAWDATA/OB/benchmark/sci-seq/GRCh38/
split_root=/share/home/hanxiaoping/RAWDATA/OB/benchmark/split-seq/rawdata/300cell/GRCh38/
split_root2=/share/home/hanxiaoping/RAWDATA/OB/benchmark/split-seq/rawdata/3000cell/GRCh38/
DroNc_root=/share/home/hanxiaoping/RAWDATA/OB/benchmark/DroNc-seq/


###MATQ
#java -jar $QoRT_root/QoRTs.jar QC --singleEnded --numThreads 8 $MATQ_root/star_gene_exon_tagged.bam /share/home/hanxiaoping/tools/STAR_Reference_Human/Homo_sapiens.GRCh38.fa.gtf \
# output_MATQ 
 
 
###VASA
#samtools view -@ 16 -H $VASA_root/star_gene_exon_tagged.bam > SAM_header
#samtools view -@ 16 $VASA_root/star_gene_exon_tagged.bam | shuf -n 9000000 > filtered_SAM_body
#cat SAM_header filtered_SAM_body > VASA_downsample.sam
#samtools view -@ 16 -b VASA_downsample.sam > VASA_downsample.bam &&rm filtered_SAM_body VASA_downsample.sam SAM_header

#java -jar $QoRT_root/QoRTs.jar QC --singleEnded --numThreads 12 VASA_downsample.bam /share/home/hanxiaoping/tools/STAR_Reference_Human/Homo_sapiens.GRCh38.fa.gtf output_VASA

###UUseq
#samtools view -@ 16 -H $UUseq_root/star_gene_exon_tagged.bam > SAM_header
#samtools view $UUseq_root/star_gene_exon_tagged.bam | shuf -n 9000000 > filtered_SAM_body
#cat SAM_header filtered_SAM_body > UUseq_downsample.sam
#samtools view -@ 16 -b UUseq_downsample.sam > UUseq_downsample.bam &&rm filtered_SAM_body UUseq_downsample.sam SAM_header

#java -jar $QoRT_root/QoRTs.jar QC --singleEnded --numThreads 12 UUseq_downsample.bam /share/home/hanxiaoping/tools/STAR_Reference_Human/Homo_sapiens.GRCh38.fa.gtf output_UUseq

###TenX
#samtools view -@ 16 -H $TenX_root/star_gene_exon_tagged.bam > SAM_header
#samtools view $TenX_root/star_gene_exon_tagged.bam | shuf -n 9000000 > filtered_SAM_body
#cat SAM_header filtered_SAM_body > 10X_downsample.sam
#samtools view -@ 16 -b 10X_downsample.sam > 10X_downsample.bam &&rm filtered_SAM_body 10X_downsample.sam SAM_header

#java -jar $QoRT_root/QoRTs.jar QC --singleEnded --numThreads 12 10X_downsample.bam /share/home/hanxiaoping/tools/STAR_Reference_Human/Homo_sapiens.GRCh38.fa.gtf output_10X

###Smartseq3
#samtools view -@ 16 -H $Smartseq3_root/starAligned.out.bam > SAM_header
#samtools view $Smartseq3_root/starAligned.out.bam | shuf -n 9000000 > filtered_SAM_body
#cat SAM_header filtered_SAM_body > Smartseq3_downsample.sam
#samtools view -@ 16 -b Smartseq3_downsample.sam > Smartseq3_downsample.bam &&rm filtered_SAM_body Smartseq3_downsample.sam SAM_header

#java -jar $QoRT_root/QoRTs.jar QC --singleEnded --numThreads 12 Smartseq3_downsample.bam /share/home/hanxiaoping/tools/STAR_Reference_Human/Homo_sapiens.GRCh38.fa.gtf output_Smartseq3


###split-seq
samtools view -@ 16 -H $split_root/star_gene_exon_tagged.bam > SAM_header
samtools view $split_root/star_gene_exon_tagged.bam | shuf -n 9000000 > filtered_SAM_body
cat SAM_header filtered_SAM_body > split_300cell_downsample.sam
samtools view -@ 16 -b split_300cell_downsample.sam > split_300cell_downsample.bam &&rm filtered_SAM_body split_300cell_downsample.sam SAM_header

java -jar $QoRT_root/QoRTs.jar QC --singleEnded --numThreads 20 split_300cell_downsample.bam /share/home/hanxiaoping/tools/STAR_Reference_Human/Homo_sapiens.GRCh38.fa.gtf output_split_300cell

###split-seq2
samtools view -@ 16 -H $split_root2/star_gene_exon_tagged.bam > SAM_header
samtools view $split_root2/star_gene_exon_tagged.bam | shuf -n 9000000 > filtered_SAM_body
cat SAM_header filtered_SAM_body > split_3000cell_downsample.sam
samtools view -@ 16 -b split_3000cell_downsample.sam > split_3000cell_downsample.bam &&rm filtered_SAM_body split_3000cell_downsample.sam SAM_header

java -jar $QoRT_root/QoRTs.jar QC --singleEnded --numThreads 20 split_3000cell_downsample.bam /share/home/hanxiaoping/tools/STAR_Reference_Human/Homo_sapiens.GRCh38.fa.gtf output_split_3000cell

###DroNc-seq
samtools view -@ 16 -H $DroNc_root/star_gene_exon_tagged.bam > SAM_header
samtools view $DroNc_root/star_gene_exon_tagged.bam | shuf -n 9000000 > filtered_SAM_body
cat SAM_header filtered_SAM_body > DroNc_downsample.sam
samtools view -@ 16 -b DroNc_downsample.sam > DroNc_downsample.bam &&rm filtered_SAM_body DroNc_downsample.sam SAM_header

java -jar $QoRT_root/QoRTs.jar QC --singleEnded --numThreads 20 DroNc_downsample.bam /share/home/hanxiaoping/tools/STAR_Reference_Mouse/Mus_musculus.GRCm38.88.gtf output_DroNc

###sci-seq
samtools view -@ 16 -H $sci_root/star_gene_exon_tagged.bam > SAM_header
samtools view $sci_root/star_gene_exon_tagged.bam | shuf -n 9000000 > filtered_SAM_body
cat SAM_header filtered_SAM_body > sci_downsample.sam
samtools view -@ 16 -b sci_downsample.sam > sci_downsample.bam &&rm filtered_SAM_body sci_downsample.sam SAM_header

java -jar $QoRT_root/QoRTs.jar QC --singleEnded --numThreads 20 sci_downsample.bam /share/home/hanxiaoping/tools/STAR_Reference_Human/Homo_sapiens.GRCh38.fa.gtf output_sci