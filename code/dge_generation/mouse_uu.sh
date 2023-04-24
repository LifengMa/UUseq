#BSUB -q normal
#BSUB -J UU
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -R "span[ptile=10]"
#BSUB -n 10

cd `pwd`
mkdir tmp
tmpdir=tmp
sample_name=$(basename `pwd`)
dropseq_root=/share/home/hanxiaoping/tools/Drop-seq_tools-2.5.1

# fastq --> bam
java -jar ${dropseq_root}/3rdParty/picard/picard.jar FastqToSam F1=H_R1.fq.gz F2=H_R2.fq.gz \
  O=H.bam QUALITY_FORMAT=Standard SAMPLE_NAME=sample_name TMP_DIR=$tmpdir

#########################Paired################################

##-------------  Cell Barcode -------------
#add RT barcode
#  #tag barcodes1
${dropseq_root}/TagBamWithReadSequenceExtended \
BASE_RANGE=1-20 BASE_QUALITY=10 BARCODED_READ=1 TAG_BARCODED_READ=false DISCARD_READ=true TAG_NAME=PC NUM_BASES_BELOW_QUALITY=4 \
INPUT=H.bam OUTPUT=$tmpdir/H1.bam COMPRESSION_LEVEL=5
${dropseq_root}/TagBamWithReadSequenceExtended \
BASE_RANGE=1-17 BASE_QUALITY=10 BARCODED_READ=1 TAG_BARCODED_READ=true DISCARD_READ=false TAG_NAME=RT NUM_BASES_BELOW_QUALITY=4 \
INPUT=$tmpdir/H1.bam OUTPUT=$tmpdir/H2.bam COMPRESSION_LEVEL=5 && rm $tmpdir/H1.bam 




# #FilterBAM
${dropseq_root}/FilterBam TAG_REJECT=XQ INPUT=$tmpdir/H2.bam OUTPUT=$tmpdir/H3.bam && rm $tmpdir/H2.bam

# # corrected bam for one mismatch
barcodepath=/share/home/hanxiaoping/tools/UU_barcode/

python ${barcodepath}/UU_correct_1124.py $barcodepath $tmpdir/H3.bam $tmpdir/filtered.bam && rm $tmpdir/H3.bam
# echo "correct sam files done"
# ########################

java -Xmx100g -jar /share/home/hanxiaoping/tools/Drop-seq_tools-2.5.1/3rdParty/picard/picard.jar\
 SamToFastq INPUT=$tmpdir/filtered.bam FASTQ=$tmpdir/R2.fastq  READ1_TRIM=17

# # ## PolyATrimmer

cutadapt -a A{10} -j 0 -O 10 --minimum-length=20 -o $tmpdir/R2_trim.fastq \
$tmpdir/R2.fastq && rm $tmpdir/R2.fastq

/share/home/hanxiaoping/tools/seqtk/seqtk seq -Ar tmp/R2_trim.fastq > tmp/R2_polyA_trim_reverse.fastq && rm tmp/R2_trim.fastq
## Alignment STAR
/share/home/hanxiaoping/tools/STAR-2.5.2a/source/STAR \
--genomeDir /share/home/hanxiaoping/tools/STAR_Reference_Mouse/genomeDir/ \
--readFilesIn $tmpdir/R2_polyA_trim_reverse.fastq \
--outFileNamePrefix star \
--limitOutSJcollapsed 5000000 \
--outSAMtype BAM Unsorted 

## MergeBamAlignment
java -Xmx100g -jar ${dropseq_root}/3rdParty/picard/picard.jar MergeBamAlignment \
REFERENCE_SEQUENCE=/share/home/hanxiaoping/tools/STAR_Reference_Mouse/genomeDir/Mus_musculus.GRCm38.88.fasta \
UNMAPPED_BAM=$tmpdir/filtered.bam \
ALIGNED_BAM=starAligned.out.bam \
OUTPUT=merged.bam \
INCLUDE_SECONDARY_ALIGNMENTS=false \
PAIRED_RUN=false \
&& rm $tmpdir/filtered.bam

 ${dropseq_root}/TagReadWithGeneFunction I=merged.bam \
 O=star_gene_exon_tagged.bam \
 ANNOTATIONS_FILE=/share/home/hanxiaoping/tools/STAR_Reference_Mouse/genomeDir/Mus_musculus.GRCm38.88.gtf \
&& rm merged.bam

${dropseq_root}/DigitalExpression -m 8g \
I=star_gene_exon_tagged.bam \
CELL_BARCODE_TAG=XC \
MOLECULAR_BARCODE_TAG=XM \
O=_dge.txt.gz \
SUMMARY=_dge.summary.txt \
NUM_CORE_BARCODES=20000 \
LOCUS_FUNCTION_LIST=INTRONIC \
TMP_DIR=.

##  Read Summary
${dropseq_root}/BamTagHistogram I=star_gene_exon_tagged.bam O=out_cell_readcounts.txt TAG=XC 
samtools view star_gene_exon_tagged.bam |awk '{for (i=1; i<=NF; ++i) {if($i ~ "^XF:Z:"){print $i}}}' >type.txt