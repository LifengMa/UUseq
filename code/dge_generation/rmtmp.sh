#BSUB -q normal
#BSUB -J rm
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -R "span[ptile=6]"
#BSUB -n 18

cd /share/home/hanxiaoping/RAWDATA/OB/E150006510/process/
for i in {13..14} {32..35} 38 43 44 {61..63} {69..71}; do
    cd $i;
    rm starAligned.out.bam H_R1.fq.gz H_R2.fq.gz tmp/R2_polyA_trim_reverse.fastq

    cd ../
done













