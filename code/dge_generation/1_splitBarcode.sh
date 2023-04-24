#BSUB -q normal
#BSUB -J splitbarcode
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -R "span[ptile=18]"
#BSUB -n 72

/share/home/hanxiaoping/tools/splitBarcode-master/V0.1.6_release/linux/splitBarcode \
/share/home/hanxiaoping/tools/splitBarcode-master/index_HUADA.txt H_R1.fq.gz -2 H_R2.fq.gz \
-o OUTPUT -b 250 10 2 -r

