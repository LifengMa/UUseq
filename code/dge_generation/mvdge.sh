#BSUB -q normal
#BSUB -J cpdge
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -R "span[ptile=6]"
#BSUB -n 18

cd /share/home/hanxiaoping/RAWDATA/OB/E150006510/process/
str1="_dge.txt.gz"
str2="_type.txt"
str3="_out_cell_readcounts.txt"
str4="col"
mkdir dge
for i in {13..14} {32..35} 38 43 44 {61..63} {69..71}; do
    cd $i;
    cp _dge.txt.gz ../dge/$str4$i$str1
    cp type.txt ../dge/$str4$i$str2
    cp out_cell_readcounts.txt ../dge/$str4$i$str3
    cd ../
done













