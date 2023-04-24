#BSUB -q normal
#BSUB -J cpdge
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -R "span[ptile=6]"
#BSUB -n 18

cd /share/home/hanxiaoping/RAWDATA/OB/E150006510/process/
str1="_dge.summary.txt"

str4="col"
mkdir summary
for i in {13..14} {32..35} 38 43 44 {61..63} {69..71}; do
    cd $i;
    cp _dge.summary.txt ../summary/$str4$i$str1

    cd ../
done













