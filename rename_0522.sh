#BSUB -q normal
#BSUB -J rename
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -R "span[ptile=18]"
#BSUB -n 18

cd /share/home/hanxiaoping/RAWDATA/OB/E150002201/L01/process/

str1="H_R1.fq.gz"
str2="H_R2.fq.gz"
str3="_1.fq.gz"
str4="_2.fq.gz"
str5="E150002201_L01_"

for i in {45..48} {57..64} {73..80} {85..88} {121..124}; do
    cd $i
    mv $str5$i$str3 $str1
    mv $str5$i$str4 $str2
 cd ../
done












