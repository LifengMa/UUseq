#BSUB -q normal
#BSUB -J creatdir
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -R "span[ptile=16]"
#BSUB -n 16
cd /share/home/hanxiaoping/RAWDATA/OB/E150002201/L01/process/
str1="_1.fq.gz"
str2="_2.fq.gz"
str3="E150002201_L01_"
for i in {45..48} {57..64} {73..80} {85..88} {121..124}  ; do
    cd /share/home/hanxiaoping/RAWDATA/OB/E150002201/L01/process/
    mkdir $i
    mv /share/home/hanxiaoping/RAWDATA/OB/E150002201/L01/OUTPUT/$str3$i$str1 $i
    mv /share/home/hanxiaoping/RAWDATA/OB/E150002201/L01/OUTPUT/$str3$i$str2 $i
done












