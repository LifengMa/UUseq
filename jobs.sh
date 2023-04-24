#BSUB -q normal
#BSUB -J Jobs
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -R "span[ptile=1]"
#BSUB -n 1

cd 'pwd'

for i in {13..14} {32..35} 38 43 44 {61..63} {69..71}; do
    cd $i;
    cp ../mouse_uu.sh .;
    bsub < ./mouse_uu.sh;
    cd ../;
done
