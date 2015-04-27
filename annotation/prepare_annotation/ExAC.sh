#!/bin/bash
#$ -S /bin/bash
#$ -o /dev/null
#$ -e /dev/null
#$ -cwd
#$ -V
#$ -R y
#$ -pe smp 1
#$ -l scr=1G
#$ -l tmem=2G,h_vmem=1G
#$ -l h_rt=1:0:0
#$ -t 1-25
#$ -tc 25
set -u
set -x
let "SGE_TASK_ID=$SGE_TASK_ID-1"
scriptname=`basename $0`
mkdir -p ${scriptname}.qsub.out ${scriptname}.qsub.err
exec >${scriptname}.qsub.out/${scriptname}_${SGE_TASK_ID}_${JOB_ID}.out 2>${scriptname}.qsub.err/${scriptname}_${SGE_TASK_ID}_${JOB_ID}.err
args=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y )
f=${args[$SGE_TASK_ID]}
#vcf files must exist
#cat header.txt $f > chr$f.vcf
#DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
DIR=~/bin/pipelines/annotation/
for pop in AFR AMR Adj EAS FIN NFE OTH SAS
do
    vcftools --vcf chr$f.vcf --get-INFO AC_${pop} --get-INFO AN_${pop} --out chr${f}_${pop}
    python $DIR/INFO2VCF.py < chr${f}_${pop}.INFO | bgzip > chr${f}_${pop}.vcf.gz
    tabix -f -p vcf chr${f}_${pop}.vcf.gz
done

