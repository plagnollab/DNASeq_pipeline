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
#scriptname=`basename $0`
scriptname=liftOver
mkdir -p ${scriptname}.qsub.out ${scriptname}.qsub.err
exec >${scriptname}.qsub.out/${scriptname}_${SGE_TASK_ID}_${JOB_ID}.out 2>${scriptname}.qsub.err/${scriptname}_${SGE_TASK_ID}_${JOB_ID}.err
args=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y )
ch=${args[$SGE_TASK_ID]}
liftOver=/cluster/project8/vyp/AdamLevine/software/LiftOver/liftOver
map=/cluster/project8/vyp/pontikos/hg19ToHg38.over.chain.gz
for pop in AFR AMR Adj EAS FIN NFE OTH SAS
do
    vcfgz=chr${ch}_${pop}.vcf.gz
    out=${vcfgz%.vcf.gz}-b38.vcf.gz
    bed=${vcfgz%.gz}.bed
    mbed=${vcfgz%.gz}-b38.bed
    zcat $vcfgz | awk '{OFS="\t"; if (!/^#/){print $1,$2-1,$2,$3,$4,$5,$6}}' | sed 's/^/chr/' > $bed
    $liftOver $bed $map  $mbed ${bed%.bed}-b38-unlifted.bed
    #echo output `wc -l $input-b38.bed`
    #echo unlifted `wc -l $input-b38-unlifted.bed`
    cat <(zgrep '^#' $vcfgz) <(cat $mbed | sed 's/^chr//' | awk '{OFS="\t"; print $1,$3,$4,$5,$6,".",".","."}' | sort -k2 -n | grep ^$ch) | bgzip  > $out
    tabix -f -p vcf $out
done


