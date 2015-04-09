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
args=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y )
ch=${args[$SGE_TASK_ID]}
#scriptname=`basename $0`
scriptname=liftOver
mkdir -p ${scriptname}.qsub.out ${scriptname}.qsub.err
exec >${scriptname}.qsub.out/${scriptname}_ch${ch}_${JOB_ID}.out 2>${scriptname}.qsub.err/${scriptname}_ch${ch}_${JOB_ID}.err
liftOver=/cluster/project8/vyp/AdamLevine/software/LiftOver/liftOver
map=/cluster/project8/vyp/pontikos/hg19ToHg38.over.chain.gz 
#
grc37=/goon2/scratch2/vyp-scratch2/annotations/GRCh37/CADD
grc38=/goon2/scratch2/vyp-scratch2/annotations/GRCh38/CADD
invcfgz=$grc37/cadd.vcf.gz
inbed=$grc37/chr${ch}.bed
outbed=$grc38/chr${ch}.bed
unliftedbed=$grc38/chr${ch}.bed
outvcfgz=$grc38/chr${ch}.vcf.gz
#zcat $vcfgz | awk '{OFS="\t"; if (!/^#/){print $1,$2-1,$2,$3,$4,$5,$6}}' | sed 's/^/chr/' > $bed
$liftOver $inbed $map  $outbed $unliftedbed
#echo output `wc -l $input-b38.bed`
#echo unlifted `wc -l $input-b38-unlifted.bed`
# The standard header format is:
##fileformat=VCFv4.0
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
cat <(zgrep '^#' $invcfgz) <(sort -k 2 -n $outbed | sed 's/^chr//' | awk -v chr=$ch '{OFS="\t";if ($1==chr){print $1,$3,$4,$5,$6,".",".","."}}') | bgzip  > $outbed
tabix -f -p vcf $outbed


