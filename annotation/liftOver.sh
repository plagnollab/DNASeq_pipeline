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
args=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y )
ch=${args[$SGE_TASK_ID]}
exec >${scriptname}.qsub.out/${scriptname}_chr${ch}_${JOB_ID}.out 2>${scriptname}.qsub.err/${scriptname}_chr${ch}_${JOB_ID}.err

# ann is one of ExAC, esp etc
ann=$1
ann_dir=/cluster/project8/IBDAJE/VEP_custom_annotations/
liftOver=/cluster/project8/vyp/AdamLevine/software/LiftOver/liftOver
#map1=$ann_dir/chain_files/GRCh37ToHg19.over.chain.gz
#map2=$ann_dir/chain_files/hg19ToHg38.over.chain.gz 
map=$ann_dir/chain_files/GRCh37_to_GRCh38.chain.gz
grc37=$ann_dir/GRCh37/$ann
grc38=$ann_dir/GRCh38/$ann

#For UCLex and CADD we need this instead as we don't have any pop specific AF:
if [[ "$ann" == "UCLex" ]] || [[ "$ann" == "CADD" ]] 
then
    files=$grc37/chr${ch}.vcf.gz
else
    files=$grc37/chr${ch}_*.vcf.gz
fi
for vcfgz in $files
do
    [[ ! -s $vcfgz ]] && continue
    f=`basename ${vcfgz%.vcf.gz}`
    inbed=$grc37/${f}.bed
    outbed=$grc38/${f}.bed
    unliftedbed=$grc38/${f}_unlifted.bed
    outvcf=$grc38/${f}.vcf
    header=$grc38/header.txt
    outvcfgz=${outvcf}.gz
    if [[ ! -s $inbed ]]
    then
        zcat $vcfgz | awk '{OFS="\t"; if (!/^#/){print $1,$2-1,$2,$3,$4,$5,$6}}'  > $inbed
    fi
    # GRCh37 -> GRCh38
    if [[ ! -s $outbed ]]
    then
        $liftOver $inbed $map  $outbed $unliftedbed
    fi
    echo output `wc -l $outbed`
    echo unlifted `wc -l $unliftedbed`
    if [[ ! -s $header ]]
    then
        zgrep '^#' $vcfgz > $header
    fi
    sort -T /scratch0/ -k2 -n $outbed | awk -v chr=$ch '{OFS="\t";if ($1==chr){print $1,$3,$4,$5,$6,".",".","."}}' > $outbed.filter
    #filter on chromosome name to remove positions mapped to different chrom or alternative sequences
    cat $header $outbed.filter > $outvcf
    bgzip -f -c $outvcf > $outvcfgz
    tabix -f -p vcf $outvcfgz
done

