#!/bin/bash
#$ -S /bin/bash
#$ -o /dev/null
#$ -e /dev/null
#$ -cwd
#$ -V
#$ -l tmem=6G,h_vmem=6G
#$ -l h_rt=24:0:0
#$ -t 1-25
set -u
set -x
scriptname=annotate
mkdir -p ${scriptname}.qsub.out ${scriptname}.qsub.err
exec >${scriptname}.qsub.out/${scriptname}_${SGE_TASK_ID}_${JOB_ID}.out 2>${scriptname}.qsub.err/${scriptname}_${SGE_TASK_ID}_${JOB_ID}.err
args=( header 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y )
f=${args[$SGE_TASK_ID]}
chrCode=$f


function VEP() {
    DIR=~/bin/DNASeq_pipeline/annotation
    reference=1kg
    zcat chr${chrCode}.vcf.gz | python ${DIR}/multiallele_to_single_gvcf.py --GQ 10 --DP 5 > chr${chrCode}-single.vcf
    bash $DIR/run_VEP.sh --vcfin chr${chrCode}-single.vcf --chr $chrCode --reference $reference --vcfout VEP_${chrCode}.vcfout --coding_only yes --custom UCLEX,AJcontrols,AJcases,BroadAJcontrols,ImmunoBase
}


function annotate() {
    DIR=~/bin/DNASeq_pipeline/annotation
    python $DIR/processVEP.py --custom-allele-freq UCLEX AJcontrols AJcases BroadAJcontrols --custom-annotation ImmunoBase_CRO ImmunoBase_IBD ImmunoBase_UC --file VEP_${chrCode}.vcfout 
    #combine annotation, genotype and depth into a single VEP file.
    Rscript $DIR/combine.R --chr $chrCode | Rscript $DIR/annotate/annotate_expression.R | Rscript $DIR/annotate/annotate_pedigree_af.R| Rscript scripts/Group_Counts.R > VEP_${chrCode}.csv
}

function filter() {
    DIR=~/bin/DNASeq_pipeline/annotation/filters/
    cat VEP_$chrCode.csv | Rscript $DIR/af-filter.R --ajcontrols.thresh 0.05 --uclex.thresh 0.05 --exac.thresh 0.025 > af-filtered-VEP_$chrCode.csv
    wc -l  af-filtered-VEP_$chrCode.csv
    cat af-filtered-VEP_$chrCode.csv |  Rscript $DIR/csq-filter.R > af-csq-filtered-VEP_$chrCode.csv
    wc -l af-csq-filtered-VEP_$chrCode.csv
    cat af-csq-filtered-VEP_$chrCode.csv | Rscript $DIR/GO-filter.R > af-csq-go-filtered-VEP_$chrCode.csv 
    wc -l af-csq-go-filtered-VEP_$chrCode.csv 
}

function tonyfilter() {
    DIR=~/bin/DNASeq_pipeline/annotation/filters/
    cat VEP_$chrCode.csv | Rscript $DIR/af-filter.R --ajcontrols.thresh 0.05 --uclex.thresh 0.05 --exac.thresh 0.05 > af-filtered-VEP_$chrCode.csv
    wc -l  af-filtered-VEP_$chrCode.csv
    cat af-filtered-VEP_$chrCode.csv |  Rscript $DIR/csq-filter.R --cadd.thresh 10 > af-csq-filtered-VEP_$chrCode.csv
    wc -l af-csq-filtered-VEP_$chrCode.csv
    cat af-csq-filtered-VEP_$chrCode.csv | Rscript $DIR/expression-filter.R > af-csq-expression-filtered-VEP_$chrCode.csv 
    wc -l af-csq-expression-filtered-VEP_$chrCode.csv 
}


function stringent_filter() {
    DIR=~/bin/DNASeq_pipeline/annotation/filters/
    cat VEP_$chrCode.csv | Rscript $DIR/af-filter.R --ajcontrols.thresh 0.05 --uclex.thresh 0.05 --exac.thresh 0.025 > af-filtered-VEP_$chrCode.csv
    wc -l  af-filtered-VEP_$chrCode.csv
    cat af-filtered-VEP_$chrCode.csv |  Rscript $DIR/csq-filter.R > af-csq-filtered-VEP_$chrCode.csv
    wc -l af-csq-filtered-VEP_$chrCode.csv
    cat af-csq-filtered-VEP_$chrCode.csv | Rscript $DIR/GO-filter.R > af-csq-go-filtered-VEP_$chrCode.csv 
    wc -l af-csq-go-filtered-VEP_$chrCode.csv 
}



#VEP
#annotate
filter
#tonyfilter


