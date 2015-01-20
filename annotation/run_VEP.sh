#!/bin/bash

####INPUT FILE
#vcfin=/home/zchads1/cluster/UCL-exomes_v2/variants/Levine_${chr}_single.vcf
####SET CHROMOSOME
#ch1r=1
#out dir

set -u
set -e
set -x

# prints to stderr in red
function error() { >&2 echo -e "\033[31m$*\033[0m"; }
function stop() { error "$*"; exit 1; }
try() { "$@" || stop "cannot $*"; }



function usage() {
    echo --vcfin
    echo --chr
    echo --vcfout
    echo --reference
    exit 1
}


##
until [ -z "$1" ]
do
    # use a case statement to test vars. we always test $1 and shift at the end of the for block.
    case $1 in
    --vcfin)
        shift
        vcfin=$1;;
    --chr)
        shift
        chr=$1;;
    --vcfout)
        shift
        vcfout=$1;;
    --reference)
        shift
        reference=$1;;
    -* )
        echo "Unrecognized option: $1"
        usage
        exit 1;;
    esac
    shift
    if [ "$#" = "0" ]; then break; fi
done 


####CONFIGURE SOFTWARE SHORTCUTS AND PATHS
vep=/cluster/project8/vyp/AdamLevine/software/ensembl/src/ensembl-tools/scripts/variant_effect_predictor/variant_effect_predictor.pl
dir_cache=/cluster/project8/vyp/AdamLevine/software/ensembl/cache/
perl5142=/share/apps/perl-5.14.2/bin/perl
PERL5LIB=${PERL5LIB}:/cluster/project8/vyp/AdamLevine/software/ensembl/src/bioperl-1.6.1
PERL5LIB=${PERL5LIB}:/cluster/project8/vyp/AdamLevine/software/ensembl/src/ensembl/modules
PERL5LIB=${PERL5LIB}:/cluster/project8/vyp/AdamLevine/software/ensembl/src/ensembl-compara/modules
PERL5LIB=${PERL5LIB}:/cluster/project8/vyp/AdamLevine/software/ensembl/src/ensembl-variation/modules
PERL5LIB=${PERL5LIB}:/cluster/project8/vyp/AdamLevine/software/ensembl/src/ensembl-funcgen/modules
PERL5LIB=${PERL5LIB}:/cluster/project8/vyp/AdamLevine/software/ensembl/Plugins
export PERL5LIB
export PATH=$PATH:/cluster/project8/vyp/vincent/Software/tabix-0.2.5/
condel_config=/cluster/project8/vyp/AdamLevine/software/ensembl/Plugins/config/Condel/config

#fasta=/scratch2/vyp-scratch2/reference_datasets/human_reference_sequence/human_g1k_v37.fasta
#fasta=/scratch2/vyp-scratch2/reference_datasets/human_reference_sequence/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
#file:///scratch2/vyp-scratch2/reference_datasets/human_reference_sequence/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
cat $vcfin | grep '^##reference=' | cut -f2 -d'='

if [[ "$reference" == "hg38_noAlt" ]]
then
    #fasta=/scratch2/vyp-scratch2/reference_datasets/human_reference_sequence/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
    fasta=/cluster/project8/vyp/pontikos/data/reference_genomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
    chrPrefix='chr'
    assembly=GRCh38
    port=
    custom_annotation=
elif [[ "$reference" == "1kg" ]]
then
    #fasta=/scratch2/vyp-scratch2/reference_datasets/human_reference_sequence/human_g1k_v37.fasta
    fasta=/cluster/project8/vyp/pontikos/data/reference_genomes/human_g1k_v37.fasta
    assembly=GRCh37
    chrPrefix=''
    port='--port 3337'
    ####ANNOTATIONS for GRCh37
    ####Directory
    annotations_dir=/cluster/project8/vyp/AdamLevine/annotations
    ####CADD http://cadd.gs.washington.edu/home
    #This needs to be updated with the latest scores [ACTION!]
    #They also now provide a script which is worth exploring
    cadd=/scratch2/vyp-scratch2/AdamLevine/CADD/cadd.vcf.gz
    ####ADD ExAC [ACTION!]
    ####ESP frequency annotations
    EA=${annotations_dir}/esp/results/EA_AF_${chr}.vcf.gz
    AA=${annotations_dir}/esp/results/AA_AF_${chr}.vcf.gz
    ####UCLex frequencies and need to be updated [ACTION!]
    VP=/cluster/project8/vyp/AdamLevine/UCL-exomes_v2/VP_cohort/UCLfreq_${chr}.vcf.gz 
    ####1KG
    #Do not need all of these different annotations [ACTION!]
    OneKG=/cluster/project8/vyp/AdamLevine/exome/annotations/OneKG/OneKG_AF_${chr}.vcf.gz
    OneKGceugbr=/cluster/project8/vyp/AdamLevine/exome/annotations/OneKG_CEUGBR/OneKG-CEUGBR_AF_${chr}.vcf.gz
    okg=${annotations_dir}/1kg/results/1KG_AF_${chr}.vcf.gz
    okgEUR=${annotations_dir}/1kg/results/1KG_EUR-AF_${chr}.vcf.gz
    okgAFR=${annotations_dir}/1kg/results/1KG_AFR-AF_${chr}.vcf.gz
    okgAMR=${annotations_dir}/1kg/results/1KG_AMR-AF_${chr}.vcf.gz
    okgASN=${annotations_dir}/1kg/results/1KG_ASN-AF_${chr}.vcf.gz
    #these options are only available for this build at the moment
    custom_annotation="--custom $cadd,CADD,vcf,exact --custom $EA,EAf,vcf,exact --custom $AA,AAf,vcf,exact --custom $VP,VPf,vcf,exact --custom $OneKG,ONEKGf,vcf,exact --custom $OneKGceugbr,CEUBGRf,vcf,exact --custom $okg,okg,vcf,exact --custom $okgEUR,okgEUR,vcf,exact --custom $okgAMR,okgAMR,vcf,exact --custom $okgAFR,okgAFR,vcf,exact --custom $okgASN,okgASN,vcf,exact"
elif [[ "$reference" == "hg19" ]]
then
    fasta=/scratch2/vyp-scratch2/reference_datasets/human_reference_sequence/hg19_UCSC.fa
    chrPrefix='chr'
    port=
    custom_annotation=
else
    stop Unsupported reference $reference
fi

maf="--maf_esp --gmaf --maf_1kg"
#fields="--fields SYMBOL,CLIN_SIG,AA_MAF,EA_MAF,SIFT,PolyPhen,CAROL,Condel"
fields=''

$perl5142 $vep $port --ASSEMBLY $assembly --fasta $fasta --cache --dir_cache $dir_cache --input_file $vcfin --format vcf --sift b --polyphen b --symbol --coding_only --canonical --check_existing --check_alleles --plugin Carol --stats_text --no_progress --output_file $vcfout --plugin Condel,${condel_config},b --force_overwrite --pick  --fork 2 $maf $fields $custom_annotation

##--ASSEMBLY GRCh37  \
#--ASSEMBLY GRCh38  \
#--port 3337 \
#--cache \
#--dir_cache $dir_cache \
#--input_file $vcfin \
#--format vcf \
#--sift b \
#--polyphen b \
#--symbol \
#--coding_only \
#--canonical \
#--check_existing \
#--check_alleles \
#--plugin Carol \
##--plugin Condel,${condel_config},b \
##--custom $cadd,CADD,vcf,exact \
##--custom $EA,EAf,vcf,exact \
##--custom $AA,AAf,vcf,exact \
##--custom $VP,VPf,vcf,exact \
##--custom $OneKG,ONEKGf,vcf,exact \
##--custom $OneKGceugbr,CEUBGRf,vcf,exact \
##--custom $okg,okg,vcf,exact \
##--custom $okgEUR,okgEUR,vcf,exact \
##--custom $okgAMR,okgAMR,vcf,exact \
##--custom $okgAFR,okgAFR,vcf,exact \
##--custom $okgASN,okgASN,vcf,exact \
#--stats_text \
#--vcf  \
#--no_progress \
#--output_file $vcfout\


