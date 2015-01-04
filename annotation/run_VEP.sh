#!/bin/bash

####SET CHROMOSOME
chr=1

####INPUT FILE
vcfin=/home/zchads1/cluster/UCL-exomes_v2/variants/Levine_${chr}_single.vcf

####CONFIGURE SOFTWARE SHORTCUTS AND PATHS
vep=/cluster/project8/vyp/AdamLevine/software/ensembl/src/ensembl-tools/scripts/variant_effect_predictor/variant_effect_predictor.pl
dir_cache=/cluster/project8/vyp/AdamLevine/software/ensembl/cache/
perl5142=/share/apps/perl-5.14.2/bin/perl
PERL5LIB=${PERL5LIB}:/cluster/project8/vyp/AdamLevine/software/ensembl/src/bioperl-1.6.1
PERL5LIB=${PERL5LIB}:/cluster/project8/vyp/AdamLevine/software/ensembl/src/ensembl/modules
PERL5LIB=${PERL5LIB}:/cluster/project8/vyp/AdamLevine/software/ensembl/src/ensembl-compara/modules
PERL5LIB=${PERL5LIB}:/cluster/project8/vyp/AdamLevine/software/ensembl/src/ensembl-variation/modules
PERL5LIB=${PERL5LIB}:/cluster/project8/vyp/AdamLevine/software/ensembl/src/ensembl-funcgen/modules
export PERL5LIB
export PATH=$PATH:/cluster/project8/vyp/vincent/Software/tabix-0.2.5/
condel_config=/cluster/project8/vyp/AdamLevine/software/ensembl/Plugins/config/Condel/config

####ANNOTATIONS
	####Directory
annotations_dir=/cluster/project8/vyp/AdamLevine/annotations
	####CADD http://cadd.gs.washington.edu/home
	#This needs to be updated with the latest scores [ACTION!]
	#They also now provide a script which is worth exploring
cadd=/scratch2/vyp-scratch2/AdamLevine/CADD/cadd.vcf.gz
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

####ADD ExAC [ACTION!]
##Prepare Build 38 allele frequency annotations [ACTION!]

####RUN VEP
$perl5142 $vep \
--ASSEMBLY GRCh37 \ ###For Build 37, use GRCh38 for Build 38
--port 3337 \ ###For Build 37, remove for Build 38
--cache \
--dir_cache $dir_cache \
--input_file $vcfin \
--format vcf \
--sift b \
--polyphen b \
--symbol \
--coding_only \
--canonical \
--check_existing \
--check_alleles \
--plugin Carol \
--plugin Condel,${condel_config},b \
--custom $cadd,CADD,vcf,exact \
--custom $EA,EAf,vcf,exact \
--custom $AA,AAf,vcf,exact \
--custom $VP,VPf,vcf,exact \
--custom $OneKG,ONEKGf,vcf,exact \
--custom $OneKGceugbr,CEUBGRf,vcf,exact \
--custom $okg,okg,vcf,exact \
--custom $okgEUR,okgEUR,vcf,exact \
--custom $okgAMR,okgAMR,vcf,exact \
--custom $okgAFR,okgAFR,vcf,exact \
--custom $okgASN,okgASN,vcf,exact \
--stats_text \
--vcf \
--output_file VEP_${chr}.vcf \
--no_progress 
