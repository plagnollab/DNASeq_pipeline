# *************************************************************
#
# $Source: $
# $Revision: $                                                                 
# $State: $                                                                     
# $Date: $                                                      
# $Author: $  
#
# $Log: $
#
#
# *************************************************************

#!/bin/bash

#$ -S /bin/bash
#$ -o /home/zchads1/cluster/UCL-exomes_v2/variants/log
#$ -e /home/zchads1/cluster/UCL-exomes_v2/variants/log
#$ -l h_rt=12:0:0
#$ -l tmem=4G,h_vmem=4G
#$ -t 1-24

i=$(($SGE_TASK_ID - 1))
chromosomes=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y)
chr=${chromosomes[$i]}

gvcf=/home/zchads1/cluster/UCL-exomes_v2/by_chr/chr_${chr}.vcf.gz

vcftools=~/cluster/software/vcftools_0.1.12b/bin/vcftools

cd /home/zchads1/cluster/UCL-exomes_v2/variants/

$vcftools --gzvcf $gvcf \
--keep ../Levine_samples.txt \
--chr ${chr} \
--mac 1 \
--minGQ 20 \
--minDP 10 \
--max-missing 0.7 \
--recode \
--out Levine_${chr} > log/Levine_${chr}.log

gzip Levine_${chr}.recode.vcf
