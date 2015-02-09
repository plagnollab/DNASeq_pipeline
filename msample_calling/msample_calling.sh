#! /bin/bash

# prints to stderr in red
function error() { >&2 echo -e "\033[31m$*\033[0m"; }
function stop() { error "$*"; exit 1; }
try() { "$@" || stop "cannot $*"; }

fasta=/scratch2/vyp-scratch2/reference_datasets/human_reference_sequence/human_g1k_v37.fasta
bundle=/scratch2/vyp-scratch2/reference_datasets/GATK_bundle

Rscript=/cluster/project8/vyp/vincent/Software/R-3.1.2/bin/Rscript
Rbin=/cluster/project8/vyp/vincent/Software/R-3.1.2/bin/R

java=/share/apps/jdk1.7.0_45/bin/java
tmpDir=/scratch0/vyp
target=/cluster/project8/vyp/exome_sequencing_multisamples/target_region/data/merged_exome_target_cleaned.bed

GATK=/cluster/project8/vyp/vincent/Software/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar

submit=no

force=no
genotype=no
recal=no
gVCFlist=none

maxGaussians=6
#maxGaussiansIndels=5
numBad=1000
#numBadIndels=1000
GQ=20

#output=/scratch2/vyp-scratch2/vincent/GATK/cardioset_${currentUCLex}/cardioset_${currentUCLex}

until [ -z "$1" ]; do
    # use a case statement to test vars. we always test $1 and shift at the end of the for block.
    case $1 in
	--force )
	    shift
	    force=$1;;
	--target )
            shift
            target=$1;;
        --tmpDir )
	    shift
	    tmpDir=$1;;
	--genotype )
	    shift
	    genotype=$1;;
	--recal )
	    shift
	    recal=$1;;
	--annovar )
	    shift
	    annovar=$1;;
	--convertToR )
	    shift
	    convertToR=$1;;
	--gVCFlist )
	    shift
	    gVCFlist=$1;;
	--currentUCLex )
	    shift
	    currentUCLex=$1;;
        -* )
            echo "Unrecognized option: $1"
            exit 1;;
    esac
    shift
    if [ "$#" = "0" ]; then break; fi
done




############## Now options are all set
output=/scratch2/vyp-scratch2/vincent/GATK/mainset_${currentUCLex}/mainset_${currentUCLex}

### Check format of support file.
##should accept tab or space as delimiters
## but does read support tabs and delimeters?
mustBePath=`head -n1 $gVCFlist | cut -f1 -d' ' | cut -f1`
mustBeId=`head -n1 $gVCFlist | cut -f2 -d' ' | cut -f2`

if [[ "$mustBePath" != "path" ]]; then stop "The first column of the file $gVCFlist must have the name path $mustBePath"; fi
if [[ "$mustBeId" != "id" ]]; then stop "The second column of the file $gVCFlist must have the name id $mustBeId"; fi

memoSmall=5
memo=7.9


if [[ "$convertToR" == "yes" ]]; then memo=21.9; fi



mainScript=cluster/submission/calling.sh
## individual scripts of the form cluster/submission/subscript_chr${chr}.sh

echo "
#$ -o cluster/out
#$ -e cluster/error
#$ -S /bin/bash
#$ -l h_vmem=${memo}G,tmem=${memo}G
#$ -l h_rt=96:0:0
#$ -R y
#$ -pe smp 1
#$ -cwd 
#$ -t 1-23
#$ -tc 23

LISTCHROMS=(chr 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X )

CHR=\${LISTCHROMS[ \$SGE_TASK_ID ]}

sh cluster/submission/subscript_chr\${CHR}.sh

" > $mainScript


######## first clean up the individual files
for chr in `seq 1 22` X; do
    if [ -e cluster/submission/subscript_chr${chr}.sh ]; then rm cluster/submission/subscript_chr${chr}.sh; fi
done

##################################################
if [[ "$genotype" == "yes" ]]; then
    echo "Running the genotype module"
    
    for chr in `seq 1 22` X; do
	script=cluster/submission/subscript_chr${chr}.sh
	
	if [ ! -e ${output}_chr${chr}.vcf.gz.tbi ]; then 

	    echo "
$java -Xmx${memoSmall}g -jar $GATK \\
   -R $fasta \\
   -T GenotypeGVCFs \\
   -L $chr -L $target --interval_set_rule INTERSECTION --interval_padding 100  \\
   --annotation InbreedingCoeff --annotation QualByDepth --annotation HaplotypeScore --annotation MappingQualityRankSumTest --annotation ReadPosRankSumTest --annotation FisherStrand \\
   --dbsnp ${bundle}/dbsnp_137.b37.vcf \\" >> cluster/submission/subscript_chr${chr}.sh
	    
	    while read path id format; do
		if [[ "$format" == "v1" ]]; then gVCF=${path}/chr${chr}/${id}.gvcf.gz; fi
		if [[ "$format" == "v2" ]]; then gVCF=${path}/${id}-chr${chr}.gvcf.gz; fi

		echo "Including $gVCF"
		
		if [ ! -s $gVCF ]; then stop "Cannot find $gVCF"; fi
		if [ ! -s $gVCF.tbi ]; then stop "Cannot find $gVCF.tbi"; fi
		
		echo "   --variant $gVCF \\" >> $script
	    done < <(tail -n +2 $gVCFlist)
	    echo "   -o ${output}_chr${chr}.vcf.gz" >> $script
	fi
    done 
    
fi


##################################################
if [[ "$recal" == "yes" ]]; then
    echo "Running the recalibration module"
    

    for chr in `seq 1 22` X; do
	
	script=cluster/submission/subscript_chr${chr}.sh
	
	if [[ ! -s ${output}_chr${chr}_filtered.vcf || "$force" == "yes" ]]; then 
	    
	#### creates the tmpDir if needed
	    tmpDir=/scratch0/GATK_chr${chr}
	    

	    
	    echo "

if [ ! -e $tmpDir ]; then mkdir $tmpDir; fi

#### extract the indels
$java  -Djava.io.tmpdir=${tmpDir} -Xmx${memoSmall}g -jar ${GATK} \
     -T SelectVariants \
     -R $fasta \
     -V ${output}_chr${chr}.vcf.gz \
     -selectType INDEL \
     -selectType MIXED \
     -o ${output}_chr${chr}_indels.vcf.gz

#### apply the filters for the indels
$java -Djava.io.tmpdir=${tmpDir} -Xmx${memoSmall}g -jar ${GATK} \
    -T VariantFiltration \
    -R $fasta \
    -V ${output}_chr${chr}_indels.vcf.gz \
    --filterExpression \"QD < 2.0 || FS > 50.0 || ReadPosRankSum < -20.0\" \
    --filterName \"FAIL\" \
    -o ${output}_chr${chr}_indels_filtered.vcf.gz


#### extract the SNPs
$java  -Djava.io.tmpdir=${tmpDir} -Xmx${memoSmall}g -jar ${GATK} \
     -T SelectVariants \
     -R $fasta \
     -V ${output}_chr${chr}.vcf.gz \
     -selectType SNP \
     -o ${output}_chr${chr}_SNPs.vcf.gz



####### first SNPs
$java -Djava.io.tmpdir=${tmpDir} -Xmx${memoSmall}g -jar ${GATK} -T VariantRecalibrator -R $fasta --input ${output}_chr${chr}_SNPs.vcf.gz --maxGaussians ${maxGaussians} --mode SNP \
             -resource:hapmap,VCF,known=false,training=true,truth=true,prior=15.0 ${bundle}/hapmap_3.3.b37.vcf  \
             -resource:omni,VCF,known=false,training=true,truth=false,prior=12.0 ${bundle}/1000G_omni2.5.b37.vcf \
             -resource:dbsnp,VCF,known=true,training=false,truth=false,prior=8.0 ${bundle}/dbsnp_137.b37.vcf \
             -an QD -an FS -an ReadPosRankSum -an InbreedingCoeff \
             -tranche 100.0 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
             --minNumBadVariants ${numBad} \
             -recalFile ${output}_chr${chr}_SNPs_combrec.recal \
             -tranchesFile ${output}_chr${chr}_SNPs_combtranch \
             -rscriptFile  ${output}_chr${chr}_recal_plots_snps.R

${Rscript} ${output}_chr${chr}_recal_plots_snps.R


#apply_recal
$java -Xmx${memoSmall}g -jar ${GATK} -T ApplyRecalibration -R $fasta \
       -o ${output}_chr${chr}_SNPs_filtered.vcf.gz \
       --ts_filter_level 99.5 \
       --recal_file ${output}_chr${chr}_SNPs_combrec.recal --tranches_file ${output}_chr${chr}_SNPs_combtranch --mode SNP \
       --input ${output}_chr${chr}_SNPs.vcf.gz


#### Now we merge SNPs and indels
$java -Djava.io.tmpdir=${tmpDir} -Xmx${memoSmall}g -jar ${GATK} \
       -T CombineVariants --assumeIdenticalSamples \
       -R $fasta \
       --variant:SNPs ${output}_chr${chr}_SNPs_filtered.vcf.gz \
       --variant:indels ${output}_chr${chr}_indels_filtered.vcf.gz \
       -genotypeMergeOptions PRIORITIZE  \
       -priority SNPs,indels \
       -o ${output}_chr${chr}_filtered.vcf

rm -rf $tmpDir

rm ${output}_chr${chr}_indels.vcf.gz ${output}_chr${chr}_SNPs.vcf.gz ${output}_chr${chr}_SNPs_filtered.vcf.gz

" >> $script
	fi

    done
    
fi


if [[ "$annovar" == "yes" ]]; then


    for chr in `seq 1 22` X; do
	
	script=cluster/submission/subscript_chr${chr}.sh
	
	if [[ ! -e ${output}_snpStats/chr${chr}.done || "$force" == "yes" ]]; then  ## this is not quite right, needs fixing because it does not account for the last step
	    
	    echo "
if [ -e ${output}_${snpStats}/chr${chr}.done ]; then rm ${output}_${snpStats}/chr${chr}.done; fi  ## this is basically a log file, to make sure the job got finished

cut -f1-8 ${output}_chr${chr}_filtered.vcf > ${output}_chr${chr}_for_annovar.vcf

/cluster/project8/vyp/vincent/Software_heavy/annovar_Feb2013/convert2annovar.pl --allallele -format vcf4 --includeinfo ${output}_chr${chr}_for_annovar.vcf > ${output}_chr${chr}_db

/cluster/project8/vyp/vincent/Software_heavy/annovar_Feb2013/summarize_annovar_VP.pl -ver1000g 1000g2012apr -verdbsnp 137 -veresp 6500si -alltranscript -buildver hg19 --genetype ensgene --remove ${output}_chr${chr}_db /cluster/project8/vyp/vincent/Software_heavy/annovar_Feb2013/humandb_hg19/

perl ~/Software/pipeline/GATK_v2/custom_filtering.pl ${output}_chr${chr}_filtered.vcf ${output}_chr${chr}_recal_filtered2.vcf ${GQ}

python /cluster/project8/vyp/vincent/Software/pipeline/GATK_v2/annovar_vcf_combine_VP.py ${output}_chr${chr}_recal_filtered2.vcf ${output}_chr${chr}_db.exome_summary.csv ${output}_chr${chr}_exome_table.csv

perl /cluster/project8/vyp/vincent/Software/pipeline/msample_calling/make_matrix_calls.pl ${output}_chr${chr}_exome_table.csv ${output} $chr

touch ${output}_snpStats/chr${chr}.done  ##here we mark that the scripts finished

" >> $script

	fi
    done
fi


if [[ "$convertToR" == "yes" ]]; then
    
    if [ ! -e ${output}_snpStats ]; then mkdir ${output}_snpStats; echo "Created ${output}_snpStats"; fi

    for chr in `seq 1 22` X; do
	
	script=cluster/submission/subscript_chr${chr}.sh
	
	if [[ ! -s ${output}_snpStats/chr${chr}_snpStats.RData || "$force" == "yes" ]]; then

	    echo "

$Rbin CMD BATCH --no-save --no-restore --chromosome=${chr} --root=${output} /cluster/project8/vyp/vincent/Software/pipeline/msample_calling/convert_to_R.R cluster/R/convert_to_R_chr${chr}.out

" >> $script

	fi

    done
fi

##############################
for chr in `seq 1 22` X; do
    ls -ltrh cluster/submission/subscript_chr${chr}.sh
done
wc -l $mainScript
