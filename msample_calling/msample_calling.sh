#! /bin/bash

# prints to stderr in red
function error() { >&2 echo -e "\033[31m$*\033[0m"; }
function stop() { error "$*"; exit 1; }
try() { "$@" || stop "cannot $*"; }

fasta=/scratch2/vyp-scratch2/reference_datasets/human_reference_sequence/human_g1k_v37.fasta
bundle=/scratch2/vyp-scratch2/reference_datasets/GATK_bundle

java=/share/apps/jdk1.7.0_45/bin/java
tmpDir=/scratch0/vyp
target=/cluster/project8/vyp/exome_sequencing_multisamples/target_region/data/merged_exome_target_cleaned.bed

GATK=/cluster/project8/vyp/vincent/Software/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar

submit=no

force=no
genotype=yes
gVCFlist=none


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


if [[ "$genotype" == "yes" ]]; then
    echo "Running the genotype module"
    
    for chr in `seq 1 22` X; do
	
	echo "testing ${output}_chr${chr}.idx"
	
	if [ ! -e ${output}_chr${chr}.idx ]; then 
	    
	    script=cluster/submission/genotype_chr${chr}.sh
	    
	    echo "
#$ -o cluster/out
#$ -e cluster/error
#$ -S /bin/bash
#$ -l h_vmem=12.9G,tmem=12.9G
#$ -l h_rt=72:0:0
#$ -R y
#$ -pe smp 1
#$ -cwd " > $script
	    
	    echo "

$java -Xmx9g -jar $GATK \\
   -R $fasta \\
   -T GenotypeGVCFs \\
   -L $chr -L $target --interval_set_rule INTERSECTION --interval_padding 100  \\
   --annotation InbreedingCoeff --annotation QualByDepth --annotation HaplotypeScore --annotation MappingQualityRankSumTest --annotation ReadPosRankSumTest --annotation FisherStrand \\
   --dbsnp ${bundle}/dbsnp_137.b37.vcf \\" >> $script
	    
	    
	    while read path id; do
		gVCF=${path}/chr${chr}/${id}.gvcf.gz
		echo "Including $gVCF"
		
		if [ ! -s $gVCF ]; then stop "Cannot find $gVCF"; fi
		if [ ! -s $gVCF.tbi ]; then stop "Cannot find $gVCF.tbi"; fi
		
		echo "   --variant $gVCF \\" >> $script
	    done < <(tail -n +2 $gVCFlist)
	    echo "   -o ${output}_chr${chr}" >> $script
	    
	    ls -ltrh $script
	    qsub $script
	fi
    done 
    
fi