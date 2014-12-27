iFolder=/scratch2/vyp-scratch2/vincent/GATK/HC/IoO
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
	--output )
	    shift
	    output=$1;;
	--genotype )
	    shift
	    genotype=$1;;
	--gVCFlist )
	    shift
	    gVCFlist=$1;;
        -* )
            echo "Unrecognized option: $1"
            exit 1;;
    esac
    shift
    if [ "$#" = "0" ]; then break; fi
done



if [ ! -d $iFolder ]; then 
    echo "$iFolder is not a directory"
    exit
fi





############## Now options are all set
gVCFFolder=${iFolder}/gVCF
output=/scratch2/vyp-scratch2/vincent/GATK/mainset_${currentUCLex}/mainset_${currentUCLex}


if [[ "$genotype" == "yes" ]]; then
    echo "Running the genotype module"

    for chr in `seq 1 22` X; do
	
	echo "testing ${output}_chr${chr}.idx"
	
	if [ ! -e ${output}_chr${chr}.idx ]; then 
	
	    script=cluster/submission/genotype_${chr}.sh

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
	    
	
	    cat $gVCFlist | while read path id; do
		gVCF=${path}/chr${chr}/${id}.gvcf.gz
		if [ ! -s $gVCF ]; then echo "Cannot find $gVCF"; fi

		echo "   --variant $gVCF \\" >> $script
	    done
	    echo "   -o ${output}_chr${chr}" >> $script
	    
	    ls -ltrh $script
	    qsub $script
	fi
    done
fi