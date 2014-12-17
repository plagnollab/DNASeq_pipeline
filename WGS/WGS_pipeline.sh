#! /bin/bash
# this script has only been tested with bash and may not work with other shells

# This script does the following three tasks:
# 1) Run alignment to generate BAM and SAM files.
#    Novoalign and samtools are used for this purpose.
# 2) Genotype calling with the GATK HaplotypeCaller.  This generates the gVCF file.
# 3) Jointly call VCFs using the GATK GenotypeGVCFs . This generates a multi-sample combined VCF file.

# prints to stderr in red
function error() { >&2 echo -e "\033[31m$*\033[0m"; }
function stop() { >&2 echo -e "\033[31m$*\033[0m"; exit 1; }

####################### Alignment using Novoalign  ###########################################################################
# The alignment creates the SAM and BAM files for GATK variant calling
function mode_align() {
    ##10 days? Perhaps more.
    nhours=240
    ncores=6
    vmem=2.3
    memory2=7
    for file in $novoalignRef
    do
	ls -lh $file
	if [ ! -e "$file"  ] && [ "$file" != "none" ]
	then 
	    stop "Error, reference file $file does not exist"
	fi
    done
    #start of while loop
    #writes a script for each line of supportFrame
    tail -n +2 $supportFrame | while read code f1 f2
    do
        mkdir -p ${oFolder}/${code}
        output=${oFolder}/${code}/${code}
        script=${mainScript%.sh}_${code}.sh
	echo "
##start of script
" > $script
    ## proceed with that sample if force is set to yes or the output does not exist
	if [[ ! -s ${output}_sorted_unique.bam.bai || "$force" == "yes" ]]
    then
	    if [ ! -e $f1 ]; then stop "$f1 does not exist"; fi
	    if [ ! -e $f2 ]; then stop "$f2 does not exist"; fi
        mkdir -p ${tempFolder}/${code} 
	    echo $script >> $mainTable
	    echo "
$novoalign -c ${ncores} -o SAM $'@RG\tID:${extraID}${code}\tSM:${extraID}${code}\tLB:${extraID}$code\tPL:ILLUMINA' --rOQ --hdrhd 3 -H -k -a -o Soft -t ${tparam} -F ${inputFormat} -f ${f1} ${f2}  -d ${novoalignRef} | ${samblaster} -e -d ${output}_disc.sam  | ${samtools} view -Sb - > ${output}.bam
${samtools} view -Sb ${output}_disc.sam | $novosort - -t ${tempFolder} -c ${ncores} -m ${memory2}G -i -o ${output}_disc_sorted.bam
$novosort -t ${tempFolder}/${code} -c ${ncores} -m ${memory2}G -i -o ${output}_sorted_unique.bam ${output}.bam
#rm ${output}_disc.sam ${output}.bam
"  >> $script
	    echo "$date" >> $script  ##to measure the duration
	    echo $script
	fi
    done
    #end of while loop
}


####################### GATK HaplotypeCaller no splitting by chromosome  #####################################################
# Instead of splitting by chromosome creates a single VCF file.
# This is only manageable for exon and generally smaller sequence datasets.
# This is currently not fully suported in the next step of doing joint-calling.
function mode_singlegvcf() {
    nhours=24
    vmem=6
    #start of while loop
    #each line of the support file is read
    #and a script each is generated
    tail -n +2 $supportFrame | while read code f1 f2
      do
        output=${oFolder}/${code}/${code}
        ## one job per chromosome to save time
        ## if the index is not there, we assume that we have to do the whole job
        if [ ! -s ${output}.gvcf.gz.tbi | "$force" == "yes" ]
        then
            script=${mainScript%.sh}_${code}.sh
            echo $script >> $mainTable
           #Call SNPs and indels simultaneously via local re-assembly of haplotypes in an active region.
            echo "
           $java -Djava.io.tmpdir=${tempFolder} -Xmx4g -jar $GATK \
           -T HaplotypeCaller -R $fasta -I ${output}_sorted_unique.bam  \
           --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 \
           -stand_call_conf 30.0 \
           -stand_emit_conf 10.0 \
           --downsample_to_coverage 200 \
           --GVCFGQBands 10 --GVCFGQBands 20 --GVCFGQBands 50 \
           -o ${output}.gvcf.gz
            " > $script
        fi
    done
}

####################### GATK HaplotypeCaller split by chromosome  ############################################################
# Take as input the sorted, unique BAM files and produces the gVCF files
# Splits by chromosome.
function mode_gvcf() {
    nhours=24
    vmem=6
    #memory2=6
    # GATK_HaplotypeCaller requires a sequence dictionary
    # Maybe the following should be submitted as interactive long job?
    #
    #[[ -e ${fasta%.fasta}.dict ]] && $java -jar $picard_CreateSequenceDictionary R=$fasta O=${fasta%.fasta}.dict
    if [ ! -e ${fasta%.fasta}.dict ]
    then 
        echo $java -jar $picard_CreateSequenceDictionary R=$fasta O=${fasta%.fasta}.dict
        $java -jar $picard_CreateSequenceDictionary R=$fasta O=${fasta%.fasta}.dict
    fi
    #same story
    if [ ! -e ${fasta}.fai ]
    then 
        echo $samtools faidx $fasta
        $samtools faidx $fasta
        file ${fasta}.fai
    fi
    mainScript=cluster/submission/makegVCF.sh
    mainTable=cluster/submission/makegVCF_table.sh
    cleanChr=(targets 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M )
    #start of while loop
    #each line of the support file is read
    #and a script each is generated
    tail -n +2 $supportFrame | while read code f1 f2
    do
	output=${oFolder}/${code}/${code}
	##one job per chromosome to save time
	for chrCode in `seq 1 25`
    do
	    chrCleanCode=${cleanChr[ $chrCode ]}
        #sometimes the dict index contains just the chrom number sometimes
        #it contains chr<number>
        #should check which scenario where are in
	    ##if the index is not there, we assume that we have to do the whole job
	    if [ ! -s ${output}_chr${chrCleanCode}.gvcf.gz.tbi ] || [ "$force" == "yes" ]
        then
           script=${mainScript%.sh}_chr${chrCode}_${code}.sh
           echo $script >> $mainTable
           #Call SNPs and indels simultaneously via local re-assembly of haplotypes in an active region.
           echo "
           $java -Djava.io.tmpdir=${tempFolder} -Xmx4g -jar $GATK \
           -T HaplotypeCaller \
           -R $fasta \
           -I ${output}_sorted_unique.bam  \
           --emitRefConfidence GVCF \
           --variant_index_type LINEAR \
           --variant_index_parameter 128000 \
           -stand_call_conf 30.0 \
           -stand_emit_conf 10.0 \
           -L ${chrCleanCode} \
           --downsample_to_coverage 200 \
           --GVCFGQBands 10 --GVCFGQBands 20 --GVCFGQBands 50 \
           -o ${output}_chr${chrCleanCode}.gvcf.gz
            " > $script
	    fi
	done
    done
    #end of while loop
}

####################### GATK CombineGVCFs split by chromosome  ############################################################
# Take as input the sorted, unique BAM files and produces the gVCF files
# Splits by chromosome.
function mode_combinegvcfs() {
    nhours=12
    vmem=4
    # GATK_HaplotypeCaller requires a sequence dictionary
    # Maybe the following should be submitted as interactive long job?
    #
    #[[ -e ${fasta%.fasta}.dict ]] && $java -jar $picard_CreateSequenceDictionary R=$fasta O=${fasta%.fasta}.dict
    if [ ! -e ${fasta%.fasta}.dict ]
    then 
        echo $java -jar $picard_CreateSequenceDictionary R=$fasta O=${fasta%.fasta}.dict
        $java -jar $picard_CreateSequenceDictionary R=$fasta O=${fasta%.fasta}.dict
    fi
    #same story
    if [ ! -e ${fasta}.fai ]
    then 
        echo $samtools faidx $fasta
        $samtools faidx $fasta
        file ${fasta}.fai
    fi
    cleanChr=(targets 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M )
	output=${oFolder}/${code}/${code}
	##one script per chromosome
	for chrCode in `seq 1 25`
    do
	    chrCleanCode=${cleanChr[ $chrCode ]}
        #sometimes the dict index contains just the chrom number sometimes
        #it contains chr<number>
        #should check which scenario where are in
	    ##if the index is not there, we assume that we have to do the whole job
        script=${mainScript%.sh}_chr${chrCode}.sh
        echo $script >> $mainTable
        VARIANTS=`echo aligned/*/*_chr${chrCode}.gvcf.gz | xargs -n1 | xargs -I f echo -n ' --variant' f`
        #Call SNPs and indels simultaneously via local re-assembly of haplotypes in an active region.
        echo "
        $java -Djava.io.tmpdir=${tempFolder} -Xmx4g -jar $GATK \
        -T CombineGVCFs \
        -R $fasta \
        -L ${chrCleanCode} \
        -o combined/chr${chr}.gvcf.gz \
        $VARIANTS
        " > $script
	done
}



####################### GATK GenotypeGVCFs  ##################################################################################
### This is the part that combines all the VCFs across samples to do the joint calling.
### This is a more practical aprroach of doing joint-calling than using the UnifiedGenotyper
### which relies on the BAM files.
function mode_jointvcf() {
    nhours=12
    vmem=4
    #memory2=noneed
    echo "listScripts" > $mainTable
    cleanChr=(targets 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M )
    for chrCode in `seq 1 25`
    do 
        ##one job per chromosome to save time
        chrCleanCode=${cleanChr[ $chrCode ]}
        output=${oFolder}/combined/combined_chr${chrCleanCode}.vcf.gz
        ##if the index is missing, or we use the "force" option
        if [ ! -s ${output}.tbi ] || [ "$force" == "yes" ]
        then 
            script=`echo $mainScript | sed -e 's/.sh$//'`_chr${chrCleanCode}.sh	
            echo "$script" >> $mainTable
            #Genotypes any number of gVCF files that were produced by the Haplotype Caller into a single joint VCF file.
            echo "
            $java -Xmx2g -jar $GATK \
               -T GenotypeGVCFs \
               -R $fasta \
               -L chr${chrCleanCode} \
               --interval_padding 100  \
               --annotation InbreedingCoeff \
               --annotation QualByDepth \
               --annotation HaplotypeScore \
               --annotation MappingQualityRankSumTest \
               --annotation ReadPosRankSumTest \
               --annotation FisherStrand \\" > $script
            # for each line in support file
            tail -n +2 $supportFrame | while read code f1 f2
            do  ### now look at each gVCF file
            output=${oFolder}/${code}/${code}
            gVCF="${output}_chr${chrCleanCode}.gvcf.gz"
            if [[ "$enforceStrict" == "yes" && ! -s $gVCF ]]
            then
                stop "Cannot find $gVCF"
            fi
            if [ -s $gVCF ]
            then 
                echo "   --variant $gVCF \\" >> $script; 
            fi
            done
            #echo "   --dbsnp ${bundle}/dbsnp_137.b37.vcf \\
            #-o ${oFolder}/combined/combined_chr${chrCleanCode}.vcf.gz" >> $script
            echo "   -o ${oFolder}/combined/combined_chr${chrCleanCode}.vcf.gz" >> $script
      fi
    done
}


###
function usage() {
    echo "syntax: $0"
    echo "--mode : [`declare -F | grep mode | sed 's/.*mode_//' | xargs -n1 echo -n ' |' | sed 's/ \|//'` ]"
    echo "--extraID "
    echo "--tempFolder : temp directory for the java picard code"
    echo "--supportFrame : critical to specify the output file"
    echo "--tparam "
    echo "--projectID"
    echo "--reference"
    echo "--force"
    echo "--enforceStrict"
    echo "--inputFormat"
    echo "--help : prints this message"
    exit 1
}


##
until [ -z "$1" ]
do
    # use a case statement to test vars. we always test $1 and shift at the end of the for block.
    case $1 in
    --extraID )
        shift
        extraID="$1_";;
     --tempFolder )   ##specify a temp directory for the java picard code
        shift
        tempFolder=$1;;
    --supportFrame )    ### critical to specify the output file
        shift
        supportFrame=$1;;
    --tparam )
        shift
        tparam=$1;;
    # the main 3 steps of the program: align or gvcf or jointvcf
    --mode)
        shift
        mode=$1;;
    --projectID)
        shift
        projectID=$1;;
    --reference)
        shift
        reference=$1;;
    --force)
        shift
        force=$1;;
    --enforceStrict)
        shift
        enforceStrict=$1;;
    --inputFormat)
        shift
        inputFormat=$1;;
    -* )
        echo "Unrecognized option: $1"
        usage
        exit 1;;
    esac
    shift
    if [ "$#" = "0" ]; then break; fi
done 


#################### the code below is generic to all modules: compute the nb of jobs and create the final submission array

# If you run this in a different environment than the UCL CS cluster
# then you need to set the env variables.
computer=CS
if [[ "$computer" == "CS" ]]
then
    Software=/cluster/project8/vyp/vincent/Software
    java=/share/apps/jdk1.7.0_45/jre/bin/java
    bundle=/scratch2/vyp-scratch2/GATK_bundle
    target=/cluster/project8/vyp/exome_sequencing_multisamples/target_region/data/merged_exome_target_cleaned.bed
    tempFolder=/scratch2/vyp-scratch2/vincent/temp/novoalign
fi

### Tools needed by this script
# Two functions of GATK will be used HaplotypeCaller and GenotypeGVCFs 
GATK=${Software}/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar
#$java -Djava.io.tmpdir=${tempFolder} -Xmx4g -jar ${GATK}
novoalign=${Software}/novocraft3/novoalign
novosort=${Software}/novocraft3/novosort
samblaster=${Software}/samblaster/samblaster
##samtools
samtools=${Software}/samtools-1.1/samtools
## Picard
picard=${Software}/picard-tools-1.100/
picard_CreateSequenceDictionary=${picard}/CreateSequenceDictionary.jar 
picard_MarkDuplicates=${picard}/MarkDuplicates.jar
picard_CalculateHsMetric=${picard}/CalculateHsMetrics.jar
picard_SamToFastq=${picard}/SamToFastq.jar

############ default values
#parameters to aligner
inputFormat=STDFQ
tparam=250

####
force=no
enforceStrict=no

# current default output folder is aligned
oFolder=aligned
# this is the default reference genome
fasta="default.fasta"
# this is used by the aligner but not sure what it does.  Vincent?
extraID=""



########################### supported reference sequences
# all reference sequences and indexes should be under:
# /scratch2/vyp-scratch2/reference_datasets/human_reference_sequence/
fasta=none
novoalignRef=none
if [[ "$reference" == "hg38_noAlt" ]]
then
    fasta=/scratch2/vyp-scratch2/reference_datasets/human_reference_sequence/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
    novoalignRef=/scratch2/vyp-scratch2/reference_datasets/human_reference_sequence/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.k15.s2.novoindex
elif [[ "$reference" == "1kg" ]]
then
    fasta=/scratch2/vyp-scratch2/reference_datasets/human_reference_sequence/human_g1k_v37.fasta
    novoalignRef=/scratch2/vyp-scratch2/reference_datasets/human_reference_sequence/human_g1k_v37.fasta.k15.s2.novoindex
elif [[ "$reference" == "hg19" ]]
then
    fasta=/scratch2/vyp-scratch2/reference_datasets/human_reference_sequence/hg19_UCSC.fa
    novoalignRef=none
else
    stop Unsupported reference $reference
fi

############################### creates folders required for qsub and writing logs
mkdir -p cluster cluster/out cluster/err cluster/submission

###############  now let us check that the reference exists
for file in $fasta
do
    ls -lh $file
    if [ ! -e "$file"  ] && [ "$file" != "none" ]
    then 
        stop "Error, reference file $file does not exist"
    fi
done


###########################################################
nhours=0
ncores=1
vmem=1
memory2=5  ##used for the sort function, seem to crash when using 10
scratch=0


### Check format of support file.
##should accept tab or space as delimiters
## but does read support tabs and delimeters?
mustBeCode=`head -1 $supportFrame | cut -f1 -d' ' | cut -f1`  
mustBeF1=`head -1 $supportFrame | cut -f2 -d' ' | cut -f2`
mustBeF2=`head -1 $supportFrame | cut -f3 -d' ' | cut -f3`
if [[ "$mustBeCode" != "code" ]]; then stop "The first column of the file $supportFrame must have the name code $mustBeCode"; fi
if [[ "$mustBeF1" != "f1" ]]; then stop "The second column of the file $supportFrame must have the name f1"; fi
if [[ "$mustBeF2" != "f2" ]]; then stop "The third column of the file $supportFrame must have the name f2"; fi

# everything will be stored under $mode
mkdir -p $mode/submission $mode/err $mode/out
mainScript=$mode/submission/run_${mode}.sh
mainTable=$mode/submission/run_${mode}_table.sh

# call function (see above)
echo mode_${mode}
mode_${mode}

### The script to be submitted to qsub ###
echo $mainTable
echo $mainTable.2
cat $mainTable | sort -u > $mainTable.2
mv $mainTable.2 $mainTable
njobs=`cat $mainTable | wc -l $mainTable`

echo "
#!/bin/bash
#$ -S /bin/bash
#$ -o /dev/null
#$ -e /dev/null
#$ -cwd
#$ -V
#$ -R y 
#$ -pe smp ${ncores}
#$ -l scr=${scratch}G
#$ -l tmem=${vmem}G,h_vmem=${vmem}G
#$ -l h_rt=${nhours}:0:0
#$ -t 1-${njobs}
#$ -tc 25
array=( \`cat \"${mainTable}\" \`) 
script=\${array[ \$SGE_TASK_ID ]} 
exec >${mode}/out/${script%.sh}.out 2>${mode}/err/${script%.sh}.err 
echo \$script 
sh \$script
" > $mainScript

echo "Main submission scripts and tables for the align module:"
wc -l $mainScript $mainTable
echo run: qsub $mainScript



