#! /bin/bash
# this script has only been tested with bash and may not work with other shells

# script exits if return value of a command is not zero
set -e
# this forces all variables to be defined
set -u
# for debugging prints out every line before executing it
set -x

# This script does the following three tasks:
# 1) Run alignment to generate BAM and SAM files.
#    Novoalign and samtools are used for this purpose.
# 2) Genotype calling with the GATK HaplotypeCaller.  This generates the gVCF file.
# 3) Jointly call VCFs using the GATK GenotypeGVCFs . This generates a multi-sample combined VCF file.

# prints to stderr in red
function error() { >&2 echo -e "\033[31m$*\033[0m"; }
function stop() { error "$*"; exit 1; }
try() { "$@" || stop "cannot $*"; }


####################### Alignment using Novoalign  ###########################################################################
# The alignment creates the SAM and BAM files for GATK variant calling
function mode_align() {
    outputdir=${projectID}/align
    mainScript=${outputdir}/scripts/align.sh
    output=${projectID}/align/data/
    mkdir -p $output
    ##10 days? Perhaps more.
    nhours=${nhours-240}
    ncores=${ncores-6}
    vmem=${vmem-2.6}
    ##used for the sort function, seem to crash when using 10
    #memory2=5
    memory2=${memory2-7}
    #
    tparam=${tparam:-250}
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
    #tail -n +2 $supportFrame | while read code f1 f2
    while read code f1 f2
    do
        ## do not create the script if .bam.bai exists and is of non-zero size
        ## unless force is set to yes 
        if [[ ! -s ${output}/${code}_sorted_unique.bam.bai ]] || [[ "$force" == "yes" ]]
        then
            echo ${output}/${code}_sorted_unique.bam.bai does not exist
            if [ ! -e $f1 ]; then stop "$f1 does not exist"; fi
            if [ ! -e $f2 ]; then stop "$f2 does not exist"; fi
            # delete contents of temp dir [vplagnol/pipelines/issues/11]
            rm -f ${tempFolder}/${code}/*
            mkdir -p ${tempFolder}/${code} 
cat >${mainScript%.sh}_${code}.sh<<EOL
$novoalign -c ${ncores} -o SAM $'@RG\tID:${extraID}${code}\tSM:${extraID}${code}\tLB:${extraID}$code\tPL:ILLUMINA' --rOQ --hdrhd 3 -H -k -a -o Soft -t ${tparam} -F ${inputFormat} -f ${f1} ${f2}  -d ${novoalignRef} | ${samblaster} -e -d ${output}/${code}_disc.sam  | ${samtools} view -Sb - > ${output}/${code}.bam
${samtools} view -Sb ${output}/${code}_disc.sam | $novosort - -t ${tempFolder}/${code} -c ${ncores} -m ${memory2}G -i -o ${output}/${code}_disc_sorted.bam
$novosort -t ${tempFolder}/${code} -c ${ncores} -m ${memory2}G -i -o ${output}/${code}_sorted_unique.bam ${output}/${code}.bam
rm ${output}/${code}_disc.sam ${output}/${code}.bam
EOL
        else
            #echo ${output}/${code}_sorted_unique.bam.bai already exists
            rm -f ${mainScript%.sh}_${code}.sh
        fi
    done < <(tail -n +2 $supportFrame)
}


####################### GATK HaplotypeCaller genome-wide  ############################################################
# Take as input the sorted, unique per sample BAM files and produces the gVCF files.
# The gVCF files will not be split by chromosome. This is only pratical for smaller
# non-genomewide datasets.
# Input: BAM files.
# Output: per sample gvcf files.
function mode_gvcf_unsplit() {
    input=${projectID}/align/data/
    output=${projectID}/gvcf_unsplit/data/
    mkdir -p $output
    nhours=${nhours-72}
    ncores=${ncores-1}
    vmem=${vmem-8}
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
    #start of while loop
    #each line of the support file is read
    #and a script each is generated
    while read code f1 f2
    do
      if [ ! -s ${output}/${code}.gvcf.gz.tbi ] || [ "$force" == "yes" ]
      then
          echo ${output}/${code}.gvcf.gz.tbi does not exist
          #if file is empty stop
          [ ! -s ${input}/${code}_sorted_unique.bam ] && stop "${input}/${code}_sorted_unique.bam does not exist" 
           #Call SNPs and indels simultaneously via local re-assembly of haplotypes in an active region.
          echo "
          $HaplotypeCaller \
           -R $fasta \
           -I ${input}/${code}_sorted_unique.bam  \
           --emitRefConfidence GVCF \
           --variant_index_type LINEAR \
           --variant_index_parameter 128000 \
           -stand_call_conf 30.0 \
           -stand_emit_conf 10.0 \
           --downsample_to_coverage 250 \
           --GVCFGQBands 10 --GVCFGQBands 20 --GVCFGQBands 60 \
           -o ${output}/${code}.gvcf.gz
            " > ${mainScript%.sh}_${code}.sh
      else
          rm -f ${mainScript%.sh}_${code}.sh
      fi
    done < <(tail -n +2 $supportFrame)
    #end of while loop
}



####################### Bedtools coverageBed of BAM files #####################################################################
# Take as input the sorted, unique per sample BAM files, and the capture bed files and produces the coverage files
# Input: BAM files. BED files
function mode_coverage() {
    input=${projectID}/align/data/
    #output=${projectID}/coverage/data/
    #outputdir=${projectID}/coverage_Macrogen
    outputdir=${projectID}/coverage_BGI
    output=${outputdir}/data/
    mainScript=${outputdir}/scripts/$mode.sh
    nhours=${nhours-4}
    ncores=${ncores-2}
    vmem=${vmem-4}
    mkdir -p $outputdir/data $outputdir/scripts $outputdir/err $outputdir/out
    # Macrogen capture file
    #bed=/goon2/project99/IBDAJE_raw/Macrogen/nochr_SureSelect_v4.bed
    # BGI capture file
    bed=/goon2/project99/IBDAJE_raw/BGI/nochr_target_region.bed
    while read code f1 f2
    do
    bam=${input}/${code}_sorted_unique.bam  
    depth=${output}/${code}_depth.txt
    #if file is empty stop
    #[ ! -s $bam ] && stop "${bam} does not exist" 
    [ ! -s $bam ] && error "${bam} does not exist" 
    if [ ! -s $depth ] || [ "$force" == "yes" ]
    then
cat >${mainScript%.sh}_${code}.sh<<EOL
$coverageBed -abam ${bam} -b ${bed} -d | tr '\t' ',' > ${depth}
for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y
do
    grep "^\$chr," ${depth} > ${depth%_depth.txt}_\${chr}_depth.txt
done
rm ${depth}
EOL
    else
        rm ${mainScript%.sh}_${code}.sh 
    fi
    done < <(tail -n +2 $supportFrame)
}


####################### Picard HsMetrics ######################################################################################
# Take as input the sorted, unique per sample BAM files, and the capture bed files and produces the coverage files
# Input: BAM files. BED files
function mode_hsmetrics() {
    input=${projectID}/align/data/
    #output=${projectID}/coverage/data/
    #outputdir=${projectID}/coverage_Macrogen
    outputdir=${projectID}/hstmetrics
    output=${outputdir}/data/
    mainScript=${outputdir}/scripts/$mode.sh
    nhours=${nhours-4}
    ncores=${ncores-1}
    vmem=${vmem-8}
    mkdir -p $outputdir/data $outputdir/scripts $outputdir/err $outputdir/out
    # Macrogen capture file
    #bed=/goon2/project99/IBDAJE_raw/Macrogen/nochr_SureSelect_v4.bed
    # BGI capture file
    bed=/goon2/project99/IBDAJE_raw/BGI/nochr_target_region.bed
    while read code f1 f2
    do
    bam=${input}/${code}_sorted_unique.bam  
    depth=${output}/${code}_depth.txt
    #if file is empty stop
    #[ ! -s $bam ] && stop "${bam} does not exist" 
    [ ! -s $bam ] && error "${bam} does not exist" 
    if [ ! -s $depth ] || [ "$force" == "yes" ]
    then
cat >${mainScript%.sh}_${code}.sh << EOL
$java -jar $picard_CalculateHsMetric R=${fasta} INPUT=${bam} OUTPUT=${depth} BAIT_INTERVALS=${bed} TARGET_INTERVALS=${bed}
EOL
    else
        rm ${mainScript%.sh}_${code}.sh 
    fi
    done < <(tail -n +2 $supportFrame)
}


####################### vcftools missingness in VCF files #####################################################################
# Take as input the VCF files and reports missingness
# Input: VCF files
function mode_missingness() {
    input=${projectID}/gvcf/data/
    output=${projectID}/missingness/data/
    mkdir -p $output
    nhours=${nhours-6}
    ncores=${ncores-1}
    vmem=${vmem-2}
    # Macrogen capture file
    bed=/cluster/project8/IBDAJE/SureSelect/from_Macrogen/SureSelect_v4.bed
    while read code f1 f2
    do
        ##one job per chromosome to save time
        for chrCode in `seq 1 $cleanChrLen`
        do
            chrCleanCode=${cleanChr[ $chrCode ]}
            #sometimes the dict index contains just the chrom number sometimes
            #it contains chr<number>
            #should check which scenario where are in
            ##if the index is not there, we assume that we have to do the whole job
            if [ ! -s ${output}/${code}_chr${chrCleanCode}.gvcf.gz.tbi ] || [ "$force" == "yes" ]
            then
              echo ${output}/${code}_chr${chrCleanCode}.gvcf.gz.tbi does not exist
              #if file is empty stop
              [ ! -s ${input}/${code}_sorted_unique.bam ] && stop "${input}/${code}_sorted_unique.bam does not exist" 
cat >${mainScript%.sh}_${code}_chr${chrCode}.sh << EOL
    vcftools --bed ${bed} --gzvcf ${gvcf} --minGQ 20 --minDP 10 --recode --out ${out} > log
    vcftools --vcf ${out}.recode.vcf --missing-indv --out Levine_${chr}_missing
EOL
            else
                rm -f ${mainScript%.sh}_${code}_chr${chrCode}.sh
            fi
        done
    done < <(tail -n +2 $supportFrame)
}



####################### GATK HaplotypeCaller split by chromosome  ############################################################
# Take as input the sorted, unique per sample BAM files and produces the gVCF files
# Splits by chromosome.
# Input: BAM files.
# Output: per sample per chromosome gvcf files.
function mode_gvcf() {
    input=${projectID}/align/data/
    outputdir=${projectID}/gvcf
    output=${outputdir}/data
    mkdir -p $outputdir/data $outputdir/out $outputdir/err $outputdir/scripts
    nhours=${nhours-24}
    #ncores=${ncores-1}
    vmem=${vmem-7.8}
    mainScript=${outputdir}/scripts/gvcf.sh
    tc=`tail -n+2 $supportFrame | wc -l`
    #script files get regenerated on every run
    rm -f ${projectID}/gvcf/scripts/*.sh
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
    #start of while loop
    #each line of the support file is read
    #and a script each is generated
    while read code f1 f2
    do
        ##one job per chromosome to save time
        for chrCode in `seq 1 $cleanChrLen`
        do
            chrCleanCode=${cleanChr[ $chrCode ]}
            #sometimes the dict index contains just the chrom number sometimes
            #it contains chr<number>
            #should check which scenario where are in
            ##if the index is not there, we assume that we have to do the whole job
            if [ ! -s ${output}/${code}_chr${chrCleanCode}.gvcf.gz.tbi ] || [ "$force" == "yes" ]
            then
              echo ${output}/${code}_chr${chrCleanCode}.gvcf.gz.tbi does not exist
              #if file is empty stop
              #[ ! -s ${input}/${code}_sorted_unique.bam ] && stop "${input}/${code}_sorted_unique.bam does not exist" 
              [ ! -s ${input}/${code}_sorted_unique.bam ] && error "${input}/${code}_sorted_unique.bam does not exist" 
              #Call SNPs and indels simultaneously via local re-assembly of haplotypes in an active region.
              echo "
              $HaplotypeCaller \
               -R $fasta \
               -I ${input}/${code}_sorted_unique.bam  \
               --emitRefConfidence GVCF \
               --variant_index_type LINEAR \
               --variant_index_parameter 128000 \
               -stand_call_conf 30.0 \
               -stand_emit_conf 10.0 \
               -L ${chrPrefix}${chrCleanCode} \
               --downsample_to_coverage 200 \
               --GVCFGQBands 10 --GVCFGQBands 20 --GVCFGQBands 50 \
               -o ${output}/${code}_chr${chrCleanCode}.gvcf.gz
              " > ${mainScript%.sh}_${code}_chr${chrCode}.sh
            else
                rm -f ${mainScript%.sh}_${code}_chr${chrCode}.sh
            fi
        done
    done < <(tail -n +2 $supportFrame)
}

####################### GATK CombineGVCFs split by chromosome  ############################################################
# Take as input the per sample gVCF files and produces the combined gvcf file.
# Splits by chromosome.
function mode_combinegvcf() {
    nhours=${nhours-12}
    ncores=${ncores-1}
    vmem=${vmem-6}
    input=${projectID}/gvcf/data/
    output=${projectID}/combinegvcf/data/
    mkdir -p $output
    ##one script per chromosome
    for chrCode in `seq 1 $cleanChrLen`
    do
        chrCleanCode=${cleanChr[ $chrCode ]}
        VARIANTS=`echo $input/*_chr${chrCleanCode}.gvcf.gz | xargs -n1 | xargs -I f echo -n ' --variant' f`
        echo "
        $CombineGVCFs \
        -R $fasta \
        -L ${chrPrefix}${chrCleanCode} \
        -o ${output}/chr${chrCleanCode}.gvcf.gz \
        $VARIANTS
        " > ${mainScript%.sh}_chr${chrCode}.sh
    done
}



####################### GATK GenotypeGVCFs  ##################################################################################
### This is the part that combines all the VCFs across samples to do the joint calling.
### This is a more practical aprroach of doing joint-calling than using the UnifiedGenotyper
### which relies on the BAM files.
function mode_jointGenotyping() {
    input=${projectID}/gvcf/data/
    outputdir=${projectID}/jointGenotyping
    output=${outputdir}/data/
    mkdir -p $outputdir/data $outputdir/out $outputdir/err $outputdir/scripts
    nhours=${nhours-12}
    ncores=${ncores-1}
    vmem=${vmem-4}
    mainScript=${outputdir}/scripts/jointGenotyping.sh
    rm -f ${projectID}/jointGenotyping/scripts/*.sh
    for chrCode in `seq 1 $cleanChrLen`
    do 
        ##one job per chromosome to save time
        chrCleanCode=${cleanChr[ $chrCode ]}
        ##if the index is missing, or we use the "force" option
        if [ ! -s ${output}/jointGenotyping_chr${chrCleanCode}.vcf.gz ] || [ "$force" == "yes" ]
        then 
            #Genotypes any number of gVCF files that were produced by the Haplotype Caller into a single joint VCF file.
            echo "
            $java -Xmx2g -jar $GATK \
               -T GenotypeGVCFs \
               -R $fasta \
               -L ${chrPrefix}${chrCleanCode} \
               --interval_padding 100  \
               --annotation InbreedingCoeff \
               --annotation QualByDepth \
               --annotation HaplotypeScore \
               --annotation MappingQualityRankSumTest \
               --annotation ReadPosRankSumTest \
               --annotation FisherStrand \\" > ${mainScript%.sh}_chr${chrCleanCode}.sh
            # for each line in support file
            # the following syntax does not respond to exit
            # because the pipe creates another subprocess
            #tail -n +2 $supportFrame | while read code f1 f2
            while read code f1 f2
            do  ### now look at each gVCF file
                gVCF="${input}/${code}_chr${chrCleanCode}.gvcf.gz"
                if [[ ! -s $gVCF ]]
                then
                    error "Cannot find $gVCF"
                    echo "cat $input/../err/gvcf_${code}_chr${chrCleanCode}_*.err | grep ERROR"
                    cat $input/../err/gvcf_${code}_chr${chrCleanCode}_*.err | grep ERROR
                    exit 1
                else
                    echo "   --variant $gVCF \\" >> ${mainScript%.sh}_chr${chrCleanCode}.sh
                fi
            done < <(tail -n +2 $supportFrame)
            echo "   -o ${output}/jointGenotyping_chr${chrCleanCode}.vcf.gz" >> ${mainScript%.sh}_chr${chrCleanCode}.sh
        else 
            # if the file already exists then delete previous script
            rm -f ${mainScript%.sh}_chr${chrCleanCode}.sh
        fi
    done
}


####################### GATK CombineGVCFs   ##################################################################################
### 
### 
### 
function mode_CombineGVCFs() {
    nhours=${nhours-30}
    ncores=${ncores-1}
    vmem=${vmem-9}
    tc=${tc-80}
    rm -f ${projectID}/CombinedGVCFs/scripts/*.sh
    if [[ ! -f $supportFile ]]
    then
        stop "support file does not exists"
    fi
    while IFS=',' read -r chrCleanCode input
    do
       #echo $chrCode $output
       VARIANTS=`echo --variant $input | sed 's/,/ --variant /g'`
       echo "
       $java -Djava.io.tmpdir=${tempFolder} -Xmx7g -jar $GATK \
       -T CombineGVCFs \
       -R $fasta \
       -L ${chrPrefix}${chrCleanCode} \
       -o ${output}/chr${chrCleanCode}.gvcf.gz \
       $VARIANTS
       " > ${mainScript%.sh}_chr${chrCleanCode}.sh
    done < <(tail -n+2 $supportFile) 
}


####################### VEP annotation      ##################################################################################
### 
### 
### 
function mode_annotation() {
    input=${projectID}/jointGenotyping/data/
    outputdir=${projectID}/annotation/
    output=${outputdir}/data/
    mkdir -p ${output} ${outputdir}/err ${outputdir}/out ${outputdir}/scripts
    nhours=${nhours-12}
    ncores=${ncores-1}
    vmem=${vmem-2}
    mainScript=${outputdir}/scripts/annotation.sh
    rm -f ${projectID}/annotation/scripts/*.sh
    for i in `seq 1 $cleanChrLen`
    do
        chrCode=${cleanChr[ $i ]}
        INPUT=${input}/jointGenotyping_chr${chrCode}.vcf.gz 
        ##if the index is missing, or we use the "force" option
        [[ ! -s ${INPUT} ]] && error "${INPUT} MISSING"
       #echo $chrCode $output
       DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../annotation/
cat >${mainScript%.sh}_chr${chrCode}.sh << EOL
zcat $INPUT | python ${DIR}/multiallele_to_single_gvcf.py > ${output}/chr${chrCode}-single.vcf
bash $DIR/run_VEP.sh --vcfin ${output}/chr${chrCode}-single.vcf --chr $chrCode --reference $reference --vcfout ${output}/VEP_${chrCode}.vcfout
python $DIR/processVEP.py ${output}/VEP_${chrCode}.vcfout 
EOL
   done
}



###
function usage() {
    echo "syntax: $0"
    echo "--mode : [`declare -F | grep mode | sed 's/.*mode_//' | awk 'BEGIN{ORS="|";}{print;}' | sed 's/|$//'`]"
    echo "--extraID : writes to the ReadGroup of the aligner "
    echo "--tempFolder : temp directory for the java picard code"
    echo "--supportFrame : critical to specify the output file"
    echo "--supportFile : for joining gVCF files"
    echo "--tparam : novoalign -t argument "
    echo "--projectID : directory under which all processing happens"
    echo "--reference"
    echo "--force"
    #echo "--enforceStrict"
    echo "--inputFormat : novoalign format [STDFQ]"
    echo "--help : prints this message"
    exit 1
}


############ default values
#parameters to aligner
inputFormat=STDFQ
####
force=no
#enforceStrict=no

# include mithchondrial DNA?  Maybe set an optional parameter
#cleanChr=(targets 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M )
cleanChr=(targets 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y ) 
cleanChrLen=${#cleanChr[@]}
#need to decrement because of header
cleanChrLen=$(( cleanChrLen-1 ))

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
    --supportFile )    
        shift
        supportFile=$1;;
    --aligner-tparam )
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
# @pontikos: what is the point of this option @vplagnol ?
    #--enforceStrict)
        #shift
        #enforceStrict=$1;;
    --inputFormat)
        shift
        inputFormat=$1;;
    --nhours)
        shift
        nhours=$1;;
    --ncores)
        shift
        ncores=$1;;
    --vmem)
        shift
        vmem=$1;;
    --aligner-memory)
        shift
        memory2=$1;;
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
HaplotypeCaller="$java -Djava.io.tmpdir=${tempFolder} -Xmx5g -jar $GATK -T HaplotypeCaller"
CombineGVCFs="$java -Djava.io.tmpdir=${tempFolder} -Xmx4g -jar $GATK -T CombineGVCFs"
novoalign=${Software}/novocraft3/novoalign
novosort=${Software}/novocraft3/novosort
#novoalign=/cluster/project8/vyp/pontikos/Software/novocraft/novoalign
#novosort=/cluster/project8/vyp/pontikos/Software/novocraft/novosort
samblaster=${Software}/samblaster/samblaster
##samtools
samtools=${Software}/samtools-1.1/samtools
##bedtools
coverageBed=${Software}/bedtools-2.17.0/bin/coverageBed
## Picard
picard=${Software}/picard-tools-1.100/
picard_CreateSequenceDictionary=${picard}/CreateSequenceDictionary.jar 
picard_MarkDuplicates=${picard}/MarkDuplicates.jar
picard_CalculateHsMetric=${picard}/CalculateHsMetrics.jar
picard_SamToFastq=${picard}/SamToFastq.jar
## vcflib
# https://github.com/ekg/vcflib
vcfbreakmulti=/cluster/project8/vyp/pontikos/Software/vcflib/bin/vcfbreakmulti

# current default output folder is aligned
oFolder=aligned
# this is the default reference genome
fasta="default.fasta"


########################### supported reference sequences
# all reference sequences and indexes should be under:
# /scratch2/vyp-scratch2/reference_datasets/human_reference_sequence/
fasta=none
novoalignRef=none
if [[ "$reference" == "hg38_noAlt" ]]
then
    fasta=/scratch2/vyp-scratch2/reference_datasets/human_reference_sequence/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
    novoalignRef=/scratch2/vyp-scratch2/reference_datasets/human_reference_sequence/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.k15.s2.novoindex
    chrPrefix='chr'
elif [[ "$reference" == "1kg" ]]
then
    fasta=/scratch2/vyp-scratch2/reference_datasets/human_reference_sequence/human_g1k_v37.fasta
    novoalignRef=/scratch2/vyp-scratch2/reference_datasets/human_reference_sequence/human_g1k_v37.fasta.k15.s2.novoindex
    chrPrefix=''
elif [[ "$reference" == "hg19" ]]
then
    fasta=/scratch2/vyp-scratch2/reference_datasets/human_reference_sequence/hg19_UCSC.fa
    novoalignRef=none
    chrPrefix='chr'
else
    stop Unsupported reference $reference
fi

###############  now let us check that the reference exists
for file in $fasta
do
    ls -lh $file
    if [ ! -e "$file"  ] && [ "$file" != "none" ]
    then 
        stop "Error, reference file $file does not exist"
    fi
done

### Check format of support file.
##should accept tab or space as delimiters
## but does read support tabs and delimeters?
if [[ "$mode" != "CombineGVCFs" ]]
then
    mustBeCode=`head -n1 $supportFrame | cut -f1 -d' ' | cut -f1`  
    mustBeF1=`head -n1 $supportFrame | cut -f2 -d' ' | cut -f2`
    mustBeF2=`head -n1 $supportFrame | cut -f3 -d' ' | cut -f3`
    if [[ "$mustBeCode" != "code" ]]; then stop "The first column of the file $supportFrame must have the name code $mustBeCode"; fi
    if [[ "$mustBeF1" != "f1" ]]; then stop "The second column of the file $supportFrame must have the name f1"; fi
    if [[ "$mustBeF2" != "f2" ]]; then stop "The third column of the file $supportFrame must have the name f2"; fi
else
    supportFrame=$supportFile
fi

############################### creates folders required for qsub and writing logs
mkdir -p $projectID
#symlink the supportFrame to a known location
supportFrameLink=${projectID}/`basename $supportFrame`
cp -f $supportFrame  $supportFrameLink
supportFrame=$supportFrameLink
touch ${projectID}/README
echo "
`date` ${mode}
"
>> ${projectID}/README 

#output dirs get created in the function which is called
#mkdir -p ${projectID}/$mode/data ${projectID}/$mode/scripts ${projectID}/$mode/out ${projectID}/$mode/err 
#output=${projectID}/$mode/data
#mainScript=${projectID}/$mode/scripts/$mode.sh
# call function (see above)
echo mode_${mode}
mode_${mode}

### The script to be submitted to qsub ###
njobs=`ls -1 ${outputdir}/scripts/${mode}_*.sh | wc -l`

echo "
#!/bin/bash
#$ -S /bin/bash
#$ -o /dev/null
#$ -e /dev/null
#$ -cwd
#$ -V
#$ -R y 
#$ -pe smp ${ncores-2}
#$ -l scr=${scratch-1}G
#$ -l tmem=${vmem}G,h_vmem=${vmem}G
#$ -l h_rt=${nhours}:0:0
#$ -t 1-${njobs}
#$ -tc ${tc-25}
set -u
set -x
mkdir -p ${outputdir}/scripts ${outputdir}/data ${outputdir}/err ${outputdir}/out
array=( header \`ls -1 ${outputdir}/scripts/${mode}_*.sh \`) 
script=\${array[ \$SGE_TASK_ID ]} 
scriptname=\`basename \${script%.sh}\`
exec >${outputdir}/out/\${scriptname}_job\${SGE_TASK_ID}_\${JOB_ID}.out 2>${outputdir}/err/\${scriptname}_job\${SGE_TASK_ID}_\${JOB_ID}.err  
echo \$script 
bash \$script

" > $mainScript

echo "Main submission script:"
cat $mainScript
echo
echo
echo run: qsub $mainScript
echo number of jobs in array: `ls -1 ${outputdir}/scripts/${mode}_*.sh | wc -l`




