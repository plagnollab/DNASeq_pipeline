###
###

computer=CS
java17=java

if [[ "$computer" == "CS" ]]
then
    Software=/cluster/project8/vyp/vincent/Software
    java17=/share/apps/jdk1.7.0_45/jre/bin/java
    bundle=/scratch2/vyp-scratch2/GATK_bundle
    target=/cluster/project8/vyp/exome_sequencing_multisamples/target_region/data/merged_exome_target_cleaned.bed
    tempFolder=/scratch2/vyp-scratch2/vincent/temp/novoalign
fi

# Two functions of GATK will be used HaplotypeCaller and GenotypeGVCFs 
GATK=${Software}/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar
novoalign=${Software}/novocraft3/novoalign
novosort=${Software}/novocraft3/novosort
samblaster=${Software}/samblaster/samblaster

##samtools
samtools=${Software}/samtools-1.1/samtools

## Picard
picardDup=${Software}/picard-tools-1.100/MarkDuplicates.jar
picardMetrics=${Software}/picard-tools-1.100/CalculateHsMetrics.jar
picardSamToFastq=${Software}/picard-tools-1.100/SamToFastq.jar

############ default values
inputFormat=STDFQ
tparam=250


####
force=no
enforceStrict=no


##### list of action parameters below
align=no
makegVCF=no
makeVCF=no


iFolder=""
oFolder=aligned
fasta="default.fasta"
extraID=""


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
	--makegVCF )
	    shift
	    makegVCF=$1;;
	--makeVCF )
	    shift
	    makeVCF=$1;;
	--align)
	    shift
	    align=$1;;
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
	    exit 1;;
    esac
    shift
    if [ "$#" = "0" ]; then break; fi
done 



extra="--rOQ --hdrhd 3 -H -k -a -o Soft -t ${tparam}"

########################### choice of reference sequence
fasta=none
novoalignRef=none
if [[ "$reference" == "hg38_noAlt" ]]; then
    fasta=/scratch2/vyp-scratch2/reference_datasets/human_reference_sequence/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
    novoalignRef=/scratch2/vyp-scratch2/reference_datasets/human_reference_sequence/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.k15.s2.novoindex
fi

if [[ "$reference" == "1kg" ]]; then
    fasta=/cluster/project8/vyp/vincent/data/reference_genomes/human_g1k_v37.fasta
    novoalignRef=/scratch2/vyp-scratch2/reference_datasets/human_reference_sequence/human_g1k_v37.fasta.k15.s2.novoindex
fi

if [[ "$reference" == "hg19" ]]; then
    fasta=/scratch2/vyp-scratch2/reference_datasets/human_reference_sequence/hg19_UCSC.fa
    novoalignRef=none
fi


############################### creates folders required for qsub and writing logs
mkdir -p cluster cluster/out cluster/err cluster/submission
#for folder in cluster cluster/out cluster/error cluster/submission; do
    #if [ ! -e $folder ]; then mkdir $folder; fi
#done

###############  now let us check that the reference exists
for file in $fasta; do
    ls -lh $file
    if [ ! -e "$file"  ] && [ "$file" != "none" ]
    then 
        echo "Error, reference file $file does not exist"
	exit
    fi
done



###########################################################
nhours=0
ncores=1
vmem=1
memory=2
memory2=5  ##used for the sort function, seem to crash when using 10
queue=queue6
scratch=0

if [[ "$summaryStats" == "yes" ]]; then ((nhours=nhours+2)); vmem=3; fi
if [[ "$makeVCF" == "yes" ]]; then ((nhours=nhours+12)); vmem=6; memory=6; fi
if [[ "$makegVCF" == "yes" ]]; then ((nhours=nhours+24)); vmem=6; memory=6; fi
if [[ "$align" == "yes" ]]; then ((nhours=nhours+240)); ncores=12; vmem=1.3; memory=3; memory2=7; fi ##10 days?


### running checks
mustBeCode=`head -1 $supportFrame | cut -f1 -d' ' | cut -f1`  ##should accept tab or space as delimiters
mustBeF1=`head -1 $supportFrame | cut -f2 -d' ' | cut -f2`
mustBeF2=`head -1 $supportFrame | cut -f3 -d' ' | cut -f3`

if [[ "$mustBeCode" != "code" ]]; then echo "The first column of the file $supportFrame must have the name code $mustBeCode"; exit; fi
if [[ "$mustBeF1" != "f1" ]]; then echo "The second column of the file $supportFrame must have the name f1"; exit; fi
if [[ "$mustBeF2" != "f2" ]]; then echo "The third column of the file $supportFrame must have the name f2"; exit; fi



########################### Now writing the script

if [[ "$align" == "yes" ]]
    then
    
    for file in $novoalignRef; do
	ls -lh $file
	if [ ! -e "$file"  ] && [ "$file" != "none" ]
	    then 
	    echo "Error, reference file $file does not exist"
	    exit
	fi
    done


    mainScript=cluster/submission/align.sh
    mainTable=cluster/submission/align_table.sh
    echo "listScripts" > $mainTable
    #start of while loop
    tail -n +2 $supportFrame | while read code f1 f2
    do
        mkdir -p ${oFolder}/${code}
        #if [ ! -e ${oFolder}/${code} ]; then mkdir ${oFolder}/${code}; fi
        output=${oFolder}/${code}/${code}
        script=`echo $mainScript | sed -e 's/.sh$//'`_${code}.sh
	echo "
##start of script

" > $script
	if [[ "$iFolder" != "" ]]
    then
	    f1=${iFolder}/${f1}
	    f2=${iFolder}/${f2}
	fi
    ## proceed with that sample if force is set to yes or the output does not exist
	if [[ ! -s ${output}_sorted_unique.bam.bai || "$force" == "yes" ]]
    then
	    if [ ! -e $f1 ]; then echo "$f1 does not exist"; exit; fi
	    if [ ! -e $f2 ]; then echo "$f2 does not exist"; exit; fi
	    if [ ! -e ${tempFolder}/${code} ]; then mkdir ${tempFolder}/${code}; fi
	    echo $script >> $mainTable
	    echo "
#$novoalign -c 11 -o SAM $'@RG\tID:${extraID}${code}\tSM:${extraID}${code}\tLB:${extraID}$code\tPL:ILLUMINA' $extra -F ${inputFormat} -f ${f1} ${f2}  -d ${novoalignRef} | ${samblaster} -e -d ${output}_disc.sam  | ${samtools} view -Sb - | $novosort - -t ${tempFolder} -c 8 -m ${memory2}G -i -o ${output}_sorted_unique.bam   ###this line is too ambitious, and I am now breaking it down into smaller pieces

$novoalign -c 11 -o SAM $'@RG\tID:${extraID}${code}\tSM:${extraID}${code}\tLB:${extraID}$code\tPL:ILLUMINA' $extra -F ${inputFormat} -f ${f1} ${f2}  -d ${novoalignRef} | ${samblaster} -e -d ${output}_disc.sam  | ${samtools} view -Sb - > ${output}.bam

${samtools} view -Sb ${output}_disc.sam | $novosort - -t ${tempFolder} -c 11 -m ${memory2}G -i -o ${output}_disc_sorted.bam

$novosort -t ${tempFolder}/${code} -c 11 -m ${memory2}G -i -o ${output}_sorted_unique.bam ${output}.bam

#rm ${output}_disc.sam ${output}.bam
"  >> $script
	    echo "$date" >> $script  ##to measure the duration
	    echo $script
	fi
    done
    #end of while loop

    #### compute the nb of jobs 
    njobs=`wc -l $mainTable | cut -f1 -d' '`
    ((njobs=njobs-1))

    ####### And now we work on the main script that combines them all
    echo "
#!/bin/bash
#$ -S /bin/bash
#$ -o cluster/out
#$ -e cluster/error
#$ -cwd
#$ -pe smp ${ncores}
#$ -l scr=${scratch}G
#$ -l tmem=${vmem}G,h_vmem=${vmem}G
#$ -l h_rt=${nhours}:0:0
#$ -tc 20
#$ -t 1-${njobs}
#$ -V
#$ -R y

array=( \`cat \"${mainTable}\" \`)

script=\${array[ \$SGE_TASK_ID ]}

echo \$script

sh \$script

" > $mainScript

    echo "Main submission scripts and tables for the align module:"
    wc -l $mainScript $mainTable
fi



# Take as input the sorted, unique BAM files and produces the gVCF files
if [[ "$makegVCF" == "yes" ]]
then

    mainScript=cluster/submission/makegVCF.sh
    mainTable=cluster/submission/makegVCF_table.sh

    cleanChr=(targets 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M )
    #sart of while loop
    tail -n +2 $supportFrame | while read code f1 f2
    do
	output=${oFolder}/${code}/${code}
	##one job per chromosome to save time
	for chrCode in `seq 1 25`
    do
	    chrCleanCode=${cleanChr[ $chrCode ]}
	    ##if the index is not there, we assume that we have to do the whole job
	    if [ ! -s ${output}_chr${chrCleanCode}.gvcf.gz.tbi | "$force" == "yes" ]
        then
           script=`echo $mainScript | sed -e 's/.sh$//'`_chr${chrCode}_${code}.sh
           echo $script >> $mainTable
           #Call SNPs and indels simultaneously via local re-assembly of haplotypes in an active region.
           echo "
           $java17 -Djava.io.tmpdir=${tempFolder} -Xmx4g -jar $GATK -T HaplotypeCaller -R $fasta -I ${output}_sorted_unique.bam  \
           --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 \
           -stand_call_conf 30.0 \
           -stand_emit_conf 10.0 \
           -L chr${chrCleanCode} \
           --downsample_to_coverage 200 \
           --GVCFGQBands 10 --GVCFGQBands 20 --GVCFGQBands 50 \
           -o ${output}_chr${chrCleanCode}.gvcf.gz
            " > $script
	    fi
	done
    done
    #end of while loop
fi


##############################################################################################################################

if [[ "$makeVCF" == "yes" ]]; then
    
    mainScript=cluster/submission/makeVCF.sh
    mainTable=cluster/submission/makeVCF_table.sh
    echo "listScripts" > $mainTable
    cleanChr=(targets 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M )
    for chrCode in `seq 1 25`
    do 
        ##one job per chromosome to save time
        chrCleanCode=${cleanChr[ $chrCode ]}
        output=${oFolder}/combined/combined_chr${chrCleanCode}.vcf.gz
        ##if the index is missing, or we use the "force" option
        if [ ! -s ${output}.tbi | "$force" == "yes" ]
        then 
	    script=`echo $mainScript | sed -e 's/.sh$//'`_chr${chrCleanCode}.sh	
	    echo "$script" >> $mainTable
        #Genotypes any number of gVCF files that were produced by the Haplotype Caller into a single joint VCF file.
	    echo "
$java17 -Xmx2g -jar $GATK \\
   -R $fasta \\
   -T GenotypeGVCFs \\
   -L chr${chrCleanCode}  --interval_padding 100  \\
   --annotation InbreedingCoeff --annotation QualByDepth --annotation HaplotypeScore --annotation MappingQualityRankSumTest --annotation ReadPosRankSumTest --annotation FisherStrand \\" > $script
	    tail -n +2 $supportFrame | while read code f1 f2
        do  ### now look at each gVCF file
            output=${oFolder}/${code}/${code}
            gVCF="${output}_chr${chrCleanCode}.gvcf.gz"
		if [[ "$enforceStrict" == "yes" && ! -s $gVCF ]]; then echo "Cannot find $gVCF"; exit; fi
		if [ -s $gVCF ]; then 
		    echo "   --variant $gVCF \\" >> $script; 
		fi
	    done
	    
	#echo "   --dbsnp ${bundle}/dbsnp_137.b37.vcf \\
   #-o ${oFolder}/combined/combined_chr${chrCleanCode}.vcf.gz" >> $script
	    
	    echo "   -o ${oFolder}/combined/combined_chr${chrCleanCode}.vcf.gz" >> $script
	fi    
    done

    #### compute the nb of jobs 
    njobs=`wc -l $mainTable | cut -f1 -d' '`
    ((njobs=njobs-1))


    echo "
#!/bin/bash
#$ -S /bin/bash
#$ -o cluster/out
#$ -e cluster/error
#$ -cwd
#$ -pe smp 1
#$ -l scr=${scratch}G
#$ -l tmem=4G,h_vmem=4G
#$ -l h_rt=24:0:0
#$ -tc 25
#$ -t 1-25
#$ -V
#$ -R y

array=( \`cat \"${mainTable}\" \`)

script=\${array[ \$SGE_TASK_ID ]}

echo \$script

sh \$script

" > $mainScript


    echo "Main submission scripts and tables for the align module:"
    wc -l $mainScript $mainTable
fi


