

computer=CS
java17=java

if [[ "$computer" == "CS" ]]; then
    Software=/cluster/project8/vyp/vincent/Software
    java17=/share/apps/jdk1.7.0_25/jre/bin/java
    bundle=/scratch2/vyp-scratch2/GATK_bundle
    target=/cluster/project8/vyp/exome_sequencing_multisamples/target_region/data/merged_exome_target_cleaned.bed
    tempFolder=/scratch2/vyp-scratch2/vincent/temp/novoalign
fi




GATK=${Software}/GenomeAnalysisTK-nightly-2014-03-30/GenomeAnalysisTK.jar
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



filtDup=TRUE
clean=TRUE
saveMemory=yes
paired=TRUE
summaryStats=no
fixUnique=no
makegVCF=no

ispool="no"
poolSAM=""
poolnovoPile=""
fasta="default.fasta"
baitFile="default.bait"
extraID=""

align=no


until [ -z "$1" ]; do
	# use a case statement to test vars. we always test $1 and shift at the end of the for block.
    case $1 in
	--inputFiles)
	    shift
	    i=0
	    nfiles=$1
	    inputFiles[ $i ]=$nfiles
	    for f in `seq 1 $nfiles`; do 
		shift
		inputFiles[ $f ]=$1
		echo $f $1
	    done;;
	--extraID )
	    shift
	    extraID="$1_";;
	--filtDup )
	    shift
	    filtDup=$1;;
	--fixUnique )
	    shift
	    fixUnique=$1;;
        --tempFolder )   ##specify a temp directory for the java picard code
          shift
          tempFolder=$1;;
	--output )    ### critical to specify the output file
	    shift
	    output=$1;;
	--tparam )
	    shift
	    tparam=$1;;
	--makegVCF )
	    shift
	    makegVCF=$1;;
	--align)
	    shift
	    align=$1;;
	--projectID)
	    shift
	    projectID=$1;;
	--mainTable)
	    shift
	    mainTable=$1;;
	--format)
	    shift
	    format=$1;;
        --reference)
            shift
            reference=$1;;
	--fasta)
	    shift
	    fasta=$1;;
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


if [[ "$tempFolder" != "" ]]; then
    picardJavaTemp="TMP_DIR=$tempFolder"
fi

echo -e "
Align: ${align} (t parameter is ${tparam})
Annotate: $annotate
Clean non required SAM/BAM files: $clean
Input format: $inputFormat\n\n\n
"


extra="--rOQ --hdrhd 3 -H -k -a -o Soft -t ${tparam}"


########################### choice of reference sequence

#if [ "$type" == RNA ]; then
#    echo "Aligning RNA data"
#    extra="-t ${tparam} -o Soft -a -v 20 0 200 \"[>]([^_]*)_\""
#fi

if [[ "$baitFile" == "default.bait" ]]; then baitFile=${query}.intList; fi
code=`basename $output`

############################### creates basic folders
for folder in cluster cluster/out cluster/error cluster/submission; do
    if [ ! -e $folder ]; then mkdir $folder; fi
done

###############  now let us check that the reference exists
for file in $reference; do
    ls -lh $file
    if [ ! -e "$file" ]; then 
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

if [[ "$summaryStats" == "yes" ]]; then ((nhours=nhours+2)); vmem=3; fi
if [[ "$fixUnique" == "yes" ]]; then ((nhours=nhours+24)); vmem=14; memory=12; fi
if [[ "$align" == "yes" ]]; then ((nhours=nhours+240)); ncores=12; vmem=1.3; memory=3; memory2=7; fi ##10 days?



mainScript=cluster/submission/${projectID}_main.sh
mainTable=cluster/submission/${projectID}_table.tab







########################### Now writing the script

echo "---------------------- Main script $mainScript"

echo "
#!/bin/bash
#$ -S /bin/bash
#$ -o cluster/out
#$ -e cluster/error
#$ -cwd" > $mainScript



if [[ "$computer" == "CS" ]]; then
    
	echo "#$ -pe smp ${ncores}
#$ -l scr=20G
#$ -l tmem=${vmem}G,h_vmem=${vmem}G
#$ -l h_rt=${nhours}:0:0
#$ -tc 20
#$ -V
#$ -R y" >> $mainScript

fi


script=`echo $mainScript | sed -e 's/.sh$//'`_${code}.sh

echo "

##start of script

" > $script


echo $script >> $mainTable


if [[ "$align" == "yes" ]]; then
    
    if [ ! -e ${tempFolder}/${code} ]; then mkdir ${tempFolder}/${code}; fi

    echo "

#$novoalign -c 11 -o SAM $'@RG\tID:${extraID}${code}\tSM:${extraID}${code}\tLB:${extraID}$code\tPL:ILLUMINA' $extra -F ${inputFormat} -f ${inputFiles[ 1 ]} ${inputFiles[ 2 ]}  -d $reference | ${samblaster} -e -d ${output}_disc.sam  | ${samtools} view -Sb - | $novosort - -t ${tempFolder} -c 8 -m ${memory2}G -i -o ${output}_sorted_unique.bam   ###this line is too ambitious, and I am now breaking it down into smaller pieces

$novoalign -c 11 -o SAM $'@RG\tID:${extraID}${code}\tSM:${extraID}${code}\tLB:${extraID}$code\tPL:ILLUMINA' $extra -F ${inputFormat} -f ${inputFiles[ 1 ]} ${inputFiles[ 2 ]}  -d $reference | ${samblaster} -e -d ${output}_disc.sam  | ${samtools} view -Sb - > ${output}.bam

${samtools} view -Sb ${output}_disc.sam | $novosort - -t ${tempFolder} -c 11 -m ${memory2}G -i -o ${output}_disc_sorted.bam

$novosort -t ${tempFolder}/${code} -c 11 -m ${memory2}G -i -o ${output}_sorted_unique.bam ${output}.bam

#rm ${output}_disc.sam ${output}.bam


" >> $script
fi



if [[ "$makegVCF" == "yes" ]]; then
##           --dbsnp ${bundle}/dbsnp_137.b37.vcf \

    echo "
$java17 -Djava.io.tmpdir=${tempFolder} -Xmx7g -jar $GATK -T HaplotypeCaller -R $fasta -I ${output}_sorted_unique.bam  \
       --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 \
       -stand_call_conf 30.0 \
       -stand_emit_conf 10.0 \
       --activeRegionExtension 100 \
       --GVCFGQBands 10 --GVCFGQBands 20 --GVCFGQBands 60 \
       -o ${output}.gvcf.gz
" >> $script

fi

echo "$date" >> $script  ##to measure the duration
echo $script

