

computer=CS
java17=java

if [[ "$computer" == "CS" ]]; then
    Software=/cluster/project8/vyp/vincent/Software
    javaTemp="/scratch2/vyp-scratch2/vincent/java_temp"
    java17=/share/apps/jdk1.7.0_25/jre/bin/java
    bundle=/scratch2/vyp-scratch2/GATK_bundle
    target=/cluster/project8/vyp/exome_sequencing_multisamples/target_region/data/merged_exome_target_cleaned.bed
fi

if [[ "$computer" == "UGI" ]]; then
    Software=/ugi/home/shared/vincent/Software
fi

if [[ "$computer" == "vanHeel" ]]; then
    Software=/data_n2/vplagnol/Software
fi


if [[ "$computer" == "ATLAS" ]]; then
   Software=/illumina/pipeline/vincent/Software
fi



GATK=${Software}/GenomeAnalysisTK-nightly-2014-03-30/GenomeAnalysisTK.jar
novoalign=${Software}/novocraft/novoalign
novosort=${Software}/novocraft/novosort
scramble=${Software}/io_lib-1.13.3/progs/scramble


##samtools
samtools=${Software}/samtools-0.1.19/samtools
samtools1=${Software}/samtools-1.1/samtools
vcftools=${Software}/vcftools_0.1.8/cpp/vcftools

tabix=${Software}/tabix-0.2.3/tabix
bgzip=${Software}/tabix-0.2.3/bgzip

## Picard
picardDup=${Software}/picard-tools-1.100/MarkDuplicates.jar
picardMetrics=${Software}/picard-tools-1.100/CalculateHsMetrics.jar
picardSamToFastq=${Software}/picard-tools-1.100/SamToFastq.jar

############ default values
suffix=_sequence.txt.gz
inputFormat=ILMFQ
tparam=250

filtDup=TRUE
clean=TRUE
saveMemory=yes
paired=TRUE
fullPileup=no
summaryStats=no
fixUnique=no
reduceReads=no
makegVCF=no

PCRamplicons=FALSE
haploplex=FALSE

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
	--fullPileup )
	    fullPileup="yes";;	
	--extraID )
	    shift
	    extraID="$1_";;
	--PCRamplicons)
	    PCRamplicons=TRUE;;
	--haploplex)
	    haploplex=TRUE
	    PCRamplicons=TRUE;;
	--filtDup )
	    shift
	    filtDup=$1;;
	--fixUnique )
	    shift
	    fixUnique=$1;;
        --javaTemp )   ##specify a temp directory for the java picard code
          shift
          javaTemp=$1;;
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
	--summaryStats)
	    shift
	    summaryStats=$1;;
	--projectID)
	    shift
	    projectID=$1;;
	--mainTable)
	    shift
	    mainTable=$1;;
	--reduceReads)
	    shift
	    reduceReads=$1;;
	--format)
	    shift
	    format=$1;;
        --reference)
            shift
            reference=$1;;
	--fasta)
	    shift
	    fasta=$1;;
	--baitFile)
	    shift
	    baitFile=$1;;
	--noclean)
	    clean=FALSE;;
	--notpaired)
	    paired=FALSE;;
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


if [[ "$javaTemp" != "" ]]; then
    picardJavaTemp="TMP_DIR=$javaTemp"
fi

echo -e "
Align: ${align} (t parameter is ${tparam})
Annotate: $annotate
Clean non required SAM/BAM files: $clean
Input format: $inputFormat\n\n\n
"


extra="--rOQ --hdrhd 3 -H -k -a -o Soft -t ${tparam}"


##The extra set of parameters for novoalign
if [[ "$PCRamplicons" == "TRUE" ]]; then filtDup=FALSE; fi
if [[ "$haploplex" == "TRUE" ]]; then 
    P1=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
    P2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
    extra="--rOQ --hdrhd 3 -H -k -a $P1 $P2 -o Soft -t ${tparam}"
fi


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
if [[ "$pileup" == "yes" ]]; then ((nhours=nhours+24)); vmem=4; fi
if [[ "$reduceReads" == "yes" ]]; then ((nhours=nhours+18)); vmem=14; memory=12; fi
if [[ "$fixUnique" == "yes" ]]; then ((nhours=nhours+24)); vmem=14; memory=12; fi
if [[ "$align" == "yes" ]]; then ((nhours=nhours+36)); ncores=12; vmem=1.3; memory=3; fi



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



if [[ "$computer" == "UGI" ]]; then
    echo "                                                                                                                                                                                                      
#$ -q $queue

" >> $mainScript
fi


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
export PERL5LIB=${Software}/vcftools_0.1.8/lib:${PERL5LIB}
" > $script


echo $script >> $mainTable


if [[ "$align" == "yes" ]]; then

    nfiles=${inputFiles[0]}

    if [[ "$paired" == "TRUE" ]]; then 
	((npairs=nfiles/2))
	for i in `seq 1 $nfiles`; do
	    ls -ltrh ${inputFiles[ $i ]}
	    
	    if [ ! -e "${inputFiles[ $i ]}" ]; then 
		echo "Error, file ${inputFiles[ $i ]} does not exist"
		exit
	    fi	
	done
    fi
    
    for file in ${fasta}i; do
	if [ ! -e "${file}" ]; then echo "Error, file ${file} does not exist"; exit; fi
    done
    
    if [[ "$inputFormat" == "BAMPE" ]]; then

	echo "

$novosort -n -f -t $javaTemp -c $ncores -m ${memory}G ${inputFiles[ 1 ]} -o ${output}_sorted_read_name.bam

" >> $script
	
	#npairs=1
	#paired=TRUE
	#inputFiles[ o ]=2
	#inputFiles[ 1 ]=${output}_1.fastq
	#inputFiles[ 2 ]=${output}_2.fastq
	#inputFormat=STDFQ
	
	nfiles=1
	inputFiles[ 1 ]=${output}_sorted_read_name.bam
    fi


    ########################################## single end data
    if [[ "$paired" == "FALSE" || "$inputFormat" == "BAMPE" ]]; then
	
	if [[ "$nfiles" == "1" ]]; then

	    echo "
$novoalign -c $ncores -o SAM $'@RG\tID:${extraID}${code}\tSM:${extraID}${code}\tLB:${extraID}$code\tPL:ILLUMINA' $extra -F ${inputFormat} -f ${inputFiles[ 1 ]} -d $reference > ${output}.sam

$samtools view -bS -t ${fasta}i -o ${output}.bam ${output}.sam  ## make BAM file

$novosort -t $javaTemp -c $ncores -m ${memory}G  ${output}.bam -o ${output}_sorted.bam

" >> $script
	else   ##multiple single end files
	    
	    allBAMfiles=""
	    allBAMsofiles=""
	    allSAMfiles=""
	    
	    for i in `seq 1 $nfiles`; do
		ls -ltrh ${inputFiles[ i ]}

		echo "
$novoalign -c $ncores -o SAM $'@RG\tID:${extraID}${code}\tSM:${extraID}${code}\tLB:${extraID}$code\tPL:ILLUMINA' $extra -F ${inputFormat} -f ${inputFiles[ i ]} -d $reference > ${output}_${i}.sam

$samtools view -bS -t ${fasta}i -o ${output}_${i}.bam ${output}_${i}.sam  ## make BAM file

$samtools sort -m ${memory2}000000000 ${output}_${i}.bam ${output}_sorted_${i}  ## sort

" >> $script
    
		
		allBAMfiles="$allBAMfiles ${output}_${i}.bam"
		allBAMsofiles="$allBAMsofiles ${output}_sorted_${i}.bam"
		allSAMfiles="$allSAMfiles ${output}_${i}.sam"
	    done

	    echo "
$samtools merge -f ${output}_sorted.bam $allBAMsofiles 

rm $allBAMfiles $allSAMfiles $allBAMsofiles

" >> $script

	fi

	if [[ "$inputFormat" == "BAMPE" ]]; then
	    
echo "

rm ${output}_sorted_read_name.bam

" >> $script
	fi
    fi
    

    if [[ "$paired" == "TRUE" && ! "$inputFormat" == "BAMPE" ]]; then      #################################### now analysis of paired end data
			
	if [[ "$npairs" == "1" ]]; then   ###first option only one set of files
	  
	  echo "

$novoalign -c 11 -o SAM $'@RG\tID:${extraID}${code}\tSM:${extraID}${code}\tLB:${extraID}$code\tPL:ILLUMINA' $extra -F ${inputFormat} -f ${inputFiles[ 1 ]} ${inputFiles[ 2 ]}  -d $reference | ${samtools1} view -uSb - | $novosort - -t /scratch0/ -c $ncores -m 5G -i -o ${output}_sorted.bam 

" >> $script

      else 
	  allBAMfiles=""
	  allBAMsofiles=""
	  allSAMfiles=""
	  
	  for i in `seq 1 $npairs`; do
	      ((p1=2*i - 1))
	      ((p2=2*i - 0))
	      #echo $p1 $p2
	      echo "
$novoalign -c $ncores -o SAM $'@RG\tID:${extraID}${code}\tSM:${extraID}${code}\tLB:${extraID}$code\tPL:ILLUMINA' $extra -F ${inputFormat} -f ${inputFiles[ $p1 ]} ${inputFiles[ $p2 ]} -d $reference | $samtools view -uSb - | $novosort - -t /scratch0/ -c 1 -m 3G -i -o ${output}_sorted_${i}.bam

" >> $script
	      
	      allBAMsofiles="$allBAMsofiles ${output}_sorted_${i}.bam"
	  done
	  
	  echo "
$samtools merge -f ${output}_sorted.bam $allBAMsofiles 

rm $allBAMsofiles

${samtools} index  ${output}_sorted.bam   ##build index

" >> $script
      fi
    fi
    
    if [[ "$filtDup" == "TRUE" ]]; then
	echo "
## Now remove duplicates using PICARD
java -Xmx${memory}g -jar ${picardDup} ${picardJavaTemp} ASSUME_SORTED=true REMOVE_DUPLICATES=FALSE INPUT=${output}_sorted.bam OUTPUT=${output}_sorted_unique.bam METRICS_FILE=${output}_picard_metrics.out

${samtools1} index  ${output}_sorted_unique.bam   ##build index

rm ${output}_sorted.bam ${output}_sorted.bam.bai

" >> $script
	else
echo "
mv ${output}_sorted.bam ${output}_sorted_unique.bam
mv ${output}_sorted.bam.bai ${output}_sorted_unique.bam.bai

" >> $script	
    fi
    
    
      if [[ "$clean" == "TRUE" ]]; then
	  echo "
rm ${output}.bam ${output}.sam 
" >> $script	  
      fi

fi


if [[ "$fixUnique" == "yes" ]]; then ## a fix if I forgot to run the clonality stuff
    
    echo "
java -Xmx${memory}g -jar ${picardDup} ${picardJavaTemp} ASSUME_SORTED=true REMOVE_DUPLICATES=FALSE INPUT=${output}_sorted_unique.bam OUTPUT=${output}_sorted_unique2.bam METRICS_FILE=${output}_picard_metrics.out

mv ${output}_sorted_unique2.bam ${output}_sorted_unique.bam

${samtools} index  ${output}_sorted_unique.bam   ##build index

" >> $script

fi

if [[ "${summaryStats}" == "yes" ]]; then
    
    if [ -e ${baitFile} ]; then  ### summary stats if there is an interval list file 
	echo "
java -Xmx${memory}g -jar ${picardMetrics} BAIT_INTERVALS=${baitFile} TARGET_INTERVALS=${baitFile}  INPUT=${output}_sorted_unique.bam  OUTPUT=${output}.hybridMetrics
" >> $script
    else 
	echo "Baitfile ${baitFile} does not exist"	
    fi
    
fi
    

if [[ "$reduceReads" == "yes" ]]; then

    if [ ! -e $fasta ]; then echo "$fasta does not exist"; exit; fi
    echo "
${java17}  -Djava.io.tmpdir=${javaTemp} -Xmx${memory}g -jar $GATK -R $fasta -T ReduceReads -I ${output}_sorted_unique.bam -o ${output}_reduced.bam
" >> $script

fi

if [[ "$makegVCF" == "yes" ]]; then

echo "
$java17 -Djava.io.tmpdir=${javaTemp} -Xmx8g  -jar $GATK -T HaplotypeCaller -R $fasta -I ${output}_sorted_unique.bam \
     --dbsnp ${bundle}/dbsnp_137.b37.vcf \
     --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 \
     -stand_call_conf 30.0 \
     -stand_emit_conf 10.0 \
     -L $target \
     --activeRegionExtension 100 \
     -o ${output}.gvcf
" >> $script

fi


echo "date" >> $script  ##to measure the duration
echo $script

