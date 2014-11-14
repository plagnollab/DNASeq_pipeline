iFolder=/scratch2/vyp-scratch2/vincent/GATK/HC/IoO
fasta=/scratch2/vyp-scratch2/reference_datasets/human_reference_sequence/human_g1k_v37.fasta
bundle=/scratch2/vyp-scratch2/reference_datasets/GATK_bundle

#nExPergVCF=105 ##for GOS
nExPergVCF=122 ##Shamima
#nExPergVCF=100 ##Levine_Aug2014
#nExPergVCF=100 ##IoO

java=/share/apps/jdk1.7.0_45/bin/java
tmpDir=/scratch0/vyp
target=/cluster/project8/vyp/exome_sequencing_multisamples/target_region/data/merged_exome_target_cleaned.bed
computerChoice=none

combinedFolder=/scratch2/vyp-scratch2/vincent/GATK/HC/combinedVCFs
GATK=/cluster/project8/vyp/vincent/Software/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar
#GATK=/cluster/project8/vyp/vincent/Software/GenomeAnalysisTK-nightly-2014-08-24-g81e689c/GenomeAnalysisTK.jar

submit=no

runHC=no
mergeVCF=no
checkVCF=no
force=no
genotype=yes
gVCFlist=none

output=/scratch2/vyp-scratch2/vincent/GATK/mainset_${currentUCLex}/mainset_${currentUCLex}
#output=/scratch2/vyp-scratch2/vincent/GATK/cardioset_${currentUCLex}/cardioset_${currentUCLex}

until [ -z "$1" ]; do
    # use a case statement to test vars. we always test $1 and shift at the end of the for block.
    case $1 in
	--computerChoice )
	    shift
	    computerChoice=$1;;
	--combinedFolder )
	    shift
	    combinedFolder=$1;;
	--folderCode )
	    shift
	    folderCode=$1;;
	--submit )
	    shift
	    submit=$1;;
	--runHC )
	    shift
	    runHC=$1;;
	--checkVCF )
	    shift
	    checkVCF=$1;;
	--force )
	    shift
	    force=$1;;
	--mergeVCF )
	    shift
	    mergeVCF=$1;;
	--iFolder )
            shift
            iFolder=$1;;
	--target )
            shift
            target=$1;;
        --BAMlist )
            shift
            BAMlist=$1;;
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
BAMlist=${iFolder}/BAMlist.list
gVCFFolder=${iFolder}/gVCF



if [ ! -e $BAMlist ]; then
    echo "Cannot find the BAMlist $BAMlist"
    exit
fi

for folder in $gVCFFolder; do
    if [ ! -e $folder ]; then mkdir $folder; fi
done



if [[ "$checkVCF" == "yes" ]]; then

    quick="TRUE"

    checkFile=${iFolder}/check.log
    if [ -e $checkFile ]; then rm $checkFile; fi

    for BAM in `cat $BAMlist`; do
	echo $BAM
	BAMcode=`basename $BAM _sorted_unique.bam`
	ogVCF=${gVCFFolder}/${BAMcode}.gvcf
	
	if [[ "$quick" == "TRUE" ]]; then
	    
	    if [ ! -e ${ogVCF}.gz.tbi ]; then
		echo "Missing index for ${ogVCF}.gz" >> ${checkFile} 
	    fi
	else 
		
	    if [ ! -e ${ogVCF}.gz ]; then
		chr="NA NA"
	    else
		chr=`zcat ${ogVCF}.gz | tail -1 | cut -f1,2`
	    fi
	    echo ${ogVCF}.gz ${chr} >> ${checkFile}
	fi
    done
    wc -l $checkFile

    echo "Potential issues flagged below"
    
    if [[ "$quick" == "TRUE" ]]; then
	cat ${checkFile}
    else 
	awk '{if ($2 != "Y") print}'  ${checkFile}
    fi

fi



if [[ "$runHC" == "yes" ]]; then

    mainScript=cluster/submission/HC_${folderCode}.sh
    mainTable=cluster/submission/HC_${folderCode}.tab
    
    echo "scripts" > $mainTable
    

    
    echo "
#$ -o cluster/out
#$ -e cluster/error
#$ -S /bin/bash
#$ -l h_vmem=11G,tmem=11G
#$ -l h_rt=96:0:0
#$ -R y
#$ -pe smp 1
#$ -cwd
#$ -l scr=10G
#$ -tc 60 " > $mainScript

    if [[ "$computerChoice" != "none" ]]; then echo "#$ -l hostname=${computerChoice}-*.local" >> $mainScript; fi
    
    njobs=0
    for BAM in `cat $BAMlist`; do
	echo $BAM
	BAMcode=`basename $BAM _sorted_unique.bam`
	ogVCF=${gVCFFolder}/${BAMcode}.gvcf
	
	locDir=/scratch0/${BAMcode}

	if [[ "$force" == "yes" || ! -s ${ogVCF}.gz.tbi ]]; then
	    
	    script=cluster/submission/${folderCode}_${BAMcode}.sh
	    ((njobs=njobs+1))
	    echo $script >> $mainTable
	    
	    echo "

if [ ! -e $locDir ]; then mkdir $locDir; fi

$java -Djava.io.tmpdir=$locDir -Xmx7g  -jar $GATK -T HaplotypeCaller -R $fasta -I $BAM --doNotRunPhysicalPhasing \
     --dbsnp ${bundle}/dbsnp_137.b37.vcf \
     --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 \
     -stand_call_conf 30.0 \
     -stand_emit_conf 10.0 \
     --activeRegionExtension 100 \
     --GVCFGQBands 10 --GVCFGQBands 20 --GVCFGQBands 60 \
     -o ${ogVCF}.gz

rm -rf $locDir

" > $script
	    
	    ls -ltrh $script
	fi
    done
    

    
    echo "#$ -t 1-${njobs}

array=( \`cat \"${mainTable}\" \`)

script=\${array[ \$SGE_TASK_ID ]}

echo \$script

sh \$script

" >> $mainScript
    

    ls -ltrh $mainTable $mainScript
    wc -l $mainTable $mainScript
    if [[ "$submit" == "yes" ]]; then qsub $mainScript; fi
fi




if [[ "$mergeVCF" == "yes" ]]; then

    script=cluster/submission/comb_${folderCode}.sh
    
    VCFlist=${iFolder}/VCFlist.list
    find ${iFolder}/gVCF/ -name \*vcf.gz.tbi | sed -e 's/.tbi$//g' > $VCFlist
    
    nVCFs=$(cat $VCFlist | wc -l)
    echo "Number of VCF files: $nVCFs"
    jobID=0
    ((startline=jobID*nExPergVCF))
    ((maxline=nExPergVCF*(1+jobID) ))

    mainScript=cluster/submission/combm_${folderCode}.sh
    mainTab=cluster/submission/combm_${folderCode}.tab

    echo "scripts" > $mainTab
    
    totjobs=0
    while [[ "$startline" -lt "$nVCFs" ]]; do
	((jobID=jobID+1))
	echo "Lines $startline $maxline"
	
	for chr in `seq 1 22` X Y MT; do
	    
	    oFile=${combinedFolder}/chr${chr}/${folderCode}_${jobID}.gvcf.gz
	    ls -ltrh ${oFile}.tbi
	    if [[ ! -e ${oFile}.tbi || "$force" == "yes" ]]; then
		
		((totjobs=totjobs+1))
		
		script=cluster/submission/com_${jobID}_chr${chr}_${folderCode}.sh
		
		echo "
$java -Xmx7g -jar $GATK \\
   -R $fasta \\
   -L $chr \\
   -T CombineGVCFs \\" > $script
	
		awk '{if ( (NR > '$startline') && (NR <= '$maxline')) print}' $VCFlist | while read VCF; do
		    VCFcode=`basename $VCF .gvcf.gz`
		    ogVCF=${gVCFFolder}/${VCFcode}.gvcf.gz
		    echo "   --variant ${ogVCF} \\" >> $script
		done
		
		echo "   -o ${combinedFolder}/chr${chr}/${folderCode}_${jobID}.gvcf.gz " >> $script
		echo $script >> $mainTab	    
		echo "Script in $script"
	    fi
	done
	((maxline=nExPergVCF + maxline))
	((startline=nExPergVCF + startline))
    done


    echo "Number of jobs $jobID"


    echo "
#$ -o cluster/out
#$ -e cluster/error
#$ -S /bin/bash
#$ -l h_vmem=9G,tmem=9G
#$ -l h_rt=30:0:0
#$ -R y
#$ -pe smp 1
#$ -cwd
#$ -tc 80
#$ -t 1-${totjobs} " > $mainScript

    if [[ "$computerChoice" != "none" ]]; then echo "#$ -l hostname=${computerChoice}-*.local" >> $mainScript; fi

echo "

array=( \`cat \"${mainTab}\" \`)

script=\${array[ \$SGE_TASK_ID ]}

echo \$script

sh \$script

" >> $mainScript



    wc -l $mainScript $mainTab
fi


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
	    
	    if [[ "$computerChoice" != "none" ]]; then echo "#$ -l hostname=${computerChoice}-*.local" >> $script; fi
	    
	    echo "

$java -Xmx9g -jar $GATK \\
   -R $fasta \\
   -T GenotypeGVCFs \\
   -L $chr -L $target --interval_set_rule INTERSECTION --interval_padding 100  \\
   --annotation InbreedingCoeff --annotation QualByDepth --annotation HaplotypeScore --annotation MappingQualityRankSumTest --annotation ReadPosRankSumTest --annotation FisherStrand \\
   --dbsnp ${bundle}/dbsnp_137.b37.vcf \\" >> $script
	    
	
	    cat support/gVCF.tab | while read path id; do
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