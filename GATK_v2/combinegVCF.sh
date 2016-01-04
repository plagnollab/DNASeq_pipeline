fasta=/scratch2/vyp-scratch2/reference_datasets/human_reference_sequence/human_g1k_v37.fasta
bundle=/scratch2/vyp-scratch2/reference_datasets/GATK_bundle

#nExPergVCF=105 ##for GOS
nExPergVCF=122 ##Shamima
#nExPergVCF=100 ##Levine_Aug2014
#nExPergVCF=100 ##IoO

java=/share/apps/jdk1.7.0_45/bin/java
target=/cluster/project8/vyp/exome_sequencing_multisamples/target_region/data/merged_exome_target_cleaned.bed
computerChoice=none

combinedFolder=/scratch2/vyp-scratch2/vincent/GATK/HC/combinedVCFs ##default value
GATK=/cluster/project8/vyp/vincent/Software/GenomeAnalysisTK-3.5-0/GenomeAnalysisTK.jar


output=/scratch2/vyp-scratch2/vincent/GATK/mainset_${currentUCLex}/mainset_${currentUCLex}
#output=/scratch2/vyp-scratch2/vincent/GATK/cardioset_${currentUCLex}/cardioset_${currentUCLex}

until [ -z "$1" ]; do
    # use a case statement to test vars. we always test $1 and shift at the end of the for block.
    case $1 in
	--combinedFolder )  ## output folder that contains the merged gVCF files
	    shift
	    combinedFolder=$1;;
	--force )
	    shift
	    force=$1;;
	--folderCode )
	    shift
	    folderCode=$1;;
	--iFolder ) ## folder that contains the gVCF files
            shift
            iFolder=$1;;
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

script=cluster/submission/comb_${folderCode}.sh

VCFlist=${iFolder}/VCFlist.list
find ${iFolder}/ -name \*gvcf.gz -o -name \*g.vcf.gz | sed -e 's/.tbi$//g' > $VCFlist

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
		##VCFcode=`basename $VCF .gvcf.gz`
		##ogVCF=${iFolder}/${VCFcode}.gvcf.gz
		echo "   --variant ${VCF} \\" >> $script
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

