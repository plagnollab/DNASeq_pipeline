referenceF=/scratch2/vyp-scratch2/reference_datasets/human_reference_sequence
#referenceF=/ugi/home/shared/vincent/reference_genome

software=/cluster/project8/vyp/vincent/Software
pipeline=${software}/pipeline/calling/align_pipeline_DNA_v6.sh


### align all the data
iFolder=fastq
oFolder=aligned
if [ ! -e $oFolder ]; then mkdir $oFolder; fi

find fastq/ -maxdepth 1 -mindepth 1 -exec basename {} \; > support/listIDs.tab
myIDs=`cat support/listIDs.tab`


projectID=Hardcastle
mainScript=cluster/submission/${projectID}_main.sh
mainTable=cluster/submission/${projectID}_table.tab

njobs=0

echo "scriptNames" > $mainTable


for nID in $myIDs; do
    echo $nID
    ((njobs=njobs+1))
    folder=${iFolder}/${nID}/clean_data
    
    echo "Folder is $folder"
    nfiles=0
    inputFiles=""
    for file1 in `find $folder -name *_1.clean.fq.gz | sort`; do
	((nfiles=nfiles+2))
	file2=`echo $file1 | sed -e 's/_1.clean.fq.gz/_2.clean.fq.gz/g'`
	inputFiles="$inputFiles $file1 $file2"	
    done
    #ls -ltrh $inputFiles; exit
    inputFiles="$nfiles $inputFiles"

    fasta=${referenceF}/human_g1k_v37.fasta
    reference=${referenceF}/human_g1k_v37.fasta.k15.s2.novoindex
    output=${oFolder}/${nID}/${nID}
    
    if [ ! -e ${oFolder}/${nID} ]; then mkdir ${oFolder}/${nID}; fi
        
    align=yes
    summaryStats=no
    makegVCF=no
    extraID=Shamima_April2014

    if [ ! -s ${oFolder}/${nID}/${nID}_sorted_unique.bam ]; then 
	echo $nID
	##ILM1.8
	
	sh ${pipeline} --inputFiles ${inputFiles}  --fasta ${fasta} --reference ${reference} --output ${output} --align ${align}  --summaryStats ${summaryStats} --tparam 320 --inputFormat STDFQ  --extraID $extraID --makegVCF ${makegVCF}  --projectID ${projectID}
	
    fi


done



################# final steps


echo "#$ -t 1-${njobs}

array=( \`cat \"${mainTable}\" \`)

script=\${array[ \$SGE_TASK_ID ]}

echo \$script

sh \$script

" >> $mainScript


echo "Final script in $mainScript"

#qsub $mainScript
