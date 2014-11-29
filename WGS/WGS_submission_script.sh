############ general folders, no need to update these
referenceF=/cluster/project8/vyp/vincent/data/reference_genomes
software=/cluster/project8/vyp/vincent/Software
pipeline=${software}/pipeline/WGS/WGS_pipeline.sh
fasta=/scratch2/vyp-scratch2/reference_datasets/human_reference_sequence/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
reference=/scratch2/vyp-scratch2/reference_datasets/human_reference_sequence/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.k15.s2.novoindex



################## below location of the input and output folders, as well as the list of IDs
iFolder=/SAN/biomed/biomed2/biomed2/ingest
oFolder=/scratch2/vyp-scratch2/WGS/Hardcastle_October2014/aligned
if [ ! -e $oFolder ]; then mkdir $oFolder; fi
find /SAN/biomed/biomed2/biomed2/ingest/ -name *1.fastq.gz -exec basename {} .fastq.gz \; | sed -e 's/_L0.*\|_1$//g' > support/listIDs.tab
myIDs=`cat support/listIDs.tab | grep 14`



#################### key parameters to choose
projectID=HardcastleStep2
mainScript=cluster/submission/${projectID}_main.sh
mainTable=cluster/submission/${projectID}_table.tab
extraID=Hardcastle_October2014


############# below the actions that we want to apply
align=no
makegVCF=yes


########################## end of parameter file




njobs=0
echo "scriptNames" > $mainTable

for nID in $myIDs; do
    echo $nID
    ((njobs=njobs+1))
    folder=${iFolder}/
    
    echo "Folder is $folder"
    nfiles=0
    inputFiles=""
    for file1 in `find $folder -name ${nID}*_R1.fastq.gz | sort`; do
	((nfiles=nfiles+2))
	file2=`echo $file1 | sed -e 's/_R1.fastq.gz/_R2.fastq.gz/g'`
	inputFiles="$inputFiles $file1 $file2"	
    done
    #ls -ltrh $inputFiles; exit
    inputFiles="$nfiles $inputFiles"

    output=${oFolder}/${nID}/${nID}
    
    if [ ! -e ${oFolder}/${nID} ]; then mkdir ${oFolder}/${nID}; fi
        
    echo "Looking if ${oFolder}/${nID}/${nID}_sorted_unique.bam exists"
    if [ -s ${oFolder}/${nID}/${nID}_sorted_unique.bam.bai ]; then 
	echo $nID
	echo "${oFolder}/${nID}/${nID}_sorted_unique.bam already exists"
	##ILM1.8
	
	sh ${pipeline} --inputFiles ${inputFiles}  --fasta ${fasta} --reference ${reference} --output ${output} --align ${align}  --tparam 320  --inputFormat STDFQ  --extraID $extraID --makegVCF ${makegVCF}  --projectID ${projectID}
    fi


done

if [[ "$makegVCF" == "yes" ]]; then
    ((njobs=njobs*24))
fi

njobs=`wc -l cluster/submission/${projectID}_table.tab | cut -f 1 -d' '`
((njobs=njobs-1)) ##because there is a header line to remove

################# final steps


echo "#$ -t 1-${njobs}

array=( \`cat \"${mainTable}\" \`)

script=\${array[ \$SGE_TASK_ID ]}

echo \$script

sh \$script

" >> $mainScript


echo "Final script in $mainScript"

#qsub $mainScript
