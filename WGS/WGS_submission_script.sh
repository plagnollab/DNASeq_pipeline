############ general folders, no need to update these
referenceF=/cluster/project8/vyp/vincent/data/reference_genomes
software=/cluster/project8/vyp/vincent/Software
pipeline=${software}/pipeline/WGS/WGS_pipeline.sh


################## below location of the input and output folders, as well as the list of IDs
#iFolder=/SAN/biomed/biomed2/biomed2/ingest
projectID=/scratch2/vyp-scratch2/exomes_temp/IoO/Manchester_batch1
if [ ! -e $oFolder ]; then mkdir $oFolder; fi


#################### key parameters to choose
extraID=Manchester_December2014
reference=1kg
supportFrame=support/Manchester_batch1.tab

########################## end of parameter file

njobs=0
echo "scriptNames" > $mainTable

mkdir -p aligned/${projectID}

bash ${pipeline} --mode align --supportFrame ${supportFrame} --reference ${reference}  --tparam 320  --inputFormat STDFQ  --extraID $extraID --projectID ${projectID}

#bash ${pipeline} --mode vcf --supportFrame ${supportFrame} --reference ${reference} --extraID $extraID  --projectID ${projectID}


