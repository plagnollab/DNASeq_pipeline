############ general folders, no need to update these
referenceF=/cluster/project8/vyp/vincent/data/reference_genomes
software=/cluster/project8/vyp/vincent/Software
pipeline=${software}/pipeline/WGS/WGS_pipeline.sh


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
reference=hg38_noAlt

############# below the actions that we want to apply
align=yes
makegVCF=no
makeVCF=no

force=no

supportFrame=support/Hardcastle_October2014.tab

########################## end of parameter file




njobs=0
echo "scriptNames" > $mainTable

sh ${pipeline} --supportFrame ${supportFrame} --reference ${reference} --align ${align}  --tparam 320  --inputFormat STDFQ  --extraID $extraID --makeVCF ${makeVCF} --makegVCF ${makegVCF}  --projectID ${projectID}
