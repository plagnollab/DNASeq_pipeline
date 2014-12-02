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
align=yes
makegVCF=no
makeVCF=no

force=no

supportFrame=support/Hardcastle_October2014.tab

########################## end of parameter file




njobs=0
echo "scriptNames" > $mainTable

sh ${pipeline} --supportFrame ${supportFrame} --fasta ${fasta} --reference ${reference} --align ${align}  --tparam 320  --inputFormat STDFQ  --extraID $extraID --makeVCF ${makeVCF} --makegVCF ${makegVCF}  --projectID ${projectID}
