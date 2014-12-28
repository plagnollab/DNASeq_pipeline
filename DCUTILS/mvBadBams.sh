# if a bam file is incomplete then running samtools view on it will produce an error saying there is no EOF and probably truncated
# this script expects bam files to be in subfolder of oFolder
# if this is not set then the default value below will be used
# any bam files which are incomplete will be moved to folder oFolder/trash

############ general folders, no need to update these
software=/cluster/project8/vyp/vincent/Software
samtools=${software}/samtools-1.1/samtools

bamFolder=/goon2/project99/bipolargenomes_raw/ingest

if [ .$oFolder == . ]
then
oFolder=$bamFolder/forlab/addRG
fi

echo Using $oFolder as location for bam files

trashFolder=$oFolder/trash
if [ ! -e $trashFolder ]; then mkdir $trashFolder; fi

find $oFolder -name '*.bam' -print | while read outFile ;
do 
    rm std.err
	echo nothing > std.err
	$samtools view -H $outFile > /dev/null 2> std.err
	truncated=`fgrep truncated std.err`
	# echo outFile=$outFile, truncated=$truncated.
	if [ ."$truncated" == . ]
	then 
	echo Will keep $outFile
	else 
	echo Will remove $outFile with
	echo mv $outFile $trashFolder
	# mv $outFile $trashFolder
	fi
done
