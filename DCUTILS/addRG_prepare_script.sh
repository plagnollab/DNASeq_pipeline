# DC script to use picard to add read groups to bam files which did not have them
# checks to see target bam file does not already exist before overwriting it
# run this then submit addRG/submission/addRG.sh
# when submitted jobs have finished, use mvBadBams.sh to move any incomplete target bams to trash directory then rm -f addRG (to get rid of *out, *err files) then run this again and submit again


############ general folders, no need to update these
software=/cluster/project8/vyp/vincent/Software
AddOrReplaceReadGroups=${software}/picard-tools-1.100/AddOrReplaceReadGroups.jar
java17=/share/apps/jdk1.7.0_45/jre/bin/java

homeFolder=/cluster/project8/bipolargenomes
tempFolder=/scratch2/vyp-scratch2/vincent/temp/novoalign

nhours=10
ncores=1
vmem=6 ##DC changed from 1 to try to get makegVCF to work
memory=2
memory2=5  ##used for the sort function, seem to crash when using 10
queue=queue6
scratch=0
 
bamFolder=/goon2/project99/bipolargenomes_raw/ingest
iFolder=$bamFolder
oFolder=$homeFolder
oFolder=$bamFolder/forlab/addRG
# I am writig these back on goon2 to avoid restrictions on my quota

if [ ! -e $oFolder ]; then mkdir $oFolder; fi
supportFrame=$oFolder/addRGsupport.tab

## DC bit because we already have bam files:
if [ -e $supportFrame ] ; then rm $supportFrame ; fi
echo code f1 f2 > $supportFrame
ls $bamFolder/*/Assembly/genome/bam/*.bam | while read fullName ;
	do
	echo $fullName; 
	ID=${fullName##*/};
	ID=${ID%.*};
	outFile=$oFolder/aligned/$ID/${ID}_sorted_unique.bam; 
	if [ ! -e $oFolder/aligned ]; then mkdir $oFolder/aligned; fi;
	if [ ! -e $oFolder/aligned/$ID ]; then mkdir $oFolder/aligned/$ID; fi; \
	echo $ID $fullName $outFile >> $supportFrame
	done

cd $homeFolder
if [ ! -e $homeFolder/addRG ] ; then mkdir $homeFolder/addRG ; fi
if [ ! -e $homeFolder/addRG/submission ] ; then mkdir $homeFolder/addRG/submission ; fi
if [ ! -e $homeFolder/addRG/out ] ; then mkdir $homeFolder/addRG/out ; fi
if [ ! -e $homeFolder/addRG/error ] ; then mkdir $homeFolder/addRG/error ; fi

    mainScript=addRG/submission/addRG.sh
    mainTable=addRG/submission/addRG_table.sh
if [ -e $mainTable ] ; then rm  $mainTable ; fi
echo firstLineToBeIgnored > $mainTable
# this is because array is indexed from 0 but tasks are indexed from 1
    tail -n+2 $supportFrame | while read code inFile outFile
    do
	output=${oFolder}/${code}/${code}
	##one job per chromosome to save time
           script=`echo $mainScript | sed -e 's/.sh$//'`_${code}.sh
           echo $script >> $mainTable
           echo " 
		   if [ -e $outFile ] 
		   then 
		   echo $outFile already exists, exiting...
		   exit 
		   fi 
           $java17 -Djava.io.tmpdir=${tempFolder} -Xmx4g -jar $AddOrReplaceReadGroups \
		   I=$inFile \
		   O=$outFile \
		   LB=$code PL=illumina PU=run SM=$code \
		   VALIDATION_STRINGENCY=LENIENT \
           " > $script
    done
    #end of while loop

njobs=`cat $mainTable | wc -l`
echo "
#!/bin/bash
#$ -S /bin/bash
#$ -o addRG/out
#$ -e addRG/error
#$ -cwd
#$ -pe smp ${ncores}
#$ -l scr=${scratch}G
#$ -l tmem=${vmem}G,h_vmem=${vmem}G
#$ -l h_rt=${nhours}:0:0
#$ -tc 25
#$ -t 1-${njobs}
#$ -V
#$ -R y
array=( \`cat \"${mainTable}\" \`)
script=\${array[ \$SGE_TASK_ID ]}
root=\${script##*/};
root=\${root%.*};
date
echo \$script
sh \$script  1> addRG/out/\$root.out 2> addRG/error/\$root.err
date
" > $mainScript


    echo "Main submission scripts and tables"
    wc -l $mainScript $mainTable

