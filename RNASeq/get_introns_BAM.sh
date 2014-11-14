
slavePerl=/cluster/project4/vyp/vincent/Software/pipeline/RNASeq/get_introns_BAM_slave.pl
BAM=default
samtools=/cluster/project4/vyp/vincent/Software/samtools-0.1.18/samtools
tempFile=temp
outputFile=junctions.default
region=1:1-1000
script=default.sh
code=default
reducedOutput=TRUE

until [ -z "$1" ]; do
    # use a case statement to test vars. we always test $1 and shift at the end of the for block.
    case $1 in
	--BAM )
	    shift
	    BAM=$1;;
	--tempFile )
	    shift
	    tempFile=$1;;
	--outputFile )
	    shift
	    outputFile=$1;;
	--region )
	    shift
	    region=$1;;
	--script )
	    shift
	    script=$1;;
	--code )
	    code=$1;;
	--reducedOutput )
	    shift;;
	-* )
	    echo "Unrecognized option: $1"
	        exit 1;;
    esac
    shift
    if [ "$#" = "0" ]; then break; fi
done

if [ ! -e $BAM ]; then
    echo "$BAM does not exist"
fi

if [[ "$code" == "default" ]]; then
    code=`basename $BAM .bam`
fi 




if [[ "$reducedOutput" == "TRUE" ]]; then
    echo "
${samtools} view -h $BAM $region | awk '{if ( (\$1 !~ /^@/) &&  (\$6 ~ /N/) ) {print;}}' | $slavePerl | awk '{print \$4\"_\"\$5}'| sort | uniq -c > $outputFile
"  > $script
    else 
    echo "
${samtools} view -h $BAM $region | awk '{if ( (\$1 !~ /^@/) &&  (\$6 ~ /N/) ) {print;}}' | $slavePerl > $outputFile
" > $script
fi



echo "
echo \"Output in $outputFile\"

" >> $script

echo $script