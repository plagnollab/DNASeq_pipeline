# Takes as a command line argument the directory where the paired-end fastq files are.
# Just redirect the output to the support file passed to the WGS pipeline.
# How to run:
# bash makesupportfile.sh <dir to fastq files>
 
# Directory where your paired end fq.gz lie.
DIR=$1
if [[ "$2" == "" ]]
then
    ext=fastq.gz
else
    ext=$2
fi
paste <(echo code) <(echo f1) <(echo f2)
paste <(ls -1 $DIR/*_1.${ext} | sed "s/_1.${ext}//" | xargs -I f basename f) <(ls -1 $DIR/*_1.${ext}) <(ls -1 $DIR/*_2.${ext})
