# *************************************************************
#
#
# *************************************************************
# Directory where your paired end fq.gz lie.
DIR=$1
paste <(echo code) <(echo f1) <(echo f2)
paste <(ls -1 $DIR/*_1.fq.gz | sed 's/_1.fq.gz//' | xargs -I f basename f) <(ls -1 $DIR/*_1.fq.gz) <(ls -1 $DIR/*_2.fq.gz)
