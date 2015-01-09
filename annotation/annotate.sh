#! /bin/bash
#
# VCF filename is the argument

#get dir in which script is in
set -u
set -e
set -x

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

in=$1

#expect the chromsome to be part of the filename
chr=`echo $in | sed 's/.*_chr\(.*\).gvcf.gz/\1/'`

out=${in%.gvcf.gz}-single.vcf
python $DIR/multiallele_to_single_gvcf.py $in > $out

in=$out
out=${out%.vcf}-VEP.gvcf
bash $DIR/run_VEP.sh $in $chr > $out

#groups.txt file needs to be generated
python $DIR/extract_VEP.py groups.txt $in > $out
