# checks on BAM files to see if the output is valid.

BAM=$1

Software=/cluster/project8/vyp/vincent/Software
alias samtools=${Software}/samtools-1.1/samtools

samtools view -H $1



