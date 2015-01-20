#!/bin/bash

#$ -S /bin/bash
#$ -o /cluster/project8/vyp/AdamLevine/coverage/log
#$ -e /cluster/project8/vyp/AdamLevine/coverage/log
#$ -l h_rt=4:0:0
#$ -l tmem=4G,h_vmem=4G
#$ -t 1-85

i=$SGE_TASK_ID

cd /cluster/project8/vyp/AdamLevine/coverage/results

array=( $( cat ../Levine.txt ) )
j=$((i-1))
this=${array[$j]}

#needs coverageBed from bedtools
#coverage=~/cluster/software/bedtools2-2.20.1/bin/coverageBed

bam=/scratch2/vyp-scratch2/exomes_temp/Levine/aligned/${this}/${this}_sorted_unique.bam

coverageBed -abam ${bam} -b ../intersect.bed -d > ${this}_out.txt

# wc ${this}_out.txt | awk '{print $1}' > ${this}_counts.txt

for i in 0 1 10 20;
do
    awk -v x=${i} '$5>=x' ${this}_out.txt | wc | awk '{print $1}'  >> ${this}_counts.txt
done

awk '{ total += $5; count++ } END { print total/count }' ${this}_out.txt > ${this}_mean.txt

rm ${this}_out.txt

