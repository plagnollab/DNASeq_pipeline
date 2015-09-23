# Parameters (to be modified)

This is an example of a RUNME.md that you can copy to the base directory from where you will be running the pipeline.
Only the following parameters need to be modified.  The rest should be standard.

```
SOFTWARE=/cluster/project8/vyp/Software/
BASEDIR=/goon2/scratch2/vyp-scratch2/RoseRichardson/1327-UCL-Fam2_Report/
cp $SOFTWARE/DNASeq_pipeline/WGS/RUNME.md .
```

```
reference=1kg
projectID=RoseRichardson_Fam2
extraID=RoseRichardson_Fam2
```

# Align
Make the support file which points to the fastq files.
```
bash $SOFTWARE/DNASeq_pipeline/WGS/makesupportfile.sh $BASEDIR/Fastq fastq.gz > ${BASEDIR}/support.txt
```
Generate the alignment script:
```
bash $SOFTWARE/DNASeq_pipeline/WGS/WGS_pipeline.sh --mode align --supportFrame ${BASEDIR}/support.txt --reference ${reference} --aligner-tparam 320 --inputFormat STDFQ --projectID ${projectID} --extraID ${extraID}
```
Submit the script:
```
qsub ${projectID}/align/scripts/align.sh
```

# Create the gVCFs

```
bash $SOFTWARE/DNASeq_pipeline/WGS/WGS_pipeline.sh --mode gvcf --supportFrame ${BASEDIR}/support.txt --reference ${reference}t --inputFormat STDFQ --projectID ${projectID} --extraID ${extraID}
```
The default memory is 7.8G but if you need more memory you can specify it on the command line:
```
bash $SOFTWARE/DNASeq_pipeline/WGS/WGS_pipeline.sh --mode gvcf --supportFrame ${BASEDIR}/support.txt --reference ${reference}t --inputFormat STDFQ --projectID ${projectID} --extraID ${extraID} --mem 10
```

# CombineGVCFs for UCLEX
```
bash $SOFTWARE/DNASeq_pipeline/WGS/WGS_pipeline.sh --mode CombineGVCFs --supportFrame ${BASEDIR}/support.txt --reference ${reference} --inputFormat STDFQ --projectID ${projectID} --extraID ${extraID}
```
