```
BASEDIR=/goon2/scratch2/vyp-scratch2/RoseRichardson/1327-UCL-Fam2_Report/
```

Make the support file which points to the fastq files.

```
bash /scratch2/vyp-scratch2/Software/DNASeq_pipeline/WGS/makesupportfile.sh $BASEDIR/Fastq fastq.gz > support.txt
```

# Align


```
reference=1kg
projectID=RoseRichardson_Fam2
extraID=RoseRichardson_Fam2
supportFrame=$BASEDIR/support.txt
```

```
bash /goon2/scratch2/vyp-scratch2/Software/DNASeq_pipeline/WGS/WGS_pipeline.sh --mode align --supportFrame $supportFrame --reference $reference --aligner-tparam 320 --inputFormat STDFQ --projectID $projectID --extraID $extraID
```


