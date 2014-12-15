

# Whole Genome Sequencing Pipeline

## Example

First create the support file which points to where the fastq files lie.
> bash makesupportfile.sh  data/ fq.gz > examples/support.txt

This should look like this:

```
code	f1	f2
sim_reads	data//sim_reads_1.fq.gz	data//sim_reads_2.fq.gz
```

Then do the alignment:

```bash
bash WGS_pipeline.sh 
     --mode align
     --supportFrame examples/support.txt
     --reference 1kg
     --tparam 320
     --inputFormat STDFQ 
     --extraID example
     --projectID example
```

This creates the following directory structure:
```
cluster
     submission
     out
     err
aligned
     sim_reads
```
     
This generates a script under `cluster/submissions/align.sh` containing an SGE job array which can then be submitted to the cluster:

> qsub cluster/submissions/align.sh

## Overview


The two main scripts are:
- WGS_submission_script.sh
- WGS_pipeline.sh

*WGS_submission_script* defines the command line arguments to call
*WGS_pipeline*
with.


> bash WGS_pipeline.sh 
     --mode [align|gvcf|jointvcf]
     --supportFrame < >
     --reference [hg38_noAlt|1kg]
     --tparam 320
     --inputFormat STDFQ 
     --extraID < >
     --projectID < >


The three modes of running the pipeline are either for alignment, for single variant calling or for joint variant calling:
- align
- gvcf
- jointvcf

If ran for alignment, *WGS_pipeline* will generate the following bash
scripts:

- cluster/submission/align.sh
- cluster/submission/align_table.sh

If ran for variant calling, the following bash scripts will be generated:

- cluster/submission/makeVCF.sh
- cluster/submission/makeVCF_table.sh

The first bash script is submitted to the cluster and uses the second script to specify the SGE job array.
Some parameters are common to both modes.
Following is an explanation of those parameters.

## Options

### supportFrame

The pipeline submission script assumes that you have written a "support" data frame.
There's a convenience script to generate it:

> bash makesupport.sh  ../data/ fq.gz > ../examples/support.txt

If you want to generate manually then the format of the first row is:

> code f1 f2 

where f1 and f2 are the paths to the paired end fastq files and code is an identifier for the exome of interest.
The code is used to create a subdirectory under aligned where the output goes.

### reference

The reference build on which to do the alignment.
Two reference builds are supported:
- hg38_noAlt
- 1kg

### extraID

A very useful parameter to set is "extraID" to know the batch.

## Alignment (align mode)

This generates the scripts:
 - cluster/submission/align.sh
 - cluster/submission/align_table.sh
 
*align.sh* is submitted to the queue using qsub and reads its job from the job array in *align_table.sh*.

The alignment is done using [novoalign](http://www.novocraft.com/main/page.php?s=novoalign) on the specified [reference build](reference) which generates a SAM.
The SAM is piped to [samblaster](https://github.com/GregoryFaust/samblaster) to mark duplicates and extract discordant and split reads.
The output is then piped to [samtools]() which generates the BAM file.
Finally [novosort]() is ran on the BAM.

The default output folder is:
> aligned

## Variant calling (gvcf mode)

This generates scripts:
- cluster/submission/makeVCF.sh
- cluster/submission/makeVCF_table.sh

The variant calling is done by [GATK](http://www.broadinstitute.org/partnerships/education/broade/best-practices-variant-calling-gatk) using HaplotypeCaller.

> **Note**: We don't use GenotypeGVCFs anymore because it requires all input BAM files which is too computationally expensive.  Instead joint calling is done using the gvcf files in the next step.

## Joint variant calling (jointvcf mode)

This generates scripts:


## Dependencies

- Alignment
  - novoalign
  - samblaster
  - [samtools]()
  - [Picard]()
- Variant calling
  - GATK

