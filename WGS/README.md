

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
     --projectID Project1
```

This creates the following directory structure:
```
Project1/
    README
    support.txt
     align/
        scripts/
        data/
        out/
        err/
```
`support.txt` is just a copy of the support file that was specified.
`README` should get automatically updated but this has not been implemented yet (also not sure what information would be useful).
If all the fastq.gz files specified in support.txt are found then this generates the script `Project1/align/scripts/align.sh` containing an SGE job array, where each job does alignment for each sample.
Check the output in `Project1/align/scripts`.
This script can then be submitted to the cluster:

> qsub Project1/align/scripts/align.sh

If certain BAM files already exist, they will not be regenerated unless the `--force` option is used.  This `--force` argument applies to the modes too: gvcf, combinegvcf etc.
Once the jobs have been submitted use `qstat` to check the running status.
Each job in the array will write logs, to both `Project1/align/out/` and `Project1/align/err/` respectively.

Once the job array finishes you can check that all files exists by running the bash script again:

```bash
bash WGS_pipeline.sh 
     --mode align
     --supportFrame examples/support.txt
     --reference 1kg
     --tparam 320
     --inputFormat STDFQ 
     --projectID Project1
     --extraID Project1_
```
The job array should now be empty if all files have been successfully created.
If they are still scripts to be submitted, first check the log files the type of error.
If its a transient error like a stale NFS handle due to a lost connection then submit the jobs again.

Once the alignment is complete you can move on to the gvcf generation:
```bash
bash WGS_pipeline.sh 
     --mode gvcfs
     --supportFrame examples/support.txt
     --reference 1kg
     --tparam 320
     --inputFormat STDFQ 
     --projectID Project1
```
This will follow a similar process to what has been described here with the same dir structure created.


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

The parameter "extraID" defines the batch.
It is added to the RG header of the BAM files and is used to distinguish samples in the merged VCF files.


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

## Single sample variant calling

The variant calling is done by [GATK](http://www.broadinstitute.org/partnerships/education/broade/best-practices-variant-calling-gatk) using [HaplotypeCaller](https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php).

The output of HaplotypeCaller will need to be filtered by variant recalibration (best) or hard-filtering before use in downstream analyses.

The VCF file contains the following information:
```
##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities">
##INFO=<ID=ClippingRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref number of hard clipped bases">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
##INFO=<ID=DS,Number=0,Type=Flag,Description="Were any of the samples downsampled?">
##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">
##INFO=<ID=HaplotypeScore,Number=1,Type=Float,Description="Consistency of the site with at most two segregating haplotypes">
##INFO=<ID=InbreedingCoeff,Number=1,Type=Float,Description="Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared again
st the Hardy-Weinberg expectation">
##INFO=<ID=MLEAC,Number=A,Type=Integer,Description="Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for eac
h ALT allele, in the same order as listed">
##INFO=<ID=MLEAF,Number=A,Type=Float,Description="Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for ea
ch ALT allele, in the same order as listed">
##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">
##INFO=<ID=MQ0,Number=1,Type=Integer,Description="Total Mapping Quality Zero Reads">
##INFO=<ID=MQRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities">
##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias">
```

> **Note**: Joint calling is done using the full set of gvcf files in the next step, not usually done within a batch.

## Merging of samples

This is done by CombineGVCFs

## Joint variant calling (jointvcf mode)

Single samples VCFs called by HaplotypeCaller are combined per chromosome using GenotypeGVCFs.
After this step, the generated joint VCFs cannot be merged using CombineGVCF.


## Variant selection, filtering and recalibration

All these steps are handled by the msample script.

## Dependencies

- Alignment
  - novoalign
  - samblaster
  - [samtools]()
  - [Picard]()
- Variant calling
  - GATK

