# Multiple Sample Calling

Once the gvcf files have been generated using WGS, the samples can then be joined and eventually integrated in the
next UCL Exome release.

## Step 0 - Merging gVCFs into combined gVCF files

First the gVCFs obtained from different batches can be merged with the following script:
```
merge_gVCF.sh
```
This script takes as input a list of gVCF files and joins them per chromosome.
The joining is done by [the GATK CombineGVCFs tool](https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_variantutils_CombineGVCFs.php).
The output goes to the combinedFolder per chromosome:

```
/scratch2/vyp-scratch2/vincent/GATK/HC/combinedVCFs/
```

The next step is to do the joint genotyping and the variant recalibration.
The output goes to:
```
/scratch2/vyp-scratch2/vincent/GATK/mainset_${currentUCLex}/mainset_${currentUCLex}
```

## Step 1: GenotypegVCF 

## Step 2: Recalibration 

## Step 3: Annotation and custom filtering

## Step 4: Convert to R snpStats

