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

This step does the joint calling on multiple gVCFs:

```
$java -Xmx${memoSmall}g -jar $GATK -R $fasta
-T GenotypeGVCFs -L $chr -L $target --interval_set_rule INTERSECTION --interval_padding 100 --annotation InbreedingCoeff --annotation QualByDepth --annotation HaplotypeScore --annotation MappingQualityRankSumTest
--annotation ReadPosRankSumTest --annotation FisherStrand
 --dbsnp ${bundle}/dbsnp_137.b37.vcf
 --variant $gVCF
 -o ${output}_chr${chr}.vcf.gz
```

## Step 2: Recalibration 

Extract the indels:
```
$java  -Djava.io.tmpdir=${tmpDir} -Xmx${memoSmall}g -jar ${GATK} -R $fasta
 -T SelectVariants 
 -V ${output}_chr${chr}.vcf.gz 
 -selectType INDEL 
 -selectType MIXED 
 -o ${output}_chr${chr}_indels.vcf.gz
```

Apply the filters for the indels:
```
$java -Djava.io.tmpdir=${tmpDir} -Xmx${memoSmall}g -jar ${GATK} -R $fasta
    -T VariantFiltration 
    -V ${output}_chr${chr}_indels.vcf.gz 
    --filterExpression \"QD < 2.0 || FS > 50.0 || ReadPosRankSum < -20.0\" 
    --filterName \"FAIL\" 
    -o ${output}_chr${chr}_indels_filtered.vcf.gz
```

Extract the SNPs:
```
$java  -Djava.io.tmpdir=${tmpDir} -Xmx${memoSmall}g -jar ${GATK} -R $fasta
     -T SelectVariants
     -V ${output}_chr${chr}.vcf.gz
     -selectType SNP 
     -o ${output}_chr${chr}_SNPs.vcf.gz
```

First SNPs:
```
$java -Djava.io.tmpdir=${tmpDir} -Xmx${memoSmall}g -jar ${GATK} -R $fasta 
   -T VariantRecalibrator
   -L $chr --input ${output}_chr${chr}_SNPs.vcf.gz --maxGaussians ${maxGauLoc} --mode SNP \
   -resource:hapmap,VCF,known=false,training=true,truth=true,prior=15.0 ${bundle}/hapmap_3.3.b37.vcf  \
   -resource:omni,VCF,known=false,training=true,truth=false,prior=12.0 ${bundle}/1000G_omni2.5.b37.vcf \
   -resource:dbsnp,VCF,known=true,training=false,truth=false,prior=8.0 ${bundle}/dbsnp_137.b37.vcf \
   -an QD -an FS -an ReadPosRankSum -an InbreedingCoeff \
   -tranche 100.0 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
   --minNumBadVariants ${numBad} \
   -recalFile ${output}_chr${chr}_SNPs_combrec.recal \
   -tranchesFile ${output}_chr${chr}_SNPs_combtranch \
   -rscriptFile  ${output}_chr${chr}_recal_plots_snps.R
```
```
${Rscript} ${output}_chr${chr}_recal_plots_snps.R
```

Apply recal:
```
$java -Xmx${memoSmall}g -jar ${GATK}  -R $fasta \
   -T ApplyRecalibration
   -o ${output}_chr${chr}_SNPs_filtered.vcf.gz \
   --ts_filter_level 99.5 \
   --recal_file ${output}_chr${chr}_SNPs_combrec.recal --tranches_file ${output}_chr${chr}_SNPs_combtranch --mode SNP \
   --input ${output}_chr${chr}_SNPs.vcf.gz
```

Now we merge SNPs and indels:
```
$java -Djava.io.tmpdir=${tmpDir} -Xmx${memoSmall}g -jar ${GATK} -R $fasta \
   -T CombineVariants --assumeIdenticalSamples \
   --variant:SNPs ${output}_chr${chr}_SNPs_filtered.vcf.gz \
   --variant:indels ${output}_chr${chr}_indels_filtered.vcf.gz \
   -genotypeMergeOptions PRIORITIZE  \
   -priority SNPs,indels \
   -o ${output}_chr${chr}_filtered.vcf
```

## Step 3: Annotation and custom filtering

## Step 4: Convert to R snpStats

