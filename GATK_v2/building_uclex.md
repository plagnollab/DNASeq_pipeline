
# Location of UCLex output
```
/SAN/vyplab/UCLex/
```

# Current working directory
```
/cluster/project8/vyp/exome_sequencing_multisamples/mainset
```

# Main joint calling script
```
/cluster/project8/vyp/exome_sequencing_multisamples/mainset/cluster/submission/calling.sh
```

# Calls GenotypeGVCFs
```
/share/apps/jdk/jre/bin/java -Djava.io.tmpdir=/scratch0/ -Xmx5g -jar /cluster/project8/vyp/vincent/Software/GenomeAnalysisTK-3.5-0/GenomeAnalysisTK.jar \
   -R /cluster/scratch3/vyp-scratch2/reference_datasets/human_reference_sequence/human_g1k_v37.fasta \
   -T GenotypeGVCFs \
   -L 1 -L /cluster/project8/vyp/exome_sequencing_multisamples/target_region/data/merged_exome_target_cleaned.bed --interval_set_rule INTERSECTION --interval_padding 100  \
   --annotation InbreedingCoeff --annotation QualByDepth --annotation HaplotypeScore --annotation MappingQualityRankSumTest --annotation ReadPosRankSumTest --annotation FisherStrand \
   --dbsnp /cluster/scratch3/vyp-scratch2/reference_datasets/GATK_bundle/dbsnp_137.b37.vcf \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr1/Syrris_ARVC_exomes_1.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr1/Prion_batch2_1.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr1/Prion_batch2_5.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr1/Kelsell_1.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr1/Prion_batch2_2.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr1/Prion_batch2_3.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr1/Prion_batch2_7.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr1/Prion_batch2_6.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr1/Prion_batch2_12.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr1/Prion_batch2_9.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr1/Prion_batch2_8.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr1/Prion_batch2_11.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr1/Prion_batch2_10.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr1/Prion_batch2_4.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr1/Sisodiya_1.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr1/Sisodiya_4.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr1/Sisodiya_3.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr1/Sisodiya_2.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr1/Shamima_1.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr1/Vulliamy_1.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr1/Vulliamy_2.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr1/Vulliamy2_1.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr1/SPEED_2.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr1/SPEED_1.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr1/Prion_batch1_1.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr1/Prion_batch1_2.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr1/Prion_batch1_3.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr1/IoN_1.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr1/Davina_BAMs_1.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr1/1KG_4.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr1/1KG_3.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr1/1KG_2.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr1/1KG_1.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr1/UK10K_1.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr1/AdamLevine_1.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr1/GOSgene_4.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr1/GOSgene_1.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr1/GOSgene_2.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr1/Hardcastle_1.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr1/IoO_1.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr1/GOSgene_3.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr1/Levine_Aug2014_2.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr1/Levine_Aug2014_1.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr1/Levine_Nov2015_1.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr1/Lambiase_1.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr1/FFSIOO_1.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr1/Nejentsev_2.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr1/Nejentsev_1.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr1/Nejentsev2_1.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr1/RPFB_1.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr1/IRDC_1.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr1/IRDCg2_1.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr1/mcQuillin_1.gvcf.gz \
   --variant /cluster/project9/IBDAJE/batches_for_uclex/1410AHP_0004_115samples/Macrogen-1410AHP_0004-AJ-115-chr1.gvcf.gz \
   --variant /cluster/project9/IBDAJE/batches_for_uclex/1410AHP_0004_96samples/Macrogen-1410AHP_0004-AJ-96-chr1.gvcf.gz \
   --variant /cluster/project9/IBDAJE/batches_for_uclex/1411AHP_0005_67samples/Macrogen-1411AHP_0005-AJ-67-chr1.gvcf.gz \
   --variant /cluster/project9/IBDAJE/batches_for_uclex/1411AHP_0005_75samples/Macrogen-1411AHP_0005-AJ-75-chr1.gvcf.gz \
   --variant /cluster/project9/IBDAJE/batches_for_uclex/ICE_38samples/Macrogen-1411AHP_0005-ICE-38-chr1.gvcf.gz \
   --variant /cluster/project9/IBDAJE/batches_for_uclex/OFG_29samples/Macrogen-1411AHP_0005-OFG-29-chr1.gvcf.gz \
   --variant /cluster/project9/IBDAJE/batches_for_uclex/1411KHP_0031_68samples/Macrogen-1411KHP_0031-AJ-68-chr1.gvcf.gz \
   --variant /cluster/project9/IBDAJE/batches_for_uclex/1412KHP_0015_84samples/Macrogen-1412KHP_0015-AJ-84-chr1.gvcf.gz \
   --variant /cluster/project9/IBDAJE/batches_for_uclex/BGI_Nov2014_14samples/BGI-Nov2014-AJ-14-chr1.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/RoseRichardson/batches_for_uclex/chr1.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/pontikos/Hardcastle_Feb2016/chr1.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/exomes_temp2/IRDC/IRDC_batch7/IRDC/IRDC_batch7/data/chr1.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/exomes_temp2/IRDC/IRDC_batch1/IRDC/IRDC_batch1/data/chr1.gvcf.gz \
   --variant /cluster/project9/IBDAJE/ExomeSequences/BGI/BGI_April2016_88samples/BGI_April2016_88samples/data/chr1.gvcf.gz \
   -o /SAN/vyplab/UCLex/mainset_June2016/mainset_June2016_chr1.vcf.gz
```

# Extract the indels
```
/share/apps/jdk/jre/bin/java  -Djava.io.tmpdir=/scratch0/GATK_chr22 -Xmx5g -jar /cluster/project8/vyp/vincent/Software/GenomeAnalysisTK-3.5-0/GenomeAnalysisTK.jar      -T SelectVariants      -R /cluster/scratch3/vyp-scratch2/reference_datasets/human_reference_sequence/human_g1k_v37.fasta      -V /SAN/vyplab/UCLex/mainset_June2016/mainset_June2016_chr22.vcf.gz      -selectType INDEL      -selectType MIXED      -o /SAN/vyplab/UCLex/mainset_June2016/mainset_June2016_chr22_indels.vcf.gz
```

# Apply the filters for the indels
```
/share/apps/jdk/jre/bin/java -Djava.io.tmpdir=/scratch0/GATK_chr22 -Xmx5g -jar /cluster/project8/vyp/vincent/Software/GenomeAnalysisTK-3.5-0/GenomeAnalysisTK.jar     -T VariantFiltration     -R /cluster/scratch3/vyp-scratch2/reference_datasets/human_reference_sequence/human_g1k_v37.fasta     -V /SAN/vyplab/UCLex/mainset_June2016/mainset_June2016_chr22_indels.vcf.gz     --filterExpression "QD < 2.0 || FS > 50.0 || ReadPosRankSum < -20.0"     --filterName "FAIL"     -o /SAN/vyplab/UCLex/mainset_June2016/mainset_June2016_chr22_indels_filtered.vcf.gz
```

# Extract the SNPs
```
/share/apps/jdk/jre/bin/java  -Djava.io.tmpdir=/scratch0/GATK_chr22 -Xmx5g -jar /cluster/project8/vyp/vincent/Software/GenomeAnalysisTK-3.5-0/GenomeAnalysisTK.jar      -T SelectVariants      -R /cluster/scratch3/vyp-scratch2/reference_datasets/human_reference_sequence/human_g1k_v37.fasta      -V /SAN/vyplab/UCLex/mainset_June2016/mainset_June2016_chr22.vcf.gz      -selectType SNP      -o /SAN/vyplab/UCLex/mainset_June2016/mainset_June2016_chr22_SNPs.vcf.gz
```

# Variant recalibrator of SNPs
```
/share/apps/jdk/jre/bin/java -Djava.io.tmpdir=/scratch0/GATK_chr22 -Xmx5g -jar /cluster/project8/vyp/vincent/Software/GenomeAnalysisTK-3.5-0/GenomeAnalysisTK.jar -T VariantRecalibrator -R /cluster/scratch3/vyp-scratch2/reference_datasets/human_reference_sequence/human_g1k_v37.fasta -L 22
--input /SAN/vyplab/UCLex/mainset_June2016/mainset_June2016_chr22_SNPs.vcf.gz
--maxGaussians 6 --mode SNP
-resource:hapmap,VCF,known=false,training=true,truth=true,prior=15.0 /cluster/scratch3/vyp-scratch2/reference_datasets/GATK_bundle/hapmap_3.3.b37.vcf               -resource:omni,VCF,known=false,training=true,truth=false,prior=12.0 /cluster/scratch3/vyp-scratch2/reference_datasets/GATK_bundle/1000G_omni2.5.b37.vcf              -resource:dbsnp,VCF,known=true,training=false,truth=false,prior=8.0 /cluster/scratch3/vyp-scratch2/reference_datasets/GATK_bundle/dbsnp_137.b37.vcf              -an QD -an FS -an ReadPosRankSum -an InbreedingCoeff              -tranche 100.0 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0              --minNumBadVariants 1000              -recalFile /SAN/vyplab/UCLex/mainset_June2016/mainset_June2016_chr22_SNPs_combrec.recal              -tranchesFile /SAN/vyplab/UCLex/mainset_June2016/mainset_June2016_chr22_SNPs_combtranch              -rscriptFile  /SAN/vyplab/UCLex/mainset_June2016/mainset_June2016_chr22_recal_plots_snps.R

/share/apps/R-3.2.2/bin/Rscript /SAN/vyplab/UCLex/mainset_June2016/mainset_June2016_chr22_recal_plots_snps.R
```

# Applying recalibration
```
/share/apps/jdk/jre/bin/java -Xmx5g -jar /cluster/project8/vyp/vincent/Software/GenomeAnalysisTK-3.5-0/GenomeAnalysisTK.jar -T ApplyRecalibration -R /cluster/scratch3/vyp-scratch2/reference_datasets/human_reference_sequence/human_g1k_v37.fasta        -o /SAN/vyplab/UCLex/mainset_June2016/mainset_June2016_chr22_SNPs_filtered.vcf.gz        --ts_filter_level 99.5        --recal_file /SAN/vyplab/UCLex/mainset_June2016/mainset_June2016_chr22_SNPs_combrec.recal --tranches_file /SAN/vyplab/UCLex/mainset_June2016/mainset_June2016_chr22_SNPs_combtranch --mode SNP        --input /SAN/vyplab/UCLex/mainset_June2016/mainset_June2016_chr22_SNPs.vcf.gz
```

Now we merge SNPs and indels:
```
/share/apps/jdk/jre/bin/java -Djava.io.tmpdir=/scratch0/GATK_chr22 -Xmx5g -jar /cluster/project8/vyp/vincent/Software/GenomeAnalysisTK-3.5-0/GenomeAnalysisTK.jar
-T CombineVariants --assumeIdenticalSamples        -R /cluster/scratch3/vyp-scratch2/reference_datasets/human_reference_sequence/human_g1k_v37.fasta        --variant:SNPs /SAN/vyplab/UCLex/mainset_June2016/mainset_June2016_chr22_SNPs_filtered.vcf.gz        --variant:indels /SAN/vyplab/UCLex/mainset_June2016/mainset_June2016_chr22_indels_filtered.vcf.gz        -genotypeMergeOptions PRIORITIZE         -priority SNPs,indels        -o /SAN/vyplab/UCLex/mainset_June2016/mainset_June2016_chr22_filtered.vcf
```

# Removing temp files
```
rm -rf /scratch0/GATK_chr22
rm /SAN/vyplab/UCLex/mainset_June2016/mainset_June2016_chr22_indels.vcf.gz /SAN/vyplab/UCLex/mainset_June2016/mainset_June2016_chr22_SNPs.vcf.gz /SAN/vyplab/UCLex/mainset_June2016/mainset_June2016_chr22_SNPs_filtered.vcf.gz
```

```
if [ -e /SAN/vyplab/UCLex/mainset_June2016/mainset_June2016_snpStats/chr21.done ]; then rm /SAN/vyplab/UCLex/mainset_June2016/mainset_June2016_snpStats/chr21.done; fi  ## this is basically a log file, to make sure the job got finished
```
Only keep the first 8 columns: CHR, POS, ID, REF, ALT, INFO, FILTER...
```
cut -f1-8 /SAN/vyplab/UCLex/mainset_June2016/mainset_June2016_chr21_filtered.vcf > /SAN/vyplab/UCLex/mainset_June2016/mainset_June2016_chr21_for_annovar.vcf
```

# Running Annovar

```
/cluster/project8/vyp/vincent/Software_heavy/annovar_Feb2013/convert2annovar.pl --allallele -format vcf4 --includeinfo /SAN/vyplab/UCLex/mainset_June2016/mainset_June2016_chr21_for_annovar.vcf > /SAN/vyplab/UCLex/mainset_June2016/mainset_June2016_chr21_db

/cluster/project8/vyp/vincent/Software_heavy/annovar_Feb2013/summarize_annovar_VP.pl -ver1000g 1000g2012apr -verdbsnp 137 -veresp 6500si -alltranscript -buildver hg19 --genetype ensgene --remove /SAN/vyplab/UCLex/mainset_June2016/mainset_June2016_chr21_db /cluster/project8/vyp/vincent/Software_heavy/annovar_Feb2013/humandb_hg19/
```

```
perl /cluster/project8/vyp/vincent/Software/DNASeq_pipeline/GATK_v2/custom_filtering.pl /SAN/vyplab/UCLex/mainset_June2016/mainset_June2016_chr21_filtered.vcf /SAN/vyplab/UCLex/mainset_June2016/mainset_June2016_chr21_recal_filtered2.vcf 20

python /cluster/project8/vyp/vincent/Software/DNASeq_pipeline/GATK_v2/annovar_vcf_combine_VP.py /SAN/vyplab/UCLex/mainset_June2016/mainset_June2016_chr21_recal_filtered2.vcf /SAN/vyplab/UCLex/mainset_June2016/mainset_June2016_chr21_db.exome_summary.csv /SAN/vyplab/UCLex/mainset_June2016/mainset_June2016_chr21_exome_table.csv

perl /cluster/project8/vyp/vincent/Software/DNASeq_pipeline/msample_calling/make_matrix_calls.pl /SAN/vyplab/UCLex/mainset_June2016/mainset_June2016_chr21_exome_table.csv /SAN/vyplab/UCLex/mainset_June2016/mainset_June2016 21

touch /SAN/vyplab/UCLex/mainset_June2016/mainset_June2016_snpStats/chr21.done  ##here we mark that the scripts finished
```

# Building the variant lists

```
cat /cluster/project8/vyp/exome_sequencing_multisamples/IRDC/cluster/submission/IRDC.sh

#!/bin/bash
#$ -S /bin/bash
#$ -o cluster/out
#$ -e cluster/error
#$ -cwd
#$ -pe smp 1
#$ -l tmem=12.9G,h_vmem=12.9G
#$ -l h_rt=20:0:0
#$ -V
#$ -R y

/cluster/project8/vyp/vincent/Software/R-3.3.0/bin/R CMD BATCH --no-save --no-restore --root=/SAN/vyplab/UCLex/mainset_June2016/mainset_June2016 --minDepth=7 --Prion.setup=FALSE /cluster/project8/vyp/vincent/Software/DNASeq_pipeline/GATK_v2/crunch_data.R cluster/R/crunch.out

/cluster/project8/vyp/vincent/Software/R-3.3.0/bin/R CMD BATCH --no-save --no-restore scripts/merge_GATK.R cluster/R/GATK.out
```

```
cat /cluster/project8/vyp/exome_sequencing_multisamples/Vulliamy/cluster/submission/Vulliamy.sh

#!/bin/bash
#$ -S /bin/bash
#$ -o cluster/out
#$ -e cluster/error
#$ -cwd
#$ -pe smp 1
#$ -l tmem=12.9G,h_vmem=12.9G
#$ -l h_rt=20:0:0
#$ -V
#$ -R y

/cluster/project8/vyp/vincent/Software/R-3.3.0/bin/R CMD BATCH --no-save --no-restore --root=/SAN/vyplab/UCLex/mainset_June2016/mainset_June2016 --minDepth=7 --Prion.setup=FALSE /cluster/project8/vyp/vincent/Software/DNASeq_pipeline/GATK_v2/crunch_data.R cluster/R/crunch.out

/cluster/project8/vyp/vincent/Software/R-3.3.0/bin/R CMD BATCH --no-save --no-restore scripts/merge_GATK.R cluster/R/GATK.out
```
