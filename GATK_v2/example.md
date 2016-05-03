
# GenotypeGVCFs
```
/share/apps/jdk/jre/bin/java -Xmx5g -jar /cluster/project8/vyp/vincent/Software/GenomeAnalysisTK-3.5-0/GenomeAnalysisTK.jar \
   -R /cluster/scratch3/vyp-scratch2/reference_datasets/human_reference_sequence/human_g1k_v37.fasta \
   -T GenotypeGVCFs \
   -L 22 -L /cluster/project8/vyp/exome_sequencing_multisamples/target_region/data/merged_exome_target_cleaned.bed --interval_set_rule INTERSECTION --interval_padding 100  \
   --annotation InbreedingCoeff --annotation QualByDepth --annotation HaplotypeScore --annotation MappingQualityRankSumTest --annotation ReadPosRankSumTest --annotation FisherStrand \
   --dbsnp /cluster/scratch3/vyp-scratch2/reference_datasets/GATK_bundle/dbsnp_137.b37.vcf \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr22/Syrris_ARVC_exomes_1.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr22/Prion_batch2_1.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr22/Prion_batch2_5.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr22/Kelsell_1.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr22/Prion_batch2_2.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr22/Prion_batch2_3.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr22/Prion_batch2_7.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr22/Prion_batch2_6.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr22/Prion_batch2_12.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr22/Prion_batch2_9.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr22/Prion_batch2_8.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr22/Prion_batch2_11.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr22/Prion_batch2_10.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr22/Prion_batch2_4.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr22/Sisodiya_1.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr22/Sisodiya_4.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr22/Sisodiya_3.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr22/Sisodiya_2.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr22/Shamima_1.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr22/Vulliamy_1.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr22/Vulliamy_2.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr22/Vulliamy2_1.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr22/SPEED_2.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr22/SPEED_1.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr22/Prion_batch1_1.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr22/Prion_batch1_2.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr22/Prion_batch1_3.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr22/IoN_1.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr22/Davina_BAMs_1.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr22/1KG_4.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr22/1KG_3.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr22/1KG_2.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr22/1KG_1.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr22/UK10K_1.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr22/AdamLevine_1.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr22/GOSgene_4.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr22/GOSgene_1.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr22/GOSgene_2.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr22/Hardcastle_1.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr22/IoO_1.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr22/GOSgene_3.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr22/Levine_Aug2014_2.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr22/Levine_Aug2014_1.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr22/Levine_Nov2015_1.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr22/Lambiase_1.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr22/FFSIOO_1.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr22/Nejentsev_2.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr22/Nejentsev_1.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr22/Nejentsev2_1.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr22/RPFB_1.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr22/IRDC_1.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/vincent/GATK/HC/combinedVCFs/chr22/mcQuillin_1.gvcf.gz \
   --variant /cluster/project8/IBDAJE/batches_for_uclex/1410AHP_0004_115samples/Macrogen-1410AHP_0004-AJ-115-chr22.gvcf.gz \
   --variant /cluster/project8/IBDAJE/batches_for_uclex/1410AHP_0004_96samples/Macrogen-1410AHP_0004-AJ-96-chr22.gvcf.gz \
   --variant /cluster/project8/IBDAJE/batches_for_uclex/1411AHP_0005_67samples/Macrogen-1411AHP_0005-AJ-67-chr22.gvcf.gz \
   --variant /cluster/project8/IBDAJE/batches_for_uclex/1411AHP_0005_75samples/Macrogen-1411AHP_0005-AJ-75-chr22.gvcf.gz \
   --variant /cluster/project8/IBDAJE/batches_for_uclex/ICE_38samples/Macrogen-1411AHP_0005-ICE-38-chr22.gvcf.gz \
   --variant /cluster/project8/IBDAJE/batches_for_uclex/OFG_29samples/Macrogen-1411AHP_0005-OFG-29-chr22.gvcf.gz \
   --variant /cluster/project8/IBDAJE/batches_for_uclex/1411KHP_0031_68samples/Macrogen-1411KHP_0031-AJ-68-chr22.gvcf.gz \
   --variant /cluster/project8/IBDAJE/batches_for_uclex/1412KHP_0015_84samples/Macrogen-1412KHP_0015-AJ-84-chr22.gvcf.gz \
   --variant /cluster/project8/IBDAJE/batches_for_uclex/BGI_Nov2014_14samples/BGI-Nov2014-AJ-14-chr22.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/RoseRichardson/batches_for_uclex/chr22.gvcf.gz \
   --variant /cluster/scratch3/vyp-scratch2/pontikos/Hardcastle_Feb2016//chr22.gvcf.gz \
   -o /cluster/scratch3/vyp-scratch2/vincent/GATK/mainset_February2016/mainset_February2016_chr22.vcf.gz
```

```
if [ ! -e /scratch0/GATK_chr22 ]; then mkdir /scratch0/GATK_chr22; fi
```

#### extract the indels
```
/share/apps/jdk1.7.0_45/bin/java  -Djava.io.tmpdir=/scratch0/GATK_chr22 -Xmx5g -jar /cluster/project8/vyp/vincent/Software/GenomeAnalysisTK-3.5-0/GenomeAnalysisTK.jar      -T SelectVariants      -R /cluster/scratch3/vyp-scratch2/reference_datasets/human_reference_sequence/human_g1k_v37.fasta      -V /cluster/scratch3/vyp-scratch2/vincent/GATK/mainset_February2016/mainset_February2016_chr22.vcf.gz      -selectType INDEL      -selectType MIXED      -o /cluster/scratch3/vyp-scratch2/vincent/GATK/mainset_February2016/mainset_February2016_chr22_indels.vcf.gz
```

#### apply the filters for the indels
```
/share/apps/jdk1.7.0_45/bin/java -Djava.io.tmpdir=/scratch0/GATK_chr22 -Xmx5g -jar /cluster/project8/vyp/vincent/Software/GenomeAnalysisTK-3.5-0/GenomeAnalysisTK.jar     -T VariantFiltration     -R /cluster/scratch3/vyp-scratch2/reference_datasets/human_reference_sequence/human_g1k_v37.fasta     -V /cluster/scratch3/vyp-scratch2/vincent/GATK/mainset_February2016/mainset_February2016_chr22_indels.vcf.gz     --filterExpression "QD < 2.0 || FS > 50.0 || ReadPosRankSum < -20.0"     --filterName "FAIL"     -o /cluster/scratch3/vyp-scratch2/vincent/GATK/mainset_February2016/mainset_February2016_chr22_indels_filtered.vcf.gz
```

#### extract the SNPs
```
/share/apps/jdk1.7.0_45/bin/java  -Djava.io.tmpdir=/scratch0/GATK_chr22 -Xmx5g -jar /cluster/project8/vyp/vincent/Software/GenomeAnalysisTK-3.5-0/GenomeAnalysisTK.jar      -T SelectVariants      -R /cluster/scratch3/vyp-scratch2/reference_datasets/human_reference_sequence/human_g1k_v37.fasta      -V /cluster/scratch3/vyp-scratch2/vincent/GATK/mainset_February2016/mainset_February2016_chr22.vcf.gz      -selectType SNP      -o /cluster/scratch3/vyp-scratch2/vincent/GATK/mainset_February2016/mainset_February2016_chr22_SNPs.vcf.gz
```


####### first SNPs
```
/share/apps/jdk1.7.0_45/bin/java -Djava.io.tmpdir=/scratch0/GATK_chr22 -Xmx5g -jar /cluster/project8/vyp/vincent/Software/GenomeAnalysisTK-3.5-0/GenomeAnalysisTK.jar -T VariantRecalibrator -R /cluster/scratch3/vyp-scratch2/reference_datasets/human_reference_sequence/human_g1k_v37.fasta -L 22 --input /cluster/scratch3/vyp-scratch2/vincent/GATK/mainset_February2016/mainset_February2016_chr22_SNPs.vcf.gz --maxGaussians 6 --mode SNP              -resource:hapmap,VCF,known=false,training=true,truth=true,prior=15.0 /cluster/scratch3/vyp-scratch2/reference_datasets/GATK_bundle/hapmap_3.3.b37.vcf               -resource:omni,VCF,known=false,training=true,truth=false,prior=12.0 /cluster/scratch3/vyp-scratch2/reference_datasets/GATK_bundle/1000G_omni2.5.b37.vcf              -resource:dbsnp,VCF,known=true,training=false,truth=false,prior=8.0 /cluster/scratch3/vyp-scratch2/reference_datasets/GATK_bundle/dbsnp_137.b37.vcf              -an QD -an FS -an ReadPosRankSum -an InbreedingCoeff              -tranche 100.0 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0              --minNumBadVariants 1000              -recalFile /cluster/scratch3/vyp-scratch2/vincent/GATK/mainset_February2016/mainset_February2016_chr22_SNPs_combrec.recal              -tranchesFile /cluster/scratch3/vyp-scratch2/vincent/GATK/mainset_February2016/mainset_February2016_chr22_SNPs_combtranch              -rscriptFile  /cluster/scratch3/vyp-scratch2/vincent/GATK/mainset_February2016/mainset_February2016_chr22_recal_plots_snps.R
/cluster/project8/vyp/vincent/Software/R-3.1.2/bin/Rscript /cluster/scratch3/vyp-scratch2/vincent/GATK/mainset_February2016/mainset_February2016_chr22_recal_plots_snps.R
```

# apply_recal
```
/share/apps/jdk1.7.0_45/bin/java -Xmx5g -jar /cluster/project8/vyp/vincent/Software/GenomeAnalysisTK-3.5-0/GenomeAnalysisTK.jar -T ApplyRecalibration -R /cluster/scratch3/vyp-scratch2/reference_datasets/human_reference_sequence/human_g1k_v37.fasta        -o /cluster/scratch3/vyp-scratch2/vincent/GATK/mainset_February2016/mainset_February2016_chr22_SNPs_filtered.vcf.gz        --ts_filter_level 99.5        --recal_file /cluster/scratch3/vyp-scratch2/vincent/GATK/mainset_February2016/mainset_February2016_chr22_SNPs_combrec.recal --tranches_file /cluster/scratch3/vyp-scratch2/vincent/GATK/mainset_February2016/mainset_February2016_chr22_SNPs_combtranch --mode SNP        --input /cluster/scratch3/vyp-scratch2/vincent/GATK/mainset_February2016/mainset_February2016_chr22_SNPs.vcf.gz
```

#### Now we merge SNPs and indels
```
/share/apps/jdk1.7.0_45/bin/java -Djava.io.tmpdir=/scratch0/GATK_chr22 -Xmx5g -jar /cluster/project8/vyp/vincent/Software/GenomeAnalysisTK-3.5-0/GenomeAnalysisTK.jar        -T CombineVariants --assumeIdenticalSamples        -R /cluster/scratch3/vyp-scratch2/reference_datasets/human_reference_sequence/human_g1k_v37.fasta        --variant:SNPs /cluster/scratch3/vyp-scratch2/vincent/GATK/mainset_February2016/mainset_February2016_chr22_SNPs_filtered.vcf.gz        --variant:indels /cluster/scratch3/vyp-scratch2/vincent/GATK/mainset_February2016/mainset_February2016_chr22_indels_filtered.vcf.gz        -genotypeMergeOptions PRIORITIZE         -priority SNPs,indels        -o /cluster/scratch3/vyp-scratch2/vincent/GATK/mainset_February2016/mainset_February2016_chr22_filtered.vcf
rm -rf /scratch0/GATK_chr22
rm /cluster/scratch3/vyp-scratch2/vincent/GATK/mainset_February2016/mainset_February2016_chr22_indels.vcf.gz /cluster/scratch3/vyp-scratch2/vincent/GATK/mainset_February2016/mainset_February2016_chr22_SNPs.vcf.gz /cluster/scratch3/vyp-scratch2/vincent/GATK/mainset_February2016/mainset_February2016_chr22_SNPs_filtered.vcf.gz
```

## this is basically a log file, to make sure the job got finished
```
if [ -e /cluster/scratch3/vyp-scratch2/vincent/GATK/mainset_February2016/mainset_February2016_snpStats/chr22.done ]; then rm /cluster/scratch3/vyp-scratch2/vincent/GATK/mainset_February2016/mainset_February2016_snpStats/chr22.done; fi

cut -f1-8 /cluster/scratch3/vyp-scratch2/vincent/GATK/mainset_February2016/mainset_February2016_chr22_filtered.vcf > /cluster/scratch3/vyp-scratch2/vincent/GATK/mainset_February2016/mainset_February2016_chr22_for_annovar.vcf

/cluster/project8/vyp/vincent/Software_heavy/annovar_Feb2013/convert2annovar.pl --allallele -format vcf4 --includeinfo /cluster/scratch3/vyp-scratch2/vincent/GATK/mainset_February2016/mainset_February2016_chr22_for_annovar.vcf > /cluster/scratch3/vyp-scratch2/vincent/GATK/mainset_February2016/mainset_February2016_chr22_db

/cluster/project8/vyp/vincent/Software_heavy/annovar_Feb2013/summarize_annovar_VP.pl -ver1000g 1000g2012apr -verdbsnp 137 -veresp 6500si -alltranscript -buildver hg19 --genetype ensgene --remove /cluster/scratch3/vyp-scratch2/vincent/GATK/mainset_February2016/mainset_February2016_chr22_db /cluster/project8/vyp/vincent/Software_heavy/annovar_Feb2013/humandb_hg19/

perl /cluster/project8/vyp/vincent/Software/DNASeq_pipeline/GATK_v2/custom_filtering.pl /cluster/scratch3/vyp-scratch2/vincent/GATK/mainset_February2016/mainset_February2016_chr22_filtered.vcf /cluster/scratch3/vyp-scratch2/vincent/GATK/mainset_February2016/mainset_February2016_chr22_recal_filtered2.vcf 20

python /cluster/project8/vyp/vincent/Software/DNASeq_pipeline/GATK_v2/annovar_vcf_combine_VP.py /cluster/scratch3/vyp-scratch2/vincent/GATK/mainset_February2016/mainset_February2016_chr22_recal_filtered2.vcf /cluster/scratch3/vyp-scratch2/vincent/GATK/mainset_February2016/mainset_February2016_chr22_db.exome_summary.csv /cluster/scratch3/vyp-scratch2/vincent/GATK/mainset_February2016/mainset_February2016_chr22_exome_table.csv

perl /cluster/project8/vyp/vincent/Software/DNASeq_pipeline/msample_calling/make_matrix_calls.pl /cluster/scratch3/vyp-scratch2/vincent/GATK/mainset_February2016/mainset_February2016_chr22_exome_table.csv /cluster/scratch3/vyp-scratch2/vincent/GATK/mainset_February2016/mainset_February2016 22

touch /cluster/scratch3/vyp-scratch2/vincent/GATK/mainset_February2016/mainset_February2016_snpStats/chr22.done  ##here we mark that the scripts finished
```

