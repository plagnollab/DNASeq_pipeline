
Current working directory:
```
/cluster/project8/vyp/exome_sequencing_multisamples/mainset
```

Joint calling:
```
/cluster/project8/vyp/exome_sequencing_multisamples/mainset/cluster/submission/calling.sh
```
calls GenotypeGVCFs:
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


