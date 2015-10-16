
# Overall pipeline

All looking great

## combine gVCF files (formerly step1)

combinegVCF.sh

## Massive joint calling (formerly step2)

msample_calling.sh

## Extract cohort specific information

crunch_data.R

## Finalize the user friendly folders

process_multiVCF.R




# Samples exclusion

sampleExclusion=/cluster/project8/vyp/exome_sequencing_multisamples/mainset/support/exclusion_lists/control_exclusion.tab

# Merge of VCFs

```
head -300 ${output}_chr1 | awk '{ if (\$1 ~ /^#/) print}' > ${output}.vcf
" >> $mainScript

    for chr in `seq 1 22` X ; do

	if [ ! -e ${output}_chr${chr} ]; then
	    echo "Missing file ${output}_chr${chr}"
	else
	    echo "awk '{ if (\$1 !~ /^#/) print}' ${output}_chr${chr} >> ${output}.vcf" >> $mainScript
	fi
    done

fi
```



# first SNPs

```
$java -Djava.io.tmpdir=${tmpDir} -Xmx${memoSmall}g -jar ${GATK} -T VariantRecalibrator -R $fasta --input ${output}.vcf --maxGaussians ${maxGaussians} --mode SNP \
             -resource:hapmap,VCF,known=false,training=true,truth=true,prior=15.0 ${bundle}/hapmap_3.3.b37.vcf  \
             -resource:omni,VCF,known=false,training=true,truth=false,prior=12.0 ${bundle}/1000G_omni2.5.b37.vcf \
             -resource:dbsnp,VCF,known=true,training=false,truth=false,prior=8.0 ${bundle}/dbsnp_137.b37.vcf \
             -an QD -an FS -an ReadPosRankSum -an InbreedingCoeff \
             -tranche 100.0 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
             --minNumBadVariants ${numBad} \
             -recalFile ${output}_SNP_combrec \
             -tranchesFile ${output}_SNP_combtranch \
             -rscriptFile ${output}_recal_plots_snps.R
```

```
${Rscript} ${output}_recal_plots_snps.R
```

```
# apply_recal
$java -Xmx${memoSmall}g -jar ${GATK} -T ApplyRecalibration -R $fasta \
       -o ${output}_recal_SNP.vcf \
       --ts_filter_level 99.5 \
       --recal_file ${output}_SNP_combrec --tranches_file ${output}_SNP_combtranch --mode SNP \
       --input ${output}.vcf \
" >> $mainScript

    
    if [[ "$indelHard" == "yes" ]]; then
	
	echo "

# extract the indels
$java  -Djava.io.tmpdir=${tmpDir} -Xmx${memoSmall}g -jar ${GATK} \
     -T SelectVariants \
     -R $fasta \
     -V ${output}_recal_SNP.vcf \
     -selectType INDEL \
     -selectType MIXED \
     -o ${output}_raw_indels.vcf

# apply the filters for the indels
$java -jar ${GATK} \
    -T VariantFiltration \
    -R $fasta \
    -V ${output}_raw_indels.vcf \
    --filterExpression \"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0\" \
    --filterName \"FAIL\" \
    -o ${output}_raw_indels_filtered.vcf


# extract the SNPs
$java  -Djava.io.tmpdir=${tmpDir} -Xmx${memoSmall}g -jar ${GATK} \
     -T SelectVariants \
     -R $fasta \
     -V ${output}_recal_SNP.vcf \
     -selectType SNP \
     -o ${output}_raw_SNPs.vcf


$java -Djava.io.tmpdir=${tmpDir} -Xmx${memoSmall}g -jar ${GATK} \
       -T CombineVariants --assumeIdenticalSamples \
       -R $fasta \
       -V ${output}_raw_SNPs.vcf \
       -V ${output}_raw_indels_filtered.vcf \
       -o ${output}_recal.vcf
```


# Variant recalibration

## then indels

```
$java -Djava.io.tmpdir=${tmpDir} -Xmx${memoSmall}g -jar ${GATK} -T VariantRecalibrator -R $fasta --input ${output}_recal_SNP.vcf --mode INDEL \
           -resource:mills,known=true,training=true,truth=true,prior=12.0 ${bundle}/Mills_and_1000G_gold_standard.indels.b37.vcf \
           -an QD -an FS -an ReadPosRankSum -an InbreedingCoeff \
           -tranche 100.0 -tranche 99.5  -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
           --minNumBadVariants ${numBadIndels} \
           --maxGaussians ${maxGaussiansIndels} \
           -recalFile ${output}_INDEL_combrec \
           -tranchesFile ${output}_INDEL_combtranch \
           -rscriptFile ${output}_recal_plots_indels.R

${Rscript} ${output}_recal_plots_indels.R

#apply_recal
$java -Xmx${memoSmall}g -jar ${GATK} -T ApplyRecalibration -R $fasta \
         --input ${output}_recal_SNP.vcf --out ${output}_recal.vcf \
         --recal_file ${output}_INDEL_combrec --tranches_file ${output}_INDEL_combtranch --mode INDEL \
         --ts_filter_level 95.0

rm -rf $tmpDir

" >> $mainScript
    fi
fi
```

# Annovar


if [[ "$annovar" == "yes" ]]; then

    echo "

cut -f1-8 ${output}_recal.vcf > ${output}_for_annovar.vcf

/cluster/project8/vyp/vincent/Software_heavy/annovar_Feb2013/convert2annovar.pl --allallele -format vcf4 --includeinfo ${output}_for_annovar.vcf > ${output}_db

/cluster/project8/vyp/vincent/Software_heavy/annovar_Feb2013/summarize_annovar_VP.pl -ver1000g 1000g2012apr -verdbsnp 137 -veresp 6500si -alltranscript -buildver hg19 --genetype ensgene --remove ${output}_db /cluster/project8/vyp/vincent/Software_heavy/annovar_Feb2013/humandb_hg19/


" >> $mainScript

fi


if [[ "$filter" == "yes" ]]; then

    echo "

perl ~/Software/pipeline/GATK_v2/custom_filtering.pl ${output}_recal.vcf ${output}_recal_filtered.vcf ${GQ}

python /cluster/project8/vyp/vincent/Software/pipeline/GATK_v2/annovar_vcf_combine_VP.py ${output}_recal_filtered.vcf ${output}_db.exome_summary.csv ${output}_exome_table.csv

" >> $mainScript

fi



if [[ "$exomeprep" == "yes" ]]; then

    keyWords=data/controlKeywords.tab
    casekeyWords=data/caseKeywords.tab
    
    if [ ! -e "${output}_split_by_chr" ]; then mkdir ${output}_by_chr; fi
    
    echo "

# Make matrix calls

perl /cluster/project8/vyp/vincent/Software/pipeline/GATK_v2/make_matrix_calls.pl ${output}_exome_table.csv ${output}

$Rbin CMD BATCH --no-save --no-restore --root=${output} /cluster/project8/vyp/vincent/Software/pipeline/GATK_v2/split_data_by_chromosomes.R cluster/R/huge.out

# PCA analysis

#$Rbin CMD BATCH --no-save --no-restore scripts/PCA/PCA_analysis.R cluster/R/PCA.out



if [[ "$finalCrunch" == "yes" ]]; then

    keyWords=data/controlKeywords.tab
    casekeyWords=data/caseKeywords.tab

    echo "

# Crunch 

/cluster/project8/vyp/vincent/Software/pipeline/GATK_v2/crunch_controls.pl

${crunchpl} ${output}_exome_table.csv $keyWords $casekeyWords ${output}_exome_crunched.csv data/sampleList_exome.tab none no  ##include all samples

# Annovar VCF combine

python /cluster/project8/vyp/vincent/Software/pipeline/GATK_v2/annovar_vcf_combine_VP.py ${output}_recal_filtered.vcf ${output}_db.genome_summary.csv ${output}_genome_table.csv

# Crunch controls

/cluster/project8/vyp/vincent/Software/pipeline/GATK_v2/crunch_controls.pl

${crunchpl} ${output}_genome_table.csv $keyWords $casekeyWords ${output}_genome_crunched.csv data/sampleList_genome.tab none no  ##include all samples
