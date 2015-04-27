# Variant annotation


## Summary
The annotation of both coding and non-coding variants is done using the [Variant Effect Predictor](http://www.ensembl.org/info/docs/tools/vep/index.html)(VEP).
All variants are annotated with [1000 Genomes](www.1000genomes.org) allele frequencies.
Additionally, coding variants are annotated with allele frequencies from [ExAC](http://exac.broadinstitute.org/).
Variants are further annotated with Combined Annotation Dependent Depletion (CADD)<sup>[1][CADD]</sup>, Combined Annotation scoRing toOL (CAROL) and CONsensus DELeteriousness (Condel) consequence scores, to assess their functional impact.

## Example

```
#!/bin/bash
#$ -S /bin/bash
#$ -o /dev/null
#$ -e /dev/null
#$ -cwd
#$ -V
#$ -l tmem=6G,h_vmem=6G
#$ -l h_rt=24:0:0
#$ -t 1-25
set -u
set -x
scriptname=annotate
mkdir -p ${scriptname}.qsub.out ${scriptname}.qsub.err
exec >${scriptname}.qsub.out/${scriptname}_${SGE_TASK_ID}_${JOB_ID}.out 2>${scriptname}.qsub.err/${scriptname}_${SGE_TASK_ID}_${JOB_ID}.err
args=( header 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y )
f=${args[$SGE_TASK_ID]}
chrCode=$f


function VEP() {
    DIR=~/bin/DNASeq_pipeline/annotation
    reference=1kg
    zcat chr${chrCode}.vcf.gz | python ${DIR}/multiallele_to_single_gvcf.py --GQ 10 --DP 5 > chr${chrCode}-single.vcf
    bash $DIR/run_VEP.sh --vcfin chr${chrCode}-single.vcf --chr $chrCode --reference $reference --vcfout VEP_${chrCode}.vcfout --coding_only yes --custom UCLEX,AJcontrols,AJcases,BroadAJcontrols,ImmunoBase
}


function annotate() {
    DIR=~/bin/DNASeq_pipeline/annotation
    python $DIR/processVEP.py --custom-allele-freq UCLEX AJcontrols AJcases BroadAJcontrols --custom-annotation ImmunoBase_CRO ImmunoBase_IBD ImmunoBase_UC --file VEP_${chrCode}.vcfout
    #combine annotation, genotype and depth into a single VEP file.
    Rscript $DIR/combine.R --chr $chrCode | Rscript $DIR/annotate/annotate_expression.R | Rscript $DIR/annotate/annotate_pedigree_af.R| Rscript scripts/Group_Counts.R > VEP_${chrCode}.csv
}

function filter() {
    DIR=~/bin/DNASeq_pipeline/annotation/filters/
    cat VEP_$chrCode.csv | Rscript $DIR/af-filter.R --ajcontrols.thresh 0.05 --uclex.thresh 0.05 --exac.thresh 0.025 > af-filtered-VEP_$chrCode.csv
    wc -l  af-filtered-VEP_$chrCode.csv
    cat af-filtered-VEP_$chrCode.csv |  Rscript $DIR/csq-filter.R > af-csq-filtered-VEP_$chrCode.csv
    wc -l af-csq-filtered-VEP_$chrCode.csv
    cat af-csq-filtered-VEP_$chrCode.csv | Rscript $DIR/GO-filter.R > af-csq-go-filtered-VEP_$chrCode.csv
    wc -l af-csq-go-filtered-VEP_$chrCode.csv
}

VEP
annotate
filter
```


## Code

The main script for the annotation is ```run_VEP.sh``` which expects as input a single allelic VCF file (only one alternative allele allowed per line) and outputs a VCF file containing the annotations as well as the VEP statistics in HTML format.
The output VCF file is then parsed by ```processVEP.sh``` which formats the annotation in to separate columns and recodes the genotypes as 0 (WT), 1 (HET), 2 (HOM) or NA for missing.  This script creates three separate files from one input file:
*  annotation file
*  genotype file
*  genotype quality file.

Finally, the output of ```processVEP.sh``` is processed by ```statsVEP.sh``` which does some sanity checks, generates summary stats and plots.

There are two subdirectories:

* annotate
* filters

## Additional annotation

VEP has builtin annotation but can also take custom annotation.
There is also a number of [plugins available](https://github.com/ensembl-variation/VEP_plugins)

## Consequence annotation scores

Variants in coding regions can have a functional impact.
Beside the built-in Consequence field which is output by VEP there are additional plugins to score the likely function impact of variants:
* [CADD](http://cadd.gs.washington.edu/): Combined Annotation Dependent Depletion<sup>[1][CADD]</sup> 
* [CONDEL](http://bg.upf.edu/fannsdb/): CONsensus DELeteriousness score that combines various tools (MutationAssessor, FATHMM).
* [POLYPHEN](http://genetics.bwh.harvard.edu/pph2/): Polymorphism Phenotyping v2) is a tool which predicts possible impact of an amino acid substitution on the structure and function of a human protein using straightforward physical and comparative considerations.
* [SIFT](http://sift.jcvi.org/): predicts whether an amino acid substitution affects protein function
* [CAROL](https://www.sanger.ac.uk/resources/software/carol/): Combined Annotation scoRing toOL (CAROL) combines information PolyPhen-2 and SIFT

[CADD]: http://www.nature.com/ng/journal/v46/n3/full/ng.2892.html  "Kircher M, Witten DM, Jain P, O'Roak BJ, Cooper GM, Shendure J. A general framework for estimating the relative pathogenicity of human genetic variants. Nat Genet. 2014 Feb 2."

## Frequency annotation

* [1000 genomes](http://www.1000genomes.org/)
* [ExAC](http://exac.broadinstitute.org/)
* UCLex

### Filters

The first step of filtering is on the allele frequency.
Filtering can also be done on consequence.
Finally filtering is usualldy done on expression data.

### Consequence filtering

Option 1:
```
[frameshift, stop, start, splice] OR [missense AND [Condel deleterious OR CAROL deleterious OR CADD>20]]
```

Option 2:
```
CADD>x (x=10?) AND [frameshift, stop, start, splice, deleterious missense in CSQ]
```



