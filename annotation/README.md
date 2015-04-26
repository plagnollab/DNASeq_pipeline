# Variant annotation


## Summary
The annotation of both coding and non-coding variants is done using the [Variant Effect Predictor](http://www.ensembl.org/info/docs/tools/vep/index.html)(VEP).
All variants are annotated with [1000 Genomes](www.1000genomes.org) allele frequencies.
Additionally, coding variants are annotated with allele frequencies from [ExAC](http://exac.broadinstitute.org/).
Variants are further annotated with Combined Annotation Dependent Depletion (CADD)<sup>[1][CADD]</sup>, Combined Annotation scoRing toOL (CAROL) and CONsensus DELeteriousness (Condel) consequence scores, to assess their functional impact.

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
Filtering can also be done on consequence
Finally filtering is usualldy done on expression data



