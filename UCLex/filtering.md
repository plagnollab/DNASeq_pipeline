
# Filtered Annotation Files

```process_multiVCF.R``` generates filtered variants datasets from the UCLex calls and the variants annotated by ANNOVAR:

* all_variants/
* chrX
* comp_het_candidates.csv
* compound_hets/
* hom_candidates.csv
* hom_mapping/
* homozygous_variants/
* known_genes/
* known_genes.csv 
* rare_variants/
* signatures.csv
* support
* time_stamp.txt

# Annotation Headers

The headers for these files are explained below:

# Basic
```
Chr: Chromosome of variant
Start: Start coordinate of variant
End: End coordinate of variant
Ref: Reference allele for variant
Obs: Observed alllele for variant
Otherinfo: Other information
```

```
Samples
HUGO
Description
IRDC_batch1_LDS_4001_001_543
is.indel
```

## GATK variant quality 
```
FILTER: VQSR filter
QUAL
```

## Gene

```
Func: Function of the variant -- exonic, intronic, UTR, etc.
Gene: Gene name for variant
ensemblID
ensemblID.bis
ExonicFunc: If the variant is exonic, synonymous, non-synonymous, indel, etc.
AAChange: If exonic, variant change in nucleotide and protein format
Conserved
SegDup: Indicates if the variant is located in a segmental duplication region
```

## Frequency

```
ESP5400_ALL: MAF in Exome Sequencing Project dataset (5,400 exomes) for all populations
1000g2012feb_ALL: MAF in 1000Genomes February 2012 release
cg69: allele frequency in 69 human subjects sequenced by Complete Genomics
dbSNP137: RS# from the dbSNP database
```

```
AVSIFT: whole-exome SIFT scores for non-synonymous variants (obselete and should not be uesd any more)
SIFT Pathogenicity score: closer to 0 is more damaging
```

## LJB* (dbNSFP) non-synonymous variants annotation

The [LJB* databases](https://sites.google.com/site/jpopgen/dbNSFP) (for historical reasons, it is named as ljb rather than dbNSFP in ANNOVAR) include SIFT scores, PolyPhen2 HDIV scores, PolyPhen2 HVAR scores, LRT scores, MutationTaster scores, MutationAssessor score, FATHMM scores, GERP++ scores, PhyloP scores and SiPhy scores. As of October 2015, the lastest ljb database is dbnsfp30a. Previously, this dataset is referred to as ljb26, ljb23, ljb2, ljb, which caused confusions among many users. Starting from version 3.0a, we will adopt the dbnsfp keyword for future updates. Additionally, users should start to use table_annovar.pl to calculate scores for non-synonymous variants, since we no longer provide individual scores for individual algorithms.

```
LJB_PhyloP: Pathogenicity score from dbNSFP: conserved > 0.95, not conserved < 0.95
LJB_PhyloP_Pred: Pathogenicity call from dbNSFP: C - conserved, N - not conserved
LJB_SIFT: Pathogenicity score from dbNSFP: tolerated < 0.95, deleterious > 0.95
LJB_SIFT_Pred: Pathogenicity call from dbNSFP: T - tolerated, D - deleterious
LJB_PolyPhen2: Pathogenicity score from dbNSFP: probably damaging > 0.85, possibly damaging 0.85-0.15, benign < 0.15
LJB_PolyPhen2_Pred: Pathogenicity call: D - probably damaging, P - possibly damaging, B - benign
LJB_LRT: Pathogenicity probability score from dbNSFP: closer to 1 is more likely to be damaging -- see below
LJB_LRT_Pred: Pathogenicity score from dbNSFP: D - deleterious fulfills the following: (i) from a codon defined by LRT as significantly constrained (LRTorio0.001 and oo1), (ii) from a site with Z10 eutherian mammals alignments, and (iii) the alternative AA is not presented in any of the eutherian mammals
N - otherwise neutral
LJB_MutationTaster: Pathogenicity probability score from dbNSFP: closer to 1 is more likely to be damaging -- see below
LJB_MutationTaster_Pred: Pathogenicity score from dbNSFP: automatically calculated categories: ‘‘disease_causing_automatic,’’ ‘‘disease_causing,’’ ‘‘polymorphism,’’ and ‘‘polymorphism_automatic,’’ which we coded as ‘‘A,’’ ‘‘D,’’ ‘‘N,’’ and ‘‘P,’’
LJB_GERP++: Nucleotide conservation score from dbNSFP GERP: Higher number is more conserved, > 0 is generally conserved
```

```
Omim: Online Mendelian Inheritance in Man
```

```
signature
clean.signature
indel.length
dup.region
somewhat.rare
rare
novel
```

## Splicing info
```
exonic.splicing
splicing
core.splicing
HUGO.no.splice.info
```

```
lof: Loss of Function
non.syn
remove.bad.transcripts
```


## Case/control information from UCLex

```
freq.cases
non.missing.cases
non.ref.calls.cases
ncarriers.cases
missing.rate.cases
```

```
freq.controls
non.missing.controls
non.ref.calls.controls
```

```
non.missing.external.controls
freq.external.controls
non.ref.calls.external.controls
```

````
potential.comp.het
```


