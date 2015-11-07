
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

The headers for these files:
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
Func
ExonicFunc
HUGO
Description
non.ref.calls.cases
ncarriers.cases
missing.rate.cases
freq.controls
non.missing.controls
non.ref.calls.controls
non.missing.external.controls
freq.external.controls
AAChange
IRDC_batch1_LDS_4001_001_543
is.indel
QUAL
Gene
Conserved
```

```
ESP6500si_ALL: ESP
X1000g2012apr_ALL: 1000 genomes
dbSNP137: dbSNP
AVSIFT
```

```
LJB_PhyloP
LJB_PhyloP_Pred
LJB_SIFT
LJB_SIFT_Pred
LJB_PolyPhen2
LJB_PolyPhen2_Pred
LJB_LRT
LJB_LRT_Pred
LJB_MutationTaster
LJB_MutationTaster_Pred
LJB_GERP..
```
```
Func
Function of the variant -- exonic, intronic, UTR, etc.
Gene: Gene name for variant
ExonicFunc: If the variant is exonic, synonymous, non-synonymous, indel, etc.
AAChange: If exonic, variant change in nucleotide and protein format
Conserved
SegDup: Indicates if the variant is located in a segmental duplication region
ESP5400_ALL: MAF in Exome Sequencing Project dataset (5,400 exomes) for all populations
1000g2012feb_ALL: MAF in 1000Genomes February 2012 release
dbSNP135: RS# from the dbSNP database
AVSIFT
SIFT Pathogenicity score: closer to 0 is more damaging
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
cg69
Omim
```

```
FILTER: VQSR filter
signature
clean.signature
indel.length
ensemblID
ensemblID.bis
HUGO.no.splice.info
dup.region
somewhat.rare
rare
novel
exonic.splicing
splicing
core.splicing
lof
non.syn
remove.bad.transcripts
```

```
non.ref.calls.external.controls
freq.cases
non.missing.cases
potential.comp.het
```


