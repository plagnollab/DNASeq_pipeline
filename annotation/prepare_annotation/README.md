
# Preparing annotation for VEP

We can add custom annotation to VEP such as allele frequencies, gene expression and pretty much anything else.
The standard format for adding custom annotation is:
```
CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
22	16157496	A>G:NA	A	G	36.1	.	.
22	16157517	G>T:0.000000	G	T	311.39	.	.
22	16157603	G>C:0.077114	G	C	1850320	.	.
22	16157827	C>G:0.251214	C	G	469005	.	.
22	16157883	T>C:0.215673	T	C	201443	.	.
22	16157912	C>T:0.075910	C	T	115155	.	.
22	16157913	T>G,A:0.247173	T	G,A	275676	.	.
22	16157940	C>T:0.210526	C	T	2408.86	.	.
22	16157982	G>T:0.166667	G	T	422.17	.	.
```
# Converting between build coordinates

Often files may not be for the right build, we can convert between assemblies using the liftOver tool from UCSC.
