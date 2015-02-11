# New

The main script for the annotation is ```run_VEP.sh``` which expects as input a single allelic VCF file (only one alternative allele allowed per line) and outputs a VCF file containing the annotations as well as the VEP statistics in HTML format.
The output VCF file is then parsed by ```processVEP.sh``` which formats the annotation in to separate columns and recodes the genotypes as 0 (WT), 1 (HET), 2 (HOM) or NA for missing.  This script creates three separate files from one input file: an annotation file, a genotype file and genotype quality file.
Finally, the output of ```processVEP.sh``` is processed by ```statsVEP.sh``` which does some sanity checks, generates summary stats and plots.



# Old


See http://www.ensembl.org/info/docs/tools/vep/script/vep_other.html#assembly

If using Build 37:

```
--ASSEMBLY GRCh37 \
--port 3337
```

If using Build 38:
```
--ASSEMBLY GRCh38
```

Script to convert multiallelic variants to multiple variants each with a single alternate allele
(note that quality score information is not retained). This reads in a gz vcf.
```
multi_to_single_gvcf.py
 ```
Script to run VEP:
```
bash run_VEP.sh
```
You will see that this script points to VEP in my home directory which is obviously not going to work.
In addition I was using version 75. They are now up to version 78 but we need to check that you can still use this with Build 37. You will see that I ask it for the effect of the variant on every possible Ensembl transcript.
One can just ask for the most damaging or the canonical.
We need to review the information here http://www.ensembl.org/info/docs/tools/vep/script/vep_options.html.
 
After running the VEP I use this script to format the results.

```
extract_VEP.py
```

Like this:
```sh
python extract_VEP.py groups.txt VEP_out.vcf > Formatted.txt
```

This also extracts the effects for the correct allele (e.g. if there are C>T and C>A allele frequency annotations,
the VEP will include both but we are only interested in the one that is relevant given our particular alternate allele).
 
The groups file looks like this:
```
/cluster/project8/vyp/AdamLevine/UCL-exomes_v2/groups.txt
```

This calculates the number of WT, HET, HOM, etc. individuals for each variant within a group
(e.g. all samples or those with a particular phenotype, etc.). I find this useful but it is not necessary.
 
The VEP_out.vcf is the output from the VEP.
You can see examples here:
```
/cluster/project8/vyp/AdamLevine/UCL-exomes_v2/annotate
```

Look at `VEP_1.vcf` and then `Formatted_1_with_genotypes.txt`
 
I am currently formatting the ExAC data so that we can include this.
 
Vincent, it may take us a little while to get this up and running
(we will need to install the VEP scripts in a shared directory somewhere).
In the meanwhile I am happy to run it myself on your data (if it is build 37).
You just need to point me to the VCFs.
 

 
