The following analyses can be conducted on trios.

```
fasta=/scratch2/vyp-scratch2/reference_datasets/human_reference_sequence/human_g1k_v37.fasta
#BASEDIR=/goon2/scratch2/vyp-scratch2/RoseRichardson/1327-UCL-Fam2_Report/
#ped=$BASEDIR/pedigree_details.ped
GATK="/share/apps/jdk1.7.0_45/jre/bin/java -Djava.io.tmpdir=/scratch0/ -Xmx4g -Xms4g -jar /cluster/project8/vyp/vincent/Software/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar -R $fasta -T"
```

# Phasing

```
ls -1 align/data/*sorted_unique.bam > bam.list
output=gvcf/data/chr${f}.vcf.gz

$GATK HaplotypeCaller -I bam.list -nct 12 --variant_index_type LINEAR --variant_index_parameter 128000 -stand_call_conf 30.0 -stand_emit_conf 10.0 -L ${chrCode} --downsample_to_coverage 200 --GVCFGQBands 10 --GVCFGQBands 20 --GVCFGQBands 50 -o $output
```

# Autosomal dominant
Since both parents are unaffected, do these represent de novo mutations not seen in the parents?

# Autosomal recessive
So this would be a rare variant, het in the parents but hom in the affected child.

# Compound hets analysis:
There are two interpretations:

1. The first requires phasing.
A parent is het at different locations in the gene but due to recombination in the parent both damaging variants end up on the same chromosome.
Hence child carries both damaging variants on the same chromosome, gene doesn't work.
2. The other interpretation is that this is a burden were the parents are het at different variants, and the child has a higher burden of hets in a gene.
Both versions of gene are defective.  Burden in child higher than in parents.

# LOH analysis
a deletion has occurred so that neighbouring markers all appear hom or wt while they are het in parents

# CNV analysis
works on bam file as is done by ExomeDepth:
