#! /bin/env python

from __future__ import print_function
import argparse
import os.path
import sys
import tabix
import csv

usage_example = """
"""

parser = argparse.ArgumentParser(formatter_class = argparse.RawDescriptionHelpFormatter, epilog = usage_example) 
#compulsory arguments
parser.add_argument('--file', dest='file', help = "list of variants which we are interested in", required=False, default=None)
args = parser.parse_args()

f=args.file

variants=['%s:%s-%s' % (l['Chr'],l['Start'],l['End'],) for l in csv.DictReader(file(f,'r'))]


##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as
##INFO=<ID=AC_AFR,Number=A,Type=Integer,Description="African/African American Allele Counts">
##INFO=<ID=AC_AMR,Number=A,Type=Integer,Description="American Allele Counts">
##INFO=<ID=AC_Adj,Number=A,Type=Integer,Description="Adjusted Allele Counts">
##INFO=<ID=AC_EAS,Number=A,Type=Integer,Description="East Asian Allele Counts">
##INFO=<ID=AC_FIN,Number=A,Type=Integer,Description="Finnish Allele Counts">
##INFO=<ID=AC_Hemi,Number=A,Type=Integer,Description="Adjusted Hemizygous Counts">
##INFO=<ID=AC_Het,Number=A,Type=Integer,Description="Adjusted Heterozygous Counts">
##INFO=<ID=AC_Hom,Number=A,Type=Integer,Description="Adjusted Homozygous Counts">
##INFO=<ID=AC_NFE,Number=A,Type=Integer,Description="Non-Finnish European Allele Counts">
##INFO=<ID=AC_OTH,Number=A,Type=Integer,Description="Other Allele Counts">
##INFO=<ID=AC_SAS,Number=A,Type=Integer,Description="South Asian Allele Counts">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=AN_AFR,Number=1,Type=Integer,Description="African/African American Chromosome Count">
##INFO=<ID=AN_AMR,Number=1,Type=Integer,Description="American Chromosome Count">
##INFO=<ID=AN_Adj,Number=1,Type=Integer,Description="Adjusted Chromosome Count">
##INFO=<ID=AN_EAS,Number=1,Type=Integer,Description="East Asian Chromosome Count">
##INFO=<ID=AN_FIN,Number=1,Type=Integer,Description="Finnish Chromosome Count">
##INFO=<ID=AN_NFE,Number=1,Type=Integer,Description="Non-Finnish European Chromosome Count">
##INFO=<ID=AN_OTH,Number=1,Type=Integer,Description="Other Chromosome Count">
##INFO=<ID=AN_SAS,Number=1,Type=Integer,Description="South Asian Chromosome Count">
##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt Vs. Ref base quali
##INFO=<ID=CCC,Number=1,Type=Integer,Description="Number of called chromosomes">
##INFO=<ID=ClippingRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref number
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP Membership">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
##INFO=<ID=DS,Number=0,Type=Flag,Description="Were any of the samples downsampled?">
##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">
##INFO=<ID=FS,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher's exact test to detect strand bias">
##INFO=<ID=GQ_MEAN,Number=1,Type=Float,Description="Mean of all GQ values">
##INFO=<ID=GQ_STDDEV,Number=1,Type=Float,Description="Standard deviation of all GQ values">
##INFO=<ID=HWP,Number=1,Type=Float,Description="P value from test of Hardy Weinberg Equilibrium">
##INFO=<ID=HaplotypeScore,Number=1,Type=Float,Description="Consistency of the site with at most two segregating haplot
##INFO=<ID=Hemi_AFR,Number=A,Type=Integer,Description="African/African American Hemizygous Counts">
##INFO=<ID=Hemi_AMR,Number=A,Type=Integer,Description="American Hemizygous Counts">
##INFO=<ID=Hemi_EAS,Number=A,Type=Integer,Description="East Asian Hemizygous Counts">
##INFO=<ID=Hemi_FIN,Number=A,Type=Integer,Description="Finnish Hemizygous Counts">
##INFO=<ID=Hemi_NFE,Number=A,Type=Integer,Description="Non-Finnish European Hemizygous Counts">
##INFO=<ID=Hemi_OTH,Number=A,Type=Integer,Description="Other Hemizygous Counts">
##INFO=<ID=Hemi_SAS,Number=A,Type=Integer,Description="South Asian Hemizygous Counts">
##INFO=<ID=Het_AFR,Number=A,Type=Integer,Description="African/African American Heterozygous Counts">
##INFO=<ID=Het_AMR,Number=A,Type=Integer,Description="American Heterozygous Counts">
##INFO=<ID=Het_EAS,Number=A,Type=Integer,Description="East Asian Heterozygous Counts">
##INFO=<ID=Het_FIN,Number=A,Type=Integer,Description="Finnish Heterozygous Counts">
##INFO=<ID=Het_NFE,Number=A,Type=Integer,Description="Non-Finnish European Heterozygous Counts">
##INFO=<ID=Het_OTH,Number=A,Type=Integer,Description="Other Heterozygous Counts">
##INFO=<ID=Het_SAS,Number=A,Type=Integer,Description="South Asian Heterozygous Counts">
##INFO=<ID=Hom_AFR,Number=A,Type=Integer,Description="African/African American Homozygous Counts">
##INFO=<ID=Hom_AMR,Number=A,Type=Integer,Description="American Homozygous Counts">
##INFO=<ID=Hom_EAS,Number=A,Type=Integer,Description="East Asian Homozygous Counts">
##INFO=<ID=Hom_FIN,Number=A,Type=Integer,Description="Finnish Homozygous Counts">
##INFO=<ID=Hom_NFE,Number=A,Type=Integer,Description="Non-Finnish European Homozygous Counts">
##INFO=<ID=Hom_OTH,Number=A,Type=Integer,Description="Other Homozygous Counts">
##INFO=<ID=Hom_SAS,Number=A,Type=Integer,Description="South Asian Homozygous Counts">

tb=tabix.open('/cluster/project8/IBDAJE/VEP_custom_annotations/GRCh37/ExAC/0.3/ExAC.r0.3.sites.vep.vcf.gz')


col1=['CHROM','POS','ID','REF','ALT','QUAL','FILTER']
col2=['AC', 'AC_AFR', 'AC_AMR', 'AC_Adj', 'AC_EAS', 'AC_FIN', 'AC_Het', 'AC_Hom', 'AC_NFE', 'AC_OTH', 'AC_SAS', 'AF', 'AN', 'AN_AFR', 'AN_AMR', 'AN_Adj', 'AN_EAS', 'AN_FIN', 'AN_NFE', 'AN_OTH', 'AN_SAS', 'DP', 'FS', 'GQ_MEAN', 'GQ_STDDEV', 'Het_AFR', 'Het_AMR', 'Het_EAS', 'Het_FIN', 'Het_NFE', 'Het_OTH', 'Het_SAS', 'Hom_AFR', 'Hom_AMR', 'Hom_EAS', 'Hom_FIN', 'Hom_NFE', 'Hom_OTH', 'Hom_SAS', 'InbreedingCoeff', 'VQSLOD', 'culprit']
print(','.join(col1+col2))

for v in variants:
    records=tb.querys(v)
    for r in records:
        c1=','.join(r[0:7])
        #print( [ '%s,%s' % (), ','.join() for r in records]
        d=dict( [x for x in [x.split('=') for x in r[7].split('|')[0].split(';')] if len(x)==2] )
        c2=','.join([d[k] for k in col2])
        print(c1,c2,sep=',')


