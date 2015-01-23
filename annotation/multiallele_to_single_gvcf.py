#! /bin/env python
import sys,gzip

vcf = gzip.open(sys.argv[1],"r")

#these 9 column headers are standard to all VCF files
STD_HEADERS=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']
#the samples headers depend on the number of samples in the file
#which we find out once we read the #CHROM line
SAMPLE_HEADERS=[]

for line in vcf:
    line=line.strip()
    #remove beginning
    #lines which start with '###'
    #are not tab separated
    if line.startswith('##'):
        #print(line)
        continue
    #this is tab separated line
    s=line.split("\t")
    #header line yay!: #CHROM ...
    if line.startswith('#'):
        print(line)
        headers=s
        headers[0]=headers[0].strip('#')
        #the first 9 names in the header are standard (see above)
        if (headers[0:len(STD_HEADERS)] != STD_HEADERS): raise 'hell'
        #everything else in the header is a sample name
        SAMPLE_HEADERS=headers[len(STD_HEADERS):]
        continue
    #you can now access elements by their header name
    s=dict(zip(headers,s))
    #QUAL=., FILTER=., INFO=., FORMAT=GT
    #just set everything except GT to .
    #@pontikos: doesn't VEP use the other fields though?
    for h in ['QUAL','FILTER','INFO']: s[h]='.'
    s['FORMAT']='GT'
    #split alternate alleles
    alternative_alleles=s['ALT'].split(",")
    format=s['FORMAT'].split(':')
    #extract the genotypes of the samples
    #genotype is always first element of array
    for h in SAMPLE_HEADERS: s[h]=dict(zip(format,s[h].split(":")))['GT']
    #if single alternate allele
    #then just print out as normal
    if len(alternative_alleles)<=1:
        print('\t'.join([s[h] for h in headers]))
        continue
    #otherwise split over as many lines as there are alternative alleles
    #each genotype then takes a 1 if matches the alternative or a 0 otherwise
    for i in range(len(alternative_alleles)):
        s['ALT']=alternative_alleles[i]
        print( '\t'.join( [s[h] for h in STD_HEADERS] + ['/'.join([['0','1'][int(int(g)==i)] for g in s[h].split('/')]) for h in SAMPLE_HEADERS] ) )



