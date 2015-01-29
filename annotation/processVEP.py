#! /bin/env python
from __future__ import print_function
import sys
import re

# check the INFO header of your VCF for more ouput options
# also check the INFO CSQ field
OUTPUT= ['CHROM', 'POS', 'ID', 'REF', 'ALT']
CSQ=['Allele','Gene','Feature','Feature_type','Consequence','cDNA_position','CDS_position','Protein_position','Amino_acids','Codons','Existing_variation','DISTANCE','STRAND','SYMBOL','SYMBOL_SOURCE','HGNC_ID','CANONICAL','SIFT','PolyPhen','GMAF','AFR_MAF','AMR_MAF','ASN_MAF','EUR_MAF','AA_MAF','EA_MAF','CLIN_SIG','SOMATIC','CAROL','Condel']
# allele freq (AF) in various populations
EXAC=['EXAC_AFR', 'EXAC_AMR', 'EXAC_Adj', 'EXAC_EAS', 'EXAC_FIN', 'EXAC_NFE', 'EXAC_OTH', 'EXAC_SAS']
ONEKG=['1KG_EUR', '1KG_AFR', '1KG_AMR', '1KG_ASN', 'ESP_EA', 'ESP_AA']
UCL=['UCLEX']

filename=sys.argv[1]
infile=open(filename,'r')
basename=filename.split('.')[0]

ANNOTATION_HEADER=['VID']+CSQ+EXAC+UCL+['AF','WT','HET','HOM','MISS']
annotation_file=open('-'.join([basename,'annotations.tab']), 'w+')
genotype_file=open('-'.join([basename,'genotypes.csv']), 'w+')
quality_file=open('-'.join([basename,'genotypes_depth.tab']), 'w+')
for l in infile:
    #get the format of the VEP consequence (CSQ) field
    #the names of the VEP CSQ fields are "|" delimited
    if l.startswith('##INFO=<ID=CSQ,'):
        CSQ_HEADER=re.compile('Format: (.*)">').search(l).group(1).split('|')
        continue
    #ignore all other '##' lines
    if l.startswith('##'): continue
    #line which starts with a single '#' is the header
    #When we get to the header line, print it out.
    if l.startswith('#'):
        l=l.strip('#')
        HEADER=l.strip().split('\t')
        SAMPLES=HEADER[9:]
        print(*(['VID']+SAMPLES),sep=',',file=genotype_file)
        print(*(['VID']+SAMPLES),sep='\t',file=quality_file)
        print(*(ANNOTATION_HEADER),sep='\t',file=annotation_file)
        continue
    s=l.strip().split('\t')
    # s will contain all fields for this line (including CSQ)
    s=dict(zip(HEADER, s))
    #for k in SAMPLES: s[k]=s[k].split(':')[0]
    for k in SAMPLES: s[k]=s[k]
    VARIANT_ID='_'.join([s['CHROM'],s['POS'],s['REF'],s['ALT']])
    s['VARIANT_ID']=VARIANT_ID
    #print output, anything which was not found in the line gets a '.'
    #','.join([csq['Feature']  for csq in cons if csq['Feature_type']=='Transcript']), [s[k] for k in s if 'EXAC' in k]
    #print '\t'.join( [VARIANT_ID] + [s.get(h,'.')  for h in OUTPUT] )
    # GENOTYPE
    # output of this goes to genotype.csv file
    # number of times we have the alternative allele
    def genotype(g):
        d=dict(zip(s['FORMAT'].split(':'),g.split(':')))
        geno=d['GT']
        if geno=='0/0' or geno=='./0' or geno=='0/.': 
            return 0
        elif geno=='0/1' or geno=='./1' or geno=='1/.': 
            return 1
        elif geno=='1/1': 
            return 2
        elif geno=='./.':
            return 'NA'
        else:
            print( VARIANT_ID, geno, sep=',', file=sys.stderr)
            print( l, file=sys.stderr)
            raise 'hell'
    GENOTYPES=[genotype(s.get(h,'.'))for h in SAMPLES]
    print(*([VARIANT_ID] + GENOTYPES),sep=',',file=genotype_file)
    # QUALITY
    def genotype_quality(g):
        d=dict(zip(s['FORMAT'].split(':'),g.split(':')))
        allele_depth=d['AD']
        total_depth=d['DP']
        return total_depth
    print(*([VARIANT_ID] + [genotype_quality(s.get(h,'.'))for h in SAMPLES]),sep=',',file=quality_file)
    # ANNOTATIONS
    INFO=dict([tuple(x.split('=')) for x in s['INFO'].split(';')])
    # AFs
    # only include AFs where the ref>alt match
    s.update(dict([(k,af.split(':')[1],) for k in set(INFO).intersection(EXAC) for af in INFO[k].split(',') if ('>'.join([s['REF'],s['ALT']])==af.split(':')[0])]))
    s.update(dict([(k,af.split(':')[1],) for k in set(INFO).intersection(ONEKG) for af in INFO[k].split(',') if ('>'.join([s['REF'],s['ALT']])==af.split(':')[0])]))
    # determine whether we have a deletion or an insertion
    # the CSQ Allele is set accordingly
    # deletion
    if len(s['REF']) > len(s['ALT']):
        s['Allele']='-'
    # insertion
    elif len(s['REF']) < len(s['ALT']):
        s['Allele']=s['ALT'].lstrip(s['REF'][:(len(s['REF'])-1)])
    # variant
    else:
        s['Allele']=s['ALT']
    if 'CSQ' in INFO:
        #read the comma separated list of CSQs
        cons=[ dict(zip(CSQ_HEADER,[b for b in a.split('|')])) for a in INFO['CSQ'].split(',') ]
        #only keep ones where the Allele matches the expected
        #print s['Allele']
        #print [ co['Allele'] for co in cons ]
        cons=[ co for co in cons if co['Allele']==s['Allele'] ]
        for csq in CSQ_HEADER: s[csq] = ','.join([co[csq] for co in cons])
    # calculate freq of variant in this batch
    s['MISS'] = float(GENOTYPES.count('NA')) / len(SAMPLES)
    s['HET'] = float(GENOTYPES.count(1)) / len(SAMPLES)
    s['HOM'] = float(GENOTYPES.count(2)) / len(SAMPLES)
    s['AF']=float(GENOTYPES.count(1)+GENOTYPES.count(2)*2) / 2*len(SAMPLES)
    print(*([VARIANT_ID] + [s.get(h,'.')for h in ANNOTATION_HEADER]),sep='\t',file=annotation_file)


