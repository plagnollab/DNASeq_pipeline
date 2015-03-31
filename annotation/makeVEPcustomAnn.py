#! /bin/env python
from __future__ import print_function
import sys

#input headers
HEADERS=['CHROM','POS','N_ALLELES','N_CHR']
#,'{ALLELE:FREQ}']

#output headers
OUTPUT_HEADERS=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO']

def print_line(s):
    print( *([s[h] for h in OUTPUT_HEADERS]), sep='\t' )

print('##fileformat=VCFv4.1')
print('\t'.join(OUTPUT_HEADERS))

for line in sys.stdin:
    #header line
    if line.startswith('CHROM'): continue
    line=line.strip()
    #this is tab separated line
    s=line.split("\t")
    d=dict(zip(HEADERS,s[0:len(HEADERS)]))
    an=int(d['N_ALLELES'])
    alleles=s[len(HEADERS):]
    afs=dict(map(lambda x: (str(x.split(':')[0]),float(x.split(':')[1]),) , alleles))
    d['REF']=alleles[0].split(':')[0]
    for a in alleles[1:]:
        d['ALT']=a.split(':')[0]
        d['ID']='{}>{}'.format(d['REF'],a)
        d['QUAL']=d['FILTER']=d['INFO']='.'
        print_line(d)

