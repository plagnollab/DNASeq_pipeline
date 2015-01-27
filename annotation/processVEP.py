#! /bin/env python
import sys
import re

# check the INFO header of your VCF for more ouput options
# also check the INFO CSQ field
OUTPUT= ['CHROM', 'POS', 'ID', 'Gene', 'REF', 'ALT', 'Allele', 'Feature', 'Consequence', 'Amino_acids', 'CANONICAL']
# allele freq (AF) in various populations
EXAC=['EXAC_AFR', 'EXAC_AMR', 'EXAC_Adj', 'EXAC_EAS', 'EXAC_FIN', 'EXAC_NFE', 'EXAC_OTH', 'EXAC_SAS']
ONEKG=['1KG_EUR', '1KG_AFR', '1KG_AMR', '1KG_ASN', 'ESP_EA', 'ESP_AA']
UCL=['UCLEX']

OUTPUT=OUTPUT+EXAC+ONEKG

print '##fileformat=VCFv4.1'
print '#'+'\t'.join(OUTPUT)
for l in sys.stdin:
    #get the format of the VEP consequence (CSQ) field
    #the names of the VEP CSQ fields are "|" delimited
    if l.startswith('##INFO=<ID=CSQ,'):
        CSQ_HEADER=re.compile('Format: (.*)">').search(l).group(1).split('|')
        continue
    #ignore all other '##' lines
    if l.startswith('##'): continue
    #line which starts with a single '#' is the header
    if l.startswith('#'):
        l=l.strip('#')
        HEADER=l.strip().split('\t')
        continue
    s=l.strip().split('\t')
    s=dict(zip(HEADER, s))
    INFO=dict([tuple(x.split('=')) for x in s['INFO'].split(';')])
    # AFs
    # only include AFs where the ref>alt match
    s.update(dict([(k,af.split(':')[1],) for k in set(INFO).intersection(EXAC) for af in INFO[k].split(',') if ('>'.join([s['REF'],s['ALT']])==af.split(':')[0])]))
    s.update(dict([(k,af.split(':')[1],) for k in set(INFO).intersection(ONEKG) for af in INFO[k].split(',') if ('>'.join([s['REF'],s['ALT']])==af.split(':')[0])]))
    if 'CSQ' in INFO:
        cons=[ dict(zip(CSQ_HEADER,[b for b in a.split('|')])) for a in INFO['CSQ'].split(',') ]
        for csq in CSQ_HEADER: s[csq] = ','.join([co[csq] for co in cons])
    #print output, anything which was not found in the line gets an empty string
    print '\t'.join([s.get(h,'')  for h in OUTPUT])
    #','.join([csq['Feature']  for csq in cons if csq['Feature_type']=='Transcript']), [s[k] for k in s if 'EXAC' in k]




