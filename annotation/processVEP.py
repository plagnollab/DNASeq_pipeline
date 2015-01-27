#! /bin/env python
import sys
import re

OUTPUT= ['CHROM', 'POS', 'ID', 'Gene', 'REF', 'ALT', 'Allele', 'Feature', 'Consequence', 'Amino_acids', 'CANONICAL']

EXAC=['EXAC_AFR', 'EXAC_AMR', 'EXAC_Adj', 'EXAC_EAS', 'EXAC_FIN', 'EXAC_NFE', 'EXAC_OTH', 'EXAC_SAS']
ONEKG=['1KG_EUR', '1KG_AFR', '1KG_AMR', '1KG_ASN', 'ESP_EA', 'ESP_AA']
UCL=['UCLEX']

OUTPUT=OUTPUT+EXAC+ONEKG

print '##fileformat=VCFv4.1'
print '#'+'\t'.join(OUTPUT)
for l in sys.stdin:
    if l.startswith('##INFO=<ID=CSQ,'):
        CSQ_HEADER=re.compile('Format: (.*)">').search(l).group(1).split('|')
        continue
    if l.startswith('##'): continue
    if l.startswith('#'):
        l=l.strip('#')
        HEADER=l.strip().split('\t')
        continue
    s=l.strip().split('\t')
    s=dict(zip(HEADER, s))
    INFO=dict([tuple(x.split('=')) for x in s['INFO'].split(';')])
    s.update(dict([(k,af.split(':')[1],) for k in set(INFO).intersection(EXAC) for af in INFO[k].split(',') if ('>'.join([s['REF'],s['ALT']])==af.split(':')[0])]))
    s.update(dict([(k,af.split(':')[1],) for k in set(INFO).intersection(ONEKG) for af in INFO[k].split(',') if ('>'.join([s['REF'],s['ALT']])==af.split(':')[0])]))
    if 'CSQ' not in INFO:
        s.update(dict(zip(CSQ_HEADER,['' for b in CSQ_HEADER])))
    else:
        cons=[ dict(zip(CSQ_HEADER,[b for b in a.split('|')])) for a in INFO['CSQ'].split(',') ]
        for csq in CSQ_HEADER: s[csq] = ','.join([co[csq] for co in cons])
    print '\t'.join([s.get(h,'')  for h in OUTPUT])
    #','.join([csq['Feature']  for csq in cons if csq['Feature_type']=='Transcript']), [s[k] for k in s if 'EXAC' in k]

    #print [(h,s[h],) for h in HEADER]



