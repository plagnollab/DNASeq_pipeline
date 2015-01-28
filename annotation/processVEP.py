#! /bin/env python
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

OUTPUT=OUTPUT+CSQ+EXAC+ONEKG

print '##fileformat=VCFv4.1'
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
        SAMPLES=HEADER[9:]
        print '#'+'\t'.join(OUTPUT+SAMPLES)
        continue
    s=l.strip().split('\t')
    # s will contain all fields for this line (including CSQ)
    s=dict(zip(HEADER, s))
    for k in SAMPLES: s[k]=s[k].split(':')[0]
    alternative_alleles=s['ALT'].split(',')
    for alt in alternative_alleles:
        s['ALT']=alt
        INFO=dict([tuple(x.split('=')) for x in s['INFO'].split(';')])
        # AFs
        # only include AFs where the ref>alt match
        s.update(dict([(k,af.split(':')[1],) for k in set(INFO).intersection(EXAC) for af in INFO[k].split(',') if ('>'.join([s['REF'],s['ALT']])==af.split(':')[0])]))
        s.update(dict([(k,af.split(':')[1],) for k in set(INFO).intersection(ONEKG) for af in INFO[k].split(',') if ('>'.join([s['REF'],s['ALT']])==af.split(':')[0])]))
        if 'CSQ' in INFO:
            cons=[ dict(zip(CSQ_HEADER,[b for b in a.split('|')])) for a in INFO['CSQ'].split(',') if a.split('|')[0]==s['ALT']]
            for csq in CSQ_HEADER: s[csq] = ','.join([co[csq] for co in cons])
        #print output, anything which was not found in the line gets a '.'
        print '\t'.join([s.get(h,'.')  for h in OUTPUT+SAMPLES])
        #','.join([csq['Feature']  for csq in cons if csq['Feature_type']=='Transcript']), [s[k] for k in s if 'EXAC' in k]




