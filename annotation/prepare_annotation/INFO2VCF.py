#! /bin/env python
import sys

#these 4column headers are standard to all INFO files
STD_HEADERS=['CHROM','POS','REF','ALT']
#the samples headers depend on the number of samples in the file
#which we find out once we read the #CHROM line
SAMPLE_HEADERS=[]
#header line yay!
headers=sys.stdin.readline().strip().split('\t')
#the first 4 names in the header are standard (see above)
if (headers[0:len(STD_HEADERS)] != STD_HEADERS): raise 'hell'
#the remaining headers are the info header
INFO_HEADER=headers[len(STD_HEADERS):]

#these column headers are standard to all VCF files
VCF_HEADERS=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO']

print "##fileformat=VCFv4.1"
print '#'+'\t'.join(VCF_HEADERS)

for line in sys.stdin:
    line=line.strip()
    #this is tab separated line
    s=line.split("\t")
    #you can now access elements by their header name
    s=dict(zip(headers,s))
    #add these to output
    for h in ['QUAL','FILTER','INFO']: s[h]='.'
    #split alternate alleles extract the info field which refers to the alternative
    alternative_alleles=zip(s['ALT'].split(","),s[INFO_HEADER[0]].split(','))
    total_count=float(s[INFO_HEADER[1]])
    for alt,info, in alternative_alleles:
        if total_count !=0:
            freq=float(info)/total_count
        else:
            freq=0
        s['ID'] = s['REF']+'>'+alt+':'+"{:.8f}".format(freq)
        s['ALT'] = alt
        print( '\t'.join( [s[h] for h in VCF_HEADERS] ) )





