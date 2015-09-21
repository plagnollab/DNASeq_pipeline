#! /bin/env python

from __future__ import print_function
import sys
import pysam

samfile = pysam.AlignmentFile("/goon2/scratch2/vyp-scratch2/WGS/Hardcastle_October2014/align/data/140724_H02J1_MH-901804_sorted_unique.bam",'rb')

# 'ChrX:139,400,000-139,600,000'
#for reads in samfile.pileup('chrX',139400000,139600000):

for read in samfile.fetch('chrX',139400000,139600000):
    if not read.mate_is_unmapped and samfile.getrname(mate.reference_id)!='chrX':
        mate=samfile.mate(read)
        print(read.mpos)
    else:
        print(read)


reads=[ r for r in samfile.fetch('chrX',139400000,139600000) ]

# percentage of unmapped reads
100 * float(len(filter(lambda x: x.is_unmapped, reads))) / len(reads)

# reads where the mate maps to another chromosome
for r in filter(lambda x: x.rname!=x.mrnm, reads):
    m=samfile.mate(r)
    print('%s:%d'%(samfile.getrname(r.rname),r.pos,), '%s:%d'%(samfile.getrname(m.rname), m.pos,))



