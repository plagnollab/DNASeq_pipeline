# parse vcf file and combine with annovar table
import sys
import os
import csv
import re

vcf_file = sys.argv[1]
annovar_file = sys.argv[2]
outfile = sys.argv[3]

# remove all info except for genotype
def clean_sample_for_geno(ll_moop, AN):
        new_ll = []
        for geno_info in ll_moop:
            depth = geno_info.split(":")[2]
            geno = geno_info.split(":")[0].replace("/","|")
            if (depth == ""):
                depth = 0
            if geno.startswith("."):
                new_ll.append(".|.:%s"%(depth))
                #new_ll.append("%s:%s"%(geno_info.split(":")[0].replace("/","|"),"0"))
            else:
                alleles = geno.split("|")
                if (float(alleles[0]) == AN): alleles[0] = "1"
                else: alleles[0] = "0"
                if (float(alleles[1]) == AN): alleles[1] = "1"
                else: alleles[1] = "0"
                if ( (alleles[0] == "1") & (alleles[1] == "0")):
                    alleles[0] = "0"
                    alleles[1] = "1"
                new_ll.append("%s|%s:%s"%(alleles[0], alleles[1], depth))
        return new_ll

annovar_dict = {}
annovar_multi = {}

with open(annovar_file) as f:
    read_csv_anno = csv.reader(f, delimiter=',',quotechar='"')
    ############# read the annovar file
    header = []
    for ll in read_csv_anno:
        if ll[0].find("Func") > -1:
            header = ll[ 0:28 ] 
            continue
        if len(ll) < 8 : continue
        [chr_,start,ref,obs] = [ll[28],ll[29],ll[31],ll[32]]
        key = "%s_%s_%s_%s"%(chr_, start,ref,obs);
        ######### Now deal with multi-allelic variants by modifying the key to account for the allele code
        if ll[32].find(",") > -1:
            if annovar_multi.has_key(key):
                annovar_multi[ key ] = annovar_multi[ key ] + 1
            else:
                annovar_multi[ key ] = 1
            key = "%s_%s"%(key,annovar_multi[ key ] )
        annovar_dict[ key ] = ll[ 0:28 ]   ##this key includes the _AN information for multiallelic positions
    f.close()

# get names
with file(vcf_file,'r') as read_vcf:
    for ll in read_vcf:
        if ll.startswith('#CHROM'):
            names=ll.strip().split('\t')[9:]
            break
    read_vcf.close()

with open(outfile,"w") as f:
    write_csv = csv.writer(f, delimiter=',',quotechar='"', quoting=csv.QUOTE_ALL)
    write_csv.writerow(header+['QUAL', 'FILTER'] + names)
    read_vcf = file(vcf_file,'rb')
    ##now read the VCF file and find the matching data in the annovar one
    for ll in read_vcf:
        if ll.startswith("#"): continue
        ll=ll.strip().split('\t')
        [chr_,start,ref,obs] = [ll[0],ll[1],ll[3],ll[4]]
        key = "%s_%s_%s_%s"%(chr_, start,ref,obs);
        if ',' in obs:  ##tricky bit, this is now multiallelic
            splobs = obs.split(",")
            AN = 1
            for alt in splobs:
                mkey = "%s_%s"%(key,AN)  ##the modified key that includes the AN
                llnew = ll[5:7] + clean_sample_for_geno(ll[9:], AN)
                AN = AN + 1
                if mkey in annovar_dict:
                    write_csv.writerow(annovar_dict[ mkey ]+llnew)
                else:
                    print "%s not found in annovar dictionary"%(mkey)
        ##the easy case where the is a single allele:
        else:
            llnew = ll[5:7] + clean_sample_for_geno(ll[9:], 1)
            #llnew.replace("1|0", "0|1")
            if key in annovar_dict:
                write_csv.writerow(annovar_dict[ key ]+llnew)
            else:
                print "%s not found in annovar dictionary"%(key)
    f.close() 
    read_vcf.close()

    
