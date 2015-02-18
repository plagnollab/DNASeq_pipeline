# parse vcf file and combine with annovar table
import sys
import os
import csv
import re

vcf_file = sys.argv[1]
annovar_file = sys.argv[2]
#sample_file = sys.argv[3]
outfile = sys.argv[3]

re_uclg_id='UCLG.*?(\\d+)'
    
rg = re.compile(re_uclg_id,re.IGNORECASE|re.DOTALL)


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
			if (float(alleles[0]) == AN):
		       		alleles[0] = "1"
			else:
		       		alleles[0] = "0"

			if (float(alleles[1]) == AN):
			       	alleles[1] = "1"
			else:
				alleles[1] = "0"
			if ( (alleles[0] == "1") & (alleles[1] == "0")):
				alleles[0] = "0"
				alleles[1] = "1"
				
			new_ll.append("%s|%s:%s"%(alleles[0], alleles[1], depth))
        return new_ll

    
def create_sample_dict(sample_file):
# make sample_dict

        sample_dict= {}
        
        # use a regular expression to identify the sample id and create a dictionary containing the data`
        
        
        for line in open(sample_file):
            
                if line.startswith("Sample"):
                    
                        continue
                
                ll =line.strip().split("\t")
                
                if len(ll) < 2:
                        continue
                
                #print ">>",line
                #print ll
              
                m = rg.search(ll[1])
                
                
            
                sample_dict["UCLG"+m.group(1)] = ll[0]
            
        return sample_dict


#sample_dict = create_sample_dict(sample_file)


annovar_dict = {}
vcf_dict = {}
annovar_multi = {}

names = []
######## vcf header reader
read_csv = csv.reader(open(vcf_file), delimiter='\t',quotechar='"')

for ll in read_csv:
    
        #if ll[0].find("##contig=<ID=") > -1:
                
         #       annovar_dict[ll[0].split("ID=")[1].split(",")[0]] = {}
         #       vcf_dict[ll[0].split("ID=")[1].split(",")[0]] = {}
                
        if ll[0].find("CHROM") > -1:
                
                names = []
                
                for n in ll[9:]:
                        
                        names.append(n)
                
        if ll[0].find("#") == -1:
                
                break


print "The names found in the VCF file are: ", names


read_csv_anno = csv.reader(open(annovar_file), delimiter=',',quotechar='"')


############# read the annovar file
header = []
for ll in read_csv_anno:
    
        if ll[0].find("Func") > -1:
                
		header = ll[ 0:28 ] 
		continue
        
        if len(ll) < 8 :
                continue
        
        [chr_,start,ref,obs] = [ll[28],ll[29],ll[31],ll[32]]
	key = "%s_%s_%s_%s"%(chr_, start,ref,obs);

	######### Now deal with multi-allelic variants by modifying the key to account for the allele code
	if ll[32].find(",") > -1:
		if annovar_multi.has_key (key):
			annovar_multi[ key ] = annovar_multi[ key ] + 1
		else:
			annovar_multi[ key ] = 1
		key = "%s_%s"%(key,annovar_multi[ key ] )

		
		
	annovar_dict[ key ] = ll[ 0:28 ]   ##this key includes the _AN information for multiallelic positions
        
        
        
write_csv = csv.writer(open(outfile,"w"), delimiter=',',quotechar='"', quoting=csv.QUOTE_ALL)
read_csv = csv.reader(open(vcf_file), delimiter='\t',quotechar='"')

write_csv.writerow(header+['QUAL', 'FILTER'] + names)
#write_csv.writerow(header+ names)

for ll in read_csv:  ##now read the VCF file and find the matching data in the annovar one
        
        if ll[0].find("#") > -1:
                continue
    
        [chr_,start,ref,obs] = [ll[0],ll[1],ll[3],ll[4]]
	key = "%s_%s_%s_%s"%(chr_, start,ref,obs);


	if obs.find(",") > -1:  ##tricky bit, this is now multiallelic
		splobs = obs.split(",")
		AN = 1
		for alt in splobs:
			mkey = "%s_%s"%(key,AN)  ##the modified key that includes the AN
			llnew = ll[5:7] + clean_sample_for_geno(ll[9:], AN)
       			AN = AN + 1
			try:
				write_csv.writerow(annovar_dict[ mkey ]+llnew)
			except:
				print "%s not found in annovar dictionary"%(mkey)
	else:		##the easy case where the is a single allele	
		llnew = ll[5:7] + clean_sample_for_geno(ll[9:], 1)
		#llnew.replace("1|0", "0|1")
		try:
			write_csv.writerow(annovar_dict[ key ]+llnew)
		except:
			continue
    

    
