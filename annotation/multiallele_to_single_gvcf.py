import sys,gzip

vcf = gzip.open(sys.argv[1],"r")


print "##fileformat=VCFv4.1"
for i in vcf:
    s=i.rstrip("\n").split("\t")
    if s[0][0:2]!="##": #remove beginning
	if s[0][0]=="#":
		print "\t".join(s) #header
	if s[0][0]!="#": #data
		alt=s[4].split(",") #split alternate alleles
		s[5:9]=[".",".",".","GT"] #QUAL=., FILTER=., INFO=., FORMAT=GT
		if len(alt)==1: #if single alternate allele
			genotypes=[]
			for j in s[9:]:
				 genotypes.append(j.split(":")[0]) #extract genotype
			print "\t".join(s[0:9]) +"\t"+ "\t".join(genotypes) #print data
		if len(alt)>1:
			s[2] = "." #rs not necessarily correct
			possible=range(1,len(alt)+1)            
			possible_string=[]
			for each in possible:
				possible_string.append(str(each))
			alt_count=-1
			for j in possible_string:
				alt_count=alt_count+1
				genotypes=[]
				for k in s[9:]:
					genotypes.append(k.split(":")[0].split("/")[0])
					genotypes.append(k.split(":")[0].split("/")[1])
				if j in genotypes:
					new_genotypes=[]
					for l in genotypes:
						if l==j:
							new_genotypes.append("1")
						elif l=="." or l=="0":
							new_genotypes.append(l)
						else:
							new_genotypes.append("0")
					count=0
					out=[]
					for l in new_genotypes:
						count=count+1
						if count%2==0:
							if previous=='1' and l=='0':
								out.append('0/1')
							else:
								out.append(previous+"/"+l)
						previous=l
					print "\t".join(s[0:4])+"\t"+alt[alt_count]+"\t"+"\t".join(s[5:9])+"\t"+"\t".join(out)
