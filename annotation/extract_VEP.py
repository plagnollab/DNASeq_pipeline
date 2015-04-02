import sys

#file which specifies which group a sample belongs to
groups=dict([])
for i in open(sys.argv[1],"r"):
	s=i.strip().split()
	group=s[0]
	sample=s[1]
	retrieved = groups.get(group,"empty")
	if retrieved != "empty": groups[group] = groups[group]+","+sample
	if retrieved == "empty": groups[group] = sample

def divide(numerator, denominator):
	if denominator==0:
		return "NA"
	else:
		return str(round(float(numerator)/denominator,3))

count_labels = ["WT/HET/HOM/MISS","AF","MF","MISS-F","ALLELE-AF","ALLELE-MISS","ALLELE-ALT","ALLELE_TOTAL","ALLELE_MF"]

first=0
extra_fields=[]
#genotype columns (per person)
genotype=dict([]) 
#group columns
gcolumns=dict([]) 
for i in open(sys.argv[2],"r"):
	s=i.strip().split("\t")
	if s[0][0:6]=="##INFO":
		this_info = s[0].split(",")[0].split("##INFO=<ID=")[1]
		if this_info=="CSQ":
			info_header = s[0].split("Format: ")[1].split('">')[0].split("|")
		else:
			extra_fields.append(this_info)
         continue
	if s[0]=="#CHROM":
		genotype_header = s[9:]
		beginning_header = s[0:5]
		count=0
		counts_header = []
		for j in s[9:]:
			genotype[j]=count
			count=count+1
		for j in groups.keys():
			this = []
			for k in groups[j].split(","): this.append(genotype[k])
			gcolumns[j]=this
			for k in count_labels: counts_header.append(j+"_"+k)
         continue
	if s[0][0]=="#": continue
    if first==0:
        print "\t".join(beginning_header)+"\t"+"\t".join(info_header)+"\t"+"\t".join(extra_fields)+"\t"+"\t".join(counts_header)+"\t"+"\t".join(genotype_header)
        #add genotype data
        first=1
    #Genotype counts
    all_genotypes = s[9:]
    counts = []
    for j in groups.keys():
        this_group_genotypes=[]
        for k in gcolumns[j]:
            this_group_genotypes.append(all_genotypes[k])
        wt = str(this_group_genotypes.count("0/0"))
        het = str(this_group_genotypes.count("0/1"))
        hom = str(this_group_genotypes.count("1/1"))
        miss = str(this_group_genotypes.count("./."))
        af=divide(float(int(het)+2*int(hom)),float(2*int(wt)+2*int(het)+2*int(hom)))
        mutant=divide(float(int(het)+int(hom)),float(int(wt)+int(het)+int(hom)))
        missf=divide(int(miss),int(wt)+int(het)+int(hom)+int(miss))
        #extra to enable incoporation of imputed data with only one allele available
        by_allele=[]
        allele_total=0
        allele_alt=0
        for k in this_group_genotypes:
            k_split = k.split("/")
            by_allele.append(k_split[0])
            by_allele.append(k_split[1])
            if by_allele!="./.":
                allele_total=allele_total+1
                found_alt=k_split.count("1")
                if found_alt > 0:
                    allele_alt=allele_alt+1
        allele_total=str(allele_total)
        allele_alt=str(allele_alt)
        ref_count = str(by_allele.count("0"))
        alt_count = str(by_allele.count("1"))
        miss_count = str(by_allele.count("."))
        allele_af = divide(float(int(alt_count)),float(int(ref_count)+int(alt_count)))
        allele_mf = divide(float(int(allele_alt)),float(int(allele_total)))
        allele_miss = divide(float(int(miss_count)),float(int(ref_count)+int(alt_count)+int(miss_count)))
        counts.extend([wt+"/"+het+"/"+hom+"/"+miss,af,mutant,missf,allele_af,allele_miss,allele_alt,allele_total,allele_mf])
    #Annotations
    ref = s[3]
    alt = s[4]
    csq_out=[]
    info = s[7].split(";")
    csq=[]
    csq_out=[]
    #position of 'extra' information in info (alters if CSQ present)
    extra_n=0 
    if info[0][0:3]=="CSQ":
        extra_n=1
        csq = info[0].split("CSQ=")[1].split(",")
        for j in csq:
            this = j.split("|")
            out=[]
            for k in this:
                if k=="":
                    out.append("-")
                if k!="":
                    out.append(k)
            csq_out.append(out)
    beginning = s[0:5]
    extra=dict([])
    if len(alt.split(","))==1:
        for j in info[extra_n:]:
            j_s = j.split("=")
            key= j_s[0]
            original_value=j_s[1]
            new_value=j_s[1]
            check_allele_keys=["CADD","EAf","AAf","VPf","ONEKGf","CEUBGRf","okg","okgAMR","okgAFR","okgASN","okgEUR"]
            if key in check_allele_keys:
                new_value="-"
                for k in original_value.split(","):
                    k_s = k.split(":")
                    info_alleles = k_s[0]
                    info_ref = info_alleles.split(">")[0]
                    info_alt = info_alleles.split(">")[1]
                    if info_ref==ref and info_alt==alt:
                        new_value = k_s[1]
            extra[key]=new_value
    extra_out = []
    for j in extra_fields:
        extra_out.append(extra.get(j,"-"))
    if csq_out==[]:
        print "\t".join(beginning) + "\t" + "-\t"*len(info_header) +"\t".join(extra_out)+"\t"+"\t".join(counts)+"\t"+"\t".join(all_genotypes)
    if csq_out!=[]:
        for j in csq_out:
            print "\t".join(beginning) + "\t" + "\t".join(j)+"\t"+"\t".join(extra_out)+"\t"+"\t".join(counts)+"\t"+"\t".join(all_genotypes)

