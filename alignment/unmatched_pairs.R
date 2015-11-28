library(Rsamtools)
 
#bamFile <- "C:/Users/Jessica Gardner/Desktop/LP0022129-DNA_A01.bam"
 
# reads in all files in the directory which finish in .bam
print(files <- grep('.bam$',list.files(path="Z:/WGS/data from Niko for Jess 25.08.15/JessicaGardner/",pattern=".*.bam",full.names = TRUE),value = T))
 for (f in files) {
  print(f)
  # skip file if size less than 100 (it is probably empty - only contains header information)
  if (file.info(f)$size<100) next
  # read file
  bams <- scanBam(f)
  # not sure
  b <- bam[[1]]
  # count to which chromosome the mates map
print(sort(table(b$mrnm),decreasing=T))
}
 
bamFile <- "Z:/WGS/bam files/X-LINKED RP/140714_H00TK_PF-19956/140714_H00TK_PF-19956_sorted_unique.bam"
bamFile <- "Z:/WGS/bam files/X-LINKED RP/140724_H02J1_MH-901804/140724_H02J1_MH-901804_sorted_unique.bam"
# "which" is the region you wish to pull out
# "what" are the paramaters which will be present in bam data structure
# a subset may be selected like so:
#what=c("rname", "strand", "pos", "qwidth", "seq", "mrnm", "mpos","strand","mseq")
param <- ScanBamParam(which=RangesList(chrX=IRanges(140410000, 140430000)), what=scanBamWhat()) #,what=scanBamFlag(isPaired = TRUE,isProperPair = TRUE))
bam <- scanBam(bamFile, param=param)
# you need to extract the first item of the list to access the data, I am not sure why yet
b <- bam[[1]]
 # count to which chromosome the mates map and sort
print(sort(table(b$mrnm),decreasing=T))
# mate maps to chromosome 3
i <- which(b$mrnm=='chr3')
b$pos[i]
as.character(b$seq[i])
