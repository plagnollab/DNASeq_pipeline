library(Rsamtools)
 
#bamFile <- "C:/Users/Jessica Gardner/Desktop/LP0022129-DNA_A01.bam"
 
print(files <- grep('.bam$',list.files(path="Z:/WGS/data from Niko for Jess 25.08.15/JessicaGardner/",pattern=".*.bam",full.names = TRUE),value = T))
 
 
for (f in files) {
  print(f)
  if (file.info(f)$size<100) next
  bams <- scanBam(f)
  b <- bam[[1]]
print(sort(table(b$mrnm),decreasing=T))
}
 
bamFile <- "Z:/WGS/bam files/X-LINKED RP/140714_H00TK_PF-19956/140714_H00TK_PF-19956_sorted_unique.bam"
bamFile <- "Z:/WGS/bam files/X-LINKED RP/140724_H02J1_MH-901804/140724_H02J1_MH-901804_sorted_unique.bam"
which <- RangesList(chrX=IRanges(140410000, 140430000))
#what <- c("rname", "strand", "pos", "qwidth", "seq", "mrnm", "mpos","strand","mseq")
what <- scanBamWhat()
param <- ScanBamParam(which=which, what=what)#, what=scanBamFlag(isPaired = TRUE,isProperPair = TRUE))
bam <- scanBam(bamFile, param=param)
b <- bam[[1]]
print(sort(table(b$mrnm),decreasing=T))
i <- which(b$mrnm=='chr3')
b$pos[i]
as.character(b$seq[i])
