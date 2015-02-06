#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option(c('--annotation'), help=''),
    make_option(c('--genotype'), help=''),
    make_option(c('--genotype_depth'), help='')
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

dim(ann <- read.csv(opt$annotation))
dim(geno <- read.csv(opt$genotype))
dim(genodepth <- read.csv(opt$genotype_depth))

print(summary(ann))

#d$ONEKG_AMR

ESP=c('ESP_EA', 'ESP_AA')
EXAC=c('EXAC_AFR', 'EXAC_AMR', 'EXAC_Adj', 'EXAC_EAS', 'EXAC_FIN', 'EXAC_NFE', 'EXAC_OTH', 'EXAC_SAS')
UCL=c('UCLEX')
ONEKG=c('ONEKG_EUR', 'ONEKG_AFR', 'ONEKG_AMR', 'ONEKG_ASN')

#check that numbers make sense in terms of coding variants per chromosome
table(ann$Feature_type=='Transcript')

# missing AFs
apply(ann[,c('GMAF',ESP,EXAC,UCL,ONEKG)],2,function(x) table(!is.na(x)))


# samples
length(samples <- colnames(geno)[-1])

#Three sanity checks:
#- nb of non reference coding variants in the range of 20K per sample
apply(geno[,samples],2,table)

#- Nb of rare (< 0.01 MAF) coding variants in the range of 500
table(ann$GMAF<.01)


#- Nb of rare homozygous variants (< 0.01 MAF) coding variants should be in the range of 0-10
#- per sample.
sapply(samples, function(sample) table(ann$GMAF<.01 & geno[,sample] == 2))


# are there large regions of missing annotations?

#Check the gap at end of chr6 for 1KG

# no data for chr22 exac, make sure everything makes sense

# remove redundancy for most columns

# maybe zip the tables?


#create filtered sets with coding data to push to IoO
#fix the "?", decide what is 0 and NA- ignore for now

positions <- as.numeric(sapply((strsplit(ann$VARIANT_ID,'_')),'[[',2))

print(plot.file <- gsub('-annotations.csv','.pdf',opt$annotation))

#plots
pdf(plot.file)
#transcript location
plot(NULL,xlim=range(positions),ylim=c(0,length(samples)),xlab='position',ylab='individual',main='Feature type')
for (i in 1:length(samples)) points(cbind(positions,i),col=ifelse(is.na(ann$Feature_type),'black','red'),pch=ifelse(is.na(ann$Feature_type),'.','x'))
#GMAF
plot(positions,ann$GMAF,xlab='position',ylab='GMAF',main='',pch='.')
#HOM
plot(positions,ann$HOM,xlab='position',ylab='HOM',main='',type='l')
#relatedness
plot(hclust(dist(t(geno[,samples]))),main='relatedness',xlab='hierarchical clustering of samples')
dev.off()


