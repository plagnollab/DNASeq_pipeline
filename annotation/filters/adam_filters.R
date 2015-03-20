#!/usr/bin/env Rscript


# Filtering of variants based on annotation

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tools))
suppressPackageStartupMessages(library(xtable))

option_list <- list(
    #make_option(c('--file'), help=''),
    make_option(c('--chr'), help=''),
    make_option(c('--variants'), default=NULL, help='')
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)


chr <- opt$chr

#dim( d <- read.csv(opt$file) ) 
dim( d <- read.csv(sprintf('VEP_%s-annotations.csv',chr)) ) 

cadd.thresh <- 10
pct.thresh <- .0001
af.thresh <- .05
miss.thresh <- .05

#rare 0.0001
#moderate 0.025

#message('remove CADD NA')
#print(dim(d <- d[!is.na(d$CADD),]))

message('remove EXAC NA')
print(dim(d <- d[!is.na(d$EXAC_Adj) & !is.na(d$EXAC_NFE) & !is.na(d$EXAC_AMR),]))

# this is only valid for a singleton
#message('remove all where EXAC and AF match')

# ignore the ALL_MISS-F  > .7
#tolerate no missingness
#dim( d <- d[which(d$MISS==0),] )

# if there is a stop or indel then CADD score is NA
message(sprintf('CADD score > %d or NA',cadd.thresh))
print(table( cadd.filter <- is.na(d$CADD) | (d$CADD > cadd.thresh) ))


# Filter on allele frequency 
# keep if the freq is rare (less than T or  more than 1-T)
# keep if the freq in the AJ is less than 0.05
# (but what if a common variant in the AJ is assoc?)

#message('variant is not in a transcript then ignore')
#print(dim( d <- d[-which(is.na(d$Feature_type)),] ))

message(sprintf('EXAC_NFE: rare (less than %.2f pct or more than %.2f pct) in non-finnish europeans',pct.thresh*100,1-pct.thresh*100))
print(table(exac_nfe.filter <- (d$EXAC_NFE < pct.thresh & (d$HOM>=1|d$HET>=1)) | (d$EXAC_NFE > (1-pct.thresh) & d$WT>=1)  ))


message(sprintf('EXAC_AMR: rare (less than %.2f pct or more than %.2f pct) in americans',pct.thresh*100,1-pct.thresh*100))
print(table(exac_amr.filter <- (d$EXAC_AMR < pct.thresh & (d$HOM>=1|d$HET>=1)) | (d$EXAC_AMR > (1-pct.thresh) & d$WT>=1) ))


# Filter on consequence
# missense
# We are interested in damaging variants:
# frameshift, stop or splice
message( 'start|stop|splice|frameshift|missense_variant|stop_gained in consequence field')
#csq <- unlist(lapply(strsplit(d$Consequence, ';'), function(x) if (length(x)>0 && !is.na(x)) x))
#table(gsub('\\(.*\\)', '', csq))
print(table(csq.filter <- grepl('start|stop|splice|frameshift|missense_variant|stop_gained', d$Consequence)))

# CAROL Deleterious
message('CAROL deleterious')
#carol <- unlist(lapply(strsplit(d$CAROL, ';'), function(x) if (length(x)>0 && !is.na(x)) x))
#table(gsub('\\(.*\\)', '', carol))
print(table(carol.filter <- grepl('Deleterious',d$CAROL)))

# Condel deleterious
message('Condel deleterious')
#condel <- unlist(lapply(strsplit(d$Condel, ';'), function(x) if (length(x)>0 && !is.na(x)) x))
#table(gsub('\\(.*\\)', '', condel))
print(table(condel.filter <- grepl('deleterious', d$Condel)))

# very different from freq in  AJ controls
#controls <- read.csv(sprintf('/cluster/project8/IBDAJE/batches_from_uclex/controls/VEP_%s-annotations.csv',opt$chr))
#ctl.freq <- (controls$HET+controls$HOM)/50
#freq <- (controls$HET+controls$HOM)/50
#head( d[order(d$HOM,decreasing=T),] ) 
#print(dim(d))

dim(d <- d[,-which(colnames(d) %in% c('WT','HET','HOM','MISS'))])

dim( geno <- read.csv(sprintf('VEP_%s-genotypes.csv',chr)) ) 
dim(geno.counts <- t(apply(geno[,-1], 1, function(x) c(length(which(x==0)), length(which(x==1)), length(which(x==2)), length(which(is.na(x)))) )))
colnames(geno.counts) <- c('WT','HET','HOM','MISS')
rownames(geno.counts) <- geno[,1]

d <- cbind(d, geno.counts[ d$VARIANT_ID, ])

print(dim(v <- d[ (exac_nfe.filter & exac_amr.filter & cadd.filter) & (csq.filter | carol.filter | condel.filter),]))

if ( !is.null(opt$variants) ) {
    #monogenic <- read.table('/cluster/project8/IBDAJE/variants/monogenic_v2.txt')[,1]
    monogenic <- read.table(opt$variants)[,1]
    #flag any variant that appears in a gene associated with the monogenic form of the disease
    print(dim(mono <- do.call('rbind', lapply( monogenic, function(mono) d[which(mono==d$SYMBOL),] ))))
    print(dim(mono <- mono[ grepl('start|stop|splice|frameshift|missense_variant|stop_gained', mono$Consequence), ]))
    print(dim(v <- rbind(mono,v)))
}


dim( geno.depth <- read.csv(sprintf('VEP_%s-genotypes_depth.csv',chr)) ) 
#print(dim(v <- merge(v,geno.depth,by='VARIANT_ID',all.x=TRUE,all.y=FALSE)))

v <- v[,-which(colnames(v) %in% c('DISTANCE','CADD_RAW','CADD_PHRED','SOMATIC','CLIN_SIG','Feature_type'))]
v <- v[,-grep('ONEKG',colnames(v))]
v <- v[,-grep('ESP',colnames(v))]

message(sprintf('AFR_MAF: rare (less than %.2f pct or more than %.2f pct) in africans',pct.thresh*100,1-pct.thresh*100))
print(table(afr_maf.filter <- is.na(v$AFR_MAF) | (v$AFR_MAF < pct.thresh & (v$HOM>=1|v$HET>=1)) | (v$AFR_MAF > (1-pct.thresh) & v$WT>=1) ))

message(sprintf('AF: common (more than %.2f pct)',af.thresh))
print(table(af.filter <- v$AF > af.thresh))

message(sprintf('MISS: less than %.2f pct missingness',miss.thresh))
print(samples.number <- unique(v$WT+v$HET+v$HOM+v$MISS))
print(quantile(miss <- v$MISS/samples.number, probs=seq(0,1,.1)))
print(table(miss.filter <- miss<miss.thresh))

print(dim(v <- v[miss.filter & afr_maf.filter & af.filter,]))

write.csv( v[order(v$AF,decreasing=TRUE),] , quote=FALSE, file=sprintf('filtered-%s.csv',chr),row.names=FALSE)



