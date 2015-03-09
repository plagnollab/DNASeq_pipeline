#!/usr/bin/env Rscript


# Filtering of variants based on annotation

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tools))
suppressPackageStartupMessages(library(xtable))

option_list <- list(
    #make_option(c('--file'), help=''),
    make_option(c('--chr'), help='')
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)


chr <- opt$chr

#dim( d <- read.csv(opt$file) ) 
dim( d <- read.csv(sprintf('VEP_%s-annotations.csv',chr)) ) 

cadd.thresh <- 10
pct.thresh <- .0001

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

# Tony told me that if there is a stop or indel then CADD score is NA
message(sprintf('CADD score > %d or NA',cadd.thresh))
print(table( cadd.filter <- is.na(d$CADD) | (d$CADD > cadd.thresh) ))


# Filter on allele frequency 
# keep if the freq is rare (less than T or  more than 1-T)
# keep if the freq in the AJ is less than 0.05
# (but what if a common variant in the AJ is assoc?)

#message('variant is not in a transcript then ignore')
#print(dim( d <- d[-which(is.na(d$Feature_type)),] ))

message(sprintf('rare (less than %.2f pct or more than %.2f pct) in non-finnish europeans',pct.thresh*100,1-pct.thresh*100))
print(table(exac_nfe.filter <- (d$EXAC_NFE < pct.thresh & (d$HOM==1|d$HET==1)) | (d$EXAC_NFE > (1-pct.thresh) & d$WT==1)  ))


message(sprintf('rare (less than %.2f pct or more than %.2f pct) in americans',pct.thresh*100,1-pct.thresh*100))
print(table(exac_amr.filter <- (d$EXAC_AMR < pct.thresh & (d$HOM==1|d$HET==1)) | (d$EXAC_AMR > (1-pct.thresh) & d$WT==1) ))

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

monogenic <- read.table('/cluster/project8/IBDAJE/monogenic_v2.txt')[,1]
#flag any variant that appears in a gene associated with the monogenic form of the disease
print(dim(mono <- do.call('rbind', lapply( monogenic, function(mono) d[which(mono==d$SYMBOL),] ))))
print(dim(mono <- mono[ grepl('start|stop|splice|frameshift|missense_variant|stop_gained', mono$Consequence), ]))

print(dim(v <- d[ (exac_nfe.filter & exac_amr.filter & cadd.filter) & (csq.filter | carol.filter | condel.filter),]))

print(dim(v <- rbind(mono,v)))

dim( geno.depth <- read.csv(sprintf('VEP_%s-genotypes_depth.csv',chr)) ) 

v <- merge(v,geno.depth,by='VARIANT_ID',all.x=TRUE,all.y=FALSE)

v <- v[,-which(colnames(v) %in% c('DISTANCE','CADD_RAW','CADD_PHRED','MISS','AF'))]

write.csv(v[order(v$CADD,decreasing=TRUE),], file=sprintf('filtered-%s.csv',chr),row.names=FALSE)



