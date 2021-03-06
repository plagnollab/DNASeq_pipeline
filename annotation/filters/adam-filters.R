#!/usr/bin/env Rscript
err.cat <- function(x)  cat(x, '\n', file=stderr())

### Series of filters suggested by Adam.
message('*** ADAM FILTERING ***')

# Filtering of variants based on annotation 
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tools))
suppressPackageStartupMessages(library(xtable))


option_list <- list(
    make_option(c('--pedigree'), default='DNA_pedigree_details.csv', help=''),
    make_option(c('--cadd.thresh'), default=10, help='CADD score threshold'),
    make_option(c('--pop.thresh'), default=.05, help='pop freq threshold'),
    make_option(c('--af.thresh'), default=.05, help='af freq threshold'),
    make_option(c('--miss.thresh'), default=.2, help='miss freq threshold')
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

pedigree <- opt$pedigree
cadd.thresh <- opt$cadd.thresh
pop.thresh <- opt$pop.thresh
af.thresh <- opt$af.thresh
miss.thresh <- opt$miss.thresh

d <- read.csv(file('stdin'))

message('samples')
err.cat(samples <- gsub( 'geno\\.', '', grep('geno\\.',colnames(d),value=TRUE) ))
err.cat(length(samples))

#
message('dim of pedigree')
err.cat(nrow(pedigree <- read.csv(pedigree)))
sample.affection <- pedigree[ which( pedigree$uclex.sample %in% samples ), c('uclex.sample','Affection')]

message('cases')
err.cat(cases <- sample.affection[which(sample.affection$Affection==2),'uclex.sample'])
message('number of cases')
err.cat(length(cases))

#
message('controls')
err.cat(controls <- sample.affection[which(sample.affection$Affection==1),'uclex.sample'])
message('number of controls')
err.cat(length(controls))


#rare 0.0001
#moderate 0.025

#message('remove CADD NA')
#err.cat(dim(d <- d[!is.na(d$CADD),]))

message('remove EXAC NA')
err.cat(nrow(d <- d[!is.na(d$EXAC_Adj) & !is.na(d$EXAC_NFE) & !is.na(d$EXAC_AMR),]))

# this is only valid for a singleton
#message('remove all where EXAC and AF match')

# ignore the ALL_MISS-F  > .7
#tolerate no missingness
#dim( d <- d[which(d$MISS==0),] )

# if there is a stop or indel then CADD score is NA
message(sprintf('CADD score > %d or NA',cadd.thresh))
err.cat(table( cadd.filter <- is.na(d$CADD) | (d$CADD > cadd.thresh) ))


# Filter on allele frequency 
# keep if the freq is rare (less than T or  more than 1-T)
# keep if the freq in the AJ is less than 0.05
# (but what if a common variant in the AJ is assoc?)

#message('variant is not in a transcript then ignore')
#err.cat(dim( d <- d[-which(is.na(d$Feature_type)),] ))

message(sprintf('EXAC_NFE: rare (less than %.2f pct or more than %.2f pct) in non-finnish europeans',pop.thresh*100,100-pop.thresh*100))
err.cat(table(exac_nfe.filter <- (d$EXAC_NFE < pop.thresh & (d$HOM>=1|d$HET>=1)) | (d$EXAC_NFE > (1-pop.thresh) & d$WT>=1)  ))

message(sprintf('EXAC_AMR: rare (less than %.2f pct or more than %.2f pct) in americans',pop.thresh*100,100-pop.thresh*100))
err.cat(table(exac_amr.filter <- (d$EXAC_AMR < pop.thresh & (d$HOM>=1|d$HET>=1)) | (d$EXAC_AMR > (1-pop.thresh) & d$WT>=1) ))

# must filter on onekg allele freq because of 17_26699195_C_CG

# Filter on consequence
# missense
# We are interested in damaging variants:
# frameshift, stop or splice
message( 'start|stop|splice|frameshift|missense_variant|stop_gained in consequence field')
#csq <- unlist(lapply(strsplit(d$Consequence, ';'), function(x) if (length(x)>0 && !is.na(x)) x))
#table(gsub('\\(.*\\)', '', csq))
err.cat(table(csq.filter <- grepl('start|stop|splice|frameshift|missense_variant|stop_gained', d$Consequence)))

# CAROL Deleterious
message('CAROL deleterious')
#carol <- unlist(lapply(strsplit(d$CAROL, ';'), function(x) if (length(x)>0 && !is.na(x)) x))
#table(gsub('\\(.*\\)', '', carol))
err.cat(table(carol.filter <- grepl('Deleterious',d$CAROL)))

# Condel deleterious
message('Condel deleterious')
#condel <- unlist(lapply(strsplit(d$Condel, ';'), function(x) if (length(x)>0 && !is.na(x)) x))
#table(gsub('\\(.*\\)', '', condel))
err.cat(table(condel.filter <- grepl('deleterious', d$Condel)))

# very different from freq in  AJ controls
#controls <- read.csv(sprintf('/cluster/project8/IBDAJE/batches_from_uclex/controls/VEP_%s-annotations.csv',opt$chr))
#ctl.freq <- (controls$HET+controls$HOM)/50
#freq <- (controls$HET+controls$HOM)/50
#head( d[order(d$HOM,decreasing=T),] ) 
#print(dim(d))

#dim(d <- d[,-which(colnames(d) %in% c('WT','HET','HOM','MISS'))])

#dim( geno <- read.csv(sprintf('VEP_%s-genotypes.csv',chr)) ) 
geno <- d[,grep('geno\\.',colnames(d))]

controls.geno.counts <- t(apply(geno[,paste('geno',controls,sep='.')], 1, function(x) c(length(which(x==0)), length(which(x==1)), length(which(x==2)), length(which(is.na(x)))) ))
colnames(controls.geno.counts) <- c('co.WT','co.HET','co.HOM','co.MISS')
rownames(controls.geno.counts) <- geno[,1]

cases.geno.counts <- t(apply(geno[,paste('geno',cases,sep='.')], 1, function(x) c(length(which(x==0)), length(which(x==1)), length(which(x==2)), length(which(is.na(x)))) ))
colnames(cases.geno.counts) <- c('ca.WT','ca.HET','ca.HOM','ca.MISS')
rownames(cases.geno.counts) <- geno[,1]

d <- cbind(d, controls.geno.counts, cases.geno.counts)

d$ca.AF <- (d[,'ca.HET']+2*d[,'ca.HOM'])/(2*length(cases))
d$co.AF <- (d[,'co.HET']+2*d[,'co.HOM'])/(2*length(controls))

err.cat(dim(v <- d[ (exac_nfe.filter & exac_amr.filter & cadd.filter) & (csq.filter | carol.filter | condel.filter),]))

#dim( geno.depth <- read.csv(sprintf('VEP_%s-genotypes_depth.csv',chr)) ) 
geno.depth <- d[,grep('depth\\.',colnames(d))]
#print(dim(v <- merge(v,geno.depth,by='VARIANT_ID',all.x=TRUE,all.y=FALSE)))

v <- v[,-which(colnames(v) %in% c('DISTANCE','CADD_RAW','CADD_PHRED','SOMATIC','CLIN_SIG','Feature_type'))]
#v <- v[,-grep('ONEKG',colnames(v))]
v <- v[,-grep('ESP',colnames(v))]

message(sprintf('AFR_MAF: rare (less than %.2f pct or more than %.2f pct) in africans',pop.thresh*100,100-pop.thresh*100))
err.cat(table(afr_maf.filter <- is.na(v$AFR_MAF) | (v$AFR_MAF < pop.thresh & (v$HOM>=1|v$HET>=1)) | (v$AFR_MAF > (1-pop.thresh) & v$WT>=1) ))

#v[,grep('ONEKG_', colnames(v))]

message(sprintf('AF: common (more than %.2f pct)',af.thresh))
err.cat(table(af.filter <- v$AF > af.thresh))

message(sprintf('MISS: less than %.2f pct missingness',miss.thresh))
err.cat(samples.number <- unique(v$WT+v$HET+v$HOM+v$MISS))
err.cat(quantile(miss <- v$MISS/samples.number, probs=seq(0,1,.1)))
err.cat(table(miss.filter <- miss<miss.thresh))

err.cat(dim(v <- v[miss.filter & afr_maf.filter & af.filter,]))

v$R <- v$ca.AF/v$co.AF

write.csv( v[order(v$ca.AF,v$R,decreasing=TRUE),] , quote=FALSE, file='', row.names=FALSE)



