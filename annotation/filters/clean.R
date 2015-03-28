#!/usr/bin/env Rscript
err.cat <- function(x)  cat(x, '\n', file=stderr())

### Series of filters suggested by Adam.
message('*** CLEAN ***')

# Filtering of variants based on annotation 
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tools))

message('dim of input file')
err.cat(dim( v <- read.csv(file('stdin')) ))

#v <- v[,-which(colnames(v) %in% c('DISTANCE','CADD_RAW','CADD_PHRED','SOMATIC','CLIN_SIG','Feature_type'))]
#v <- v[,-grep('ESP',colnames(v))]
v <- v[,which(!grepl('STRAND|SYMBOL_SOURCE|HGNC_ID|CANONICAL',colnames(v)))]
v <- v[,which(!grepl('geno',colnames(v)))]
v <- v[,which(!grepl('depth',colnames(v)))]
#
#cleanup GO terms (is done in the GO filter)

write.csv(v[order(v$R,decreasing=TRUE),], quote=FALSE, file='', row.names=FALSE)

