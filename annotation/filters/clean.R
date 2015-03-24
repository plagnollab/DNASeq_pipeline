#!/usr/bin/env Rscript
err.cat <- function(x)  cat(x, '\n', file=stderr())

### Series of filters suggested by Adam.
message('*** CLEAN ***')

# Filtering of variants based on annotation 
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tools))

message('dim of input file')
err.cat(dim( v <- read.csv(file('stdin')) ))

v <- v[,-which(colnames(v) %in% c('DISTANCE','CADD_RAW','CADD_PHRED','SOMATIC','CLIN_SIG','Feature_type'))]
#v <- v[,-grep('ONEKG',colnames(v))]
v <- v[,-grep('ESP',colnames(v))]
v <- v[,-grep('geno',colnames(v)]
v <- v[,-grep('depth',colnames(v)]

write.csv(v, quote=FALSE, file='', row.names=FALSE)

