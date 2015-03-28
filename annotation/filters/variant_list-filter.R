#!/usr/bin/env Rscript

err.cat <- function(x)  cat(x, '\n', file=stderr())

# Filtering of list of variants 
message('*** VARIANT LIST FILTERING ***')

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tools))
suppressPackageStartupMessages(library(xtable))

option_list <- list(
    make_option(c('--variants'), default=NULL, help='')
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

message('dim of input file')
err.cat(dim( d <- read.csv(file('stdin')) ))

#monogenic <- read.table('/cluster/project8/IBDAJE/variants/monogenic_v2.txt')[,1]
monogenic <- read.table(opt$variants)[,1]
#flag any variant that appears in a gene associated with the monogenic form of the disease
err.cat(dim(mono <- do.call('rbind', lapply( monogenic, function(mono) d[which(mono==d$SYMBOL),] ))))
err.cat(dim(mono <- mono[ grepl('start|stop|splice|frameshift|missense_variant|stop_gained', mono$Consequence), ]))
err.cat(dim(v <- rbind(mono,v)))

write.csv(v, quote=FALSE, file='', row.names=FALSE)

