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
    make_option(c('--subfamily'), default=NULL, help='CADD score threshold')
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

pedigree <- opt$pedigree
subfamily <- opt$subfamily

d <- read.csv(file('stdin'))


