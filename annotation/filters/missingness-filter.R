#!/usr/bin/env Rscript
err.cat <- function(x)  cat(x, '\n', file=stderr())

### Series of filters suggested by Adam.
message('*** MISS FILTERING ***')

# Filtering of variants based on annotation 
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tools))
suppressPackageStartupMessages(library(xtable))


option_list <- list(
    make_option(c('--miss.thresh'), default=.2, help='miss freq threshold')
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

v <- read.csv(file('stdin'))

message(sprintf('MISS: less than %.2f pct missingness',miss.thresh))
err.cat(samples.number <- unique(v$WT+v$HET+v$HOM+v$MISS))
err.cat(quantile(miss <- v$MISS/samples.number, probs=seq(0,1,.1)))
err.cat(table(miss.filter <- miss<miss.thresh))
err.cat(table(miss.filter <- miss<miss.thresh))

write.csv( v, quote=FALSE, file='', row.names=FALSE)


