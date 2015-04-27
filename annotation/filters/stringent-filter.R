#!/usr/bin/env Rscript

#Make a more stringent filtering sequence (in addition to the existing one). Perhaps the following?
# GO & expression
# Condel & CAROL deleterious
# Frequency <0.01

err.cat <- function(x)  cat(x, '\n', file=stderr())

## Filtering based on go terms.
message('*** STRINGENT FILTERING ***')

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tools))
suppressPackageStartupMessages(library(xtable))

option_list <- list(
    make_option(c('--go'), default='/cluster/project8/vyp/AdamLevine/annotations/GO/keep_refined_exact_match_with_header.csv', help='')
)


option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

go <- opt$go

err.cat(dim( d <- read.csv(file('stdin')) ))

err.cat('dim of GO terms')
err.cat( dim(go.term <- read.csv(go,header=TRUE)) )

go.regex <- paste(go.term$go,collapse='|')
indexes <- grep(go.regex,d$GO)
err.cat(length(indexes <- sort(unique(indexes))))

go.filter  <- indexes
d <- d[go.filter,]

#cleanup GO terms
for (i in 1:nrow(d)) {
     y <- unlist(strsplit(d$GO[i],'&'))
     d$GO[i] <- paste(unique(grep(paste(go.term$go,collapse='|'),y,value=TRUE)),collapse='&',sep='&')
}


write.csv(d, quote=FALSE, file='',row.names=FALSE)


