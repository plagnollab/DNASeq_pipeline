#!/usr/bin/env Rscript
err.cat <- function(x)  cat(x, '\n', file=stderr())

## Filtering based on go terms.
message('*** GO FILTERING ***')

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tools))
suppressPackageStartupMessages(library(xtable))

option_list <- list(
    make_option(c('--go'), help='')
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

go <- opt$go

err.cat(dim( d <- read.csv(file('stdin')) ))

d$GO <- gsub('positive_regulation_of_synaptic_transmission&_glutamatergic','positive_regulation_of_synaptic_transmission_glutamatergic', d$GO)

err.cat(dim(go.term <- read.table(go,sep=';',header=TRUE)))

#lapply(strsplit(d$GO,'&'), function(x) {
#lapply(strsplit(x,':')
#}

indexes <- c()
for (x in go.term$go) {
    indexes <- c(indexes, grep(x,d$GO))
}

err.cat(length(indexes <- sort(unique(indexes))))

write.csv(d[indexes,], quote=FALSE, file='',row.names=FALSE)








