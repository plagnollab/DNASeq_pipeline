#!/usr/bin/env Rscript
err.cat <- function(x)  cat(x, '\n', file=stderr())

## Filtering based on go terms.
message('*** GO FILTERING ***')

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tools))
suppressPackageStartupMessages(library(xtable))

option_list <- list(
    make_option(c('--go'), default='/cluster/project8/IBDAJE/GO.txt', help='')
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

go <- opt$go

err.cat(dim( d <- read.csv(file('stdin')) ))

d$GO <- gsub('positive_regulation_of_synaptic_transmission&_glutamatergic','positive_regulation_of_synaptic_transmission_glutamatergic', d$GO)
d$GO <- gsub('regulation_of_transcription&_DNA-templated', 'regulation_of_transcription_DNA-templated', d$GO)

err.cat(dim(go.term <- read.table(go,sep=';',header=TRUE)))

#lapply(strsplit(d$GO,'&'), function(x) {
#lapply(strsplit(x,':')
#}

indexes <- c()
for (x in go.term$go) {
    indexes <- c(indexes, grep(x,d$GO))
}

err.cat(length(indexes <- sort(unique(indexes))))

d <- d[indexes,]

#cleanup GO terms
for (i in 1:nrow(d)) {
    y <- unlist(strsplit(d$GO[i],'&'))
    d$GO[i] <- paste(unique(grep(paste(go.term$go,collapse='|'),y,value=TRUE)),collapse='&',sep='&')
}

write.csv(d, quote=FALSE, file='',row.names=FALSE)








