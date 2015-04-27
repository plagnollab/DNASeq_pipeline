#!/usr/bin/env Rscript

# OR filtering for GO terms, gene expression and IBD gene lists.

err.cat <- function(x)  cat(x, '\n', file=stderr())

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tools))
suppressPackageStartupMessages(library(xtable))

option_list <- list(
    make_option(c('--go'), default='/cluster/project8/vyp/AdamLevine/annotations/GO/keep_refined_exact_match_with_header.csv', help=''),
    make_option(c('--variants'), default='/cluster/project8/IBDAJE/variants/monogenic_v2.txt', help='')
)


option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

#go.terms <- '/cluster/project8/vyp/AdamLevine/annotations/GO/keep_refined_exact_match_with_header.csv'
go.terms <- opt$go
err.cat('dim of GO terms')
err.cat( dim(go.term <- read.csv(go.terms,header=TRUE)) )
#monogenic.variants <- '/cluster/project8/IBDAJE/variants/monogenic_v2.txt'
monogenic.variants <- opt$variants 
err.cat('dim of GO terms')
err.cat(dim(monogenic <- read.table(monogenic.variants)[,1]))

err.cat(dim( d <- read.csv(file('stdin')) ))

#d$GO <- gsub('positive_regulation_of_synaptic_transmission&_glutamatergic','positive_regulation_of_synaptic_transmission_glutamatergic', d$GO)
#d$GO <- gsub('regulation_of_transcription&_DNA-templated', 'regulation_of_transcription_DNA-templated', d$GO)

## Filtering based on go terms.
message('*** GO FILTERING ***')
go.regex <- paste(go.term$go,collapse='|')
indexes <- grep(go.regex,d$GO)
err.cat(length(indexes <- sort(unique(indexes))))
go.filter  <- grepl(go.regex,d$GO)
X <- d[go.filter,]

#cleanup GO terms
for (i in 1:nrow(X)) {
     y <- unlist(strsplit(X$GO[i],'&'))
     X$GO[i] <- paste(unique(grep(paste(go.term$go,collapse='|'),y,value=TRUE)),collapse='&',sep='&')
}
#update the GO terms to only include relevant ones
d[go.filter,'GO'] <- X$GO

# Filtering of list of variants 
message('*** VARIANT LIST FILTERING ***')

#flag any variant that appears in a gene associated with the monogenic form of the disease
monogenic.filter <- grepl(paste(monogenic,collapse='|'),d$SYMBOL)


ibd.gene.list.filter <- d[,'ibd.gene.list']

barcode.comment.filter <- grepl('bowel|immune',d[,'barcode.comment'])

xavier.tissues.filter <- d[,'xavier.tissues']!=''

err.cat(table(f <- go.filter | ibd.gene.list.filter | barcode.comment.filter | xavier.tissues.filter ))

d <- d[f,]

write.csv(d, quote=FALSE, file='',row.names=FALSE)

