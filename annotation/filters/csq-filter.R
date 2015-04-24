#!/usr/bin/env Rscript
err.cat <- function(x)  cat(x, '\n', file=stderr())

### 
message('*** CADD FILTERING ***')
d <- read.csv(file('stdin'))

# Filtering of variants based on annotation 
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tools))
suppressPackageStartupMessages(library(xtable))


option_list <- list(
    make_option(c('--cadd.thresh'), default=20, help='CADD score threshold'),
    make_option(c('--csq.filter'), default='start|stop|splice|frameshift|stop_gained', help='csq field'),
    make_option(c('--carol.filter'), default='Deleterious', help='CAROL'),
    make_option(c('--condel.filter'), default='deleterious', help='CAROL')
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)


message('samples')
err.cat(length(samples <- grep('geno\\.',colnames(d), value=TRUE)))


# if there is a stop or indel then CADD score is NA
if (!is.null(opt$cadd.thresh)) {
    cadd.thresh <- opt$cadd.thresh
    message(sprintf('CADD score > %d or NA',cadd.thresh))
    err.cat(table( cadd.filter <- is.na(d$CADD) | (d$CADD > cadd.thresh) ))
    #d <- d[cadd.filter,]
}


# Filter on consequence, these need to ORed rather than ANDed
# missense
# We are interested in damaging variants:
# frameshift, stop or splice
if (!is.null(opt$csq.filter)) {
    csq.filter <- opt$csq.filter
    message( sprintf('%s in consequence field',opt$csq.filter) )
    err.cat(table(csq.filter <- grepl(opt$csq.filter, d$Consequence)))
    #d <- d[csq.filter,]
}

# CAROL Deleterious
if (!is.null(opt$carol.filter)) {
    carol.filter <- opt$carol.filter
    message( sprintf('CAROL %s',opt$carol.filter) )
    err.cat(table(carol.filter <- grepl(opt$carol.filter,d$CAROL)))
    #d <- d[carol.filter,]
}

# Condel deleterious
if (!is.null(opt$condel.filter)) {
    condel.filter <- opt$condel.filter
    message(sprintf('Condel %s',opt$condel.filter))
    err.cat(table(condel.filter <- grepl(opt$condel.filter, d$Condel)))
    #d <- d[condel.filter,]
}

f1 <- csq.filter
# missense and (carol or condel deleterious)
f2 <- ( grepl('missense_variant',d$Consequence) & (carol.filter|condel.filter|cadd.filter) )

d <- d[ f1 | f2, ]

write.csv( d[order(d$CADD,decreasing=TRUE),] , quote=FALSE, file='', row.names=FALSE)

