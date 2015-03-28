#!/usr/bin/env Rscript
err.cat <- function(x)  cat(x, '\n', file=stderr())

### Filtering on genotype depth.
message('*** GENO.DEPTH FILTERING ***')
d <- read.csv(file('stdin'))

# Filtering of variants based on annotation 
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tools))
suppressPackageStartupMessages(library(xtable))


option_list <- list(
    make_option(c('--pedigree'), default='DNA_pedigree_details.csv', help='Pedigree containing affection status.'),
    make_option(c('--avg.geno.depth'), default=10, help='Filter out all variants with lesser avg geno.depth'),
    make_option(c('--min.geno.depth'), default=NULL, help='Filter out all variants with lesser min geno.depth')
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

avg.deno.depth <- opt$avg.geno.depth
pedigree <- opt$pedigree

d$avg.geno.depth <- rowMeans(d[,grep('depth\\.',colnames(d))])

message('samples')
err.cat(length(samples <- grep('geno\\.',colnames(d), value=TRUE)))
# pedigree
message('dim of pedigree')
err.cat(nrow(pedigree <- read.csv(pedigree)))
sample.affection <- pedigree[ which( pedigree$uclex.sample %in% samples ), c('uclex.sample','Affection')]
# cases
message('cases')
err.cat(cases <- sample.affection[which(sample.affection$Affection==2),'uclex.sample'])
message('number of cases')
err.cat(length(cases))
# controls
message('controls')
err.cat(controls <- sample.affection[which(sample.affection$Affection==1),'uclex.sample'])
message('number of controls')
err.cat(length(controls))


d$avg.geno.depth <- rowMeans(d[,grep('depth\\.',colnames(d))])
d$ca.avg.geno.depth <- d[,paste('depth',cases,sep='.')]
d$co.avg.geno.depth <- d[,paste('depth',controls,sep='.')]


write.csv( d[order(v$R,decreasing=TRUE),] , quote=FALSE, file='', row.names=FALSE)

