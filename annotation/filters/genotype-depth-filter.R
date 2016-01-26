#!/usr/bin/env Rscript

#make ped/ped-functions into library
source('/cluster/project8/IBDAJE/Families/ped/ped-functions.R')

err.cat <- function(x)  cat(x, '\n', file=stderr())

### Filtering on genotype depth.
#message('*** GENO.DEPTH FILTERING ***')
#d <- read.csv(file('stdin'))

# Filtering of variants based on annotation 
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tools))
suppressPackageStartupMessages(library(xtable))
suppressPackageStartupMessages(library(kinship2))


option_list <- list(
   make_option(c('--ped'), default='pedigree_details.csv', help=''),
   make_option(c('--input'), default='', help='Genotype file.'),
   make_option(c('--depth'), default='', help='Genotype depth file for filttering.'),
   make_option(c('--output'), default='', help='Genotype fle after filtering.'),
   make_option(c('--avg.geno.depth'), default=NULL, help='Filter out all variants with lesser avg geno.depth'),
   make_option(c('--min.geno.depth'), default=NULL, help='Filter out all variants with lesser min geno.depth')
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

#input
input <- opt$input
if (file_ext(input)=='rds') {
   d <- readRDS(input)
} else if (file_ext(input)=='csv') {
   d <- read.csv(input,check.names=FALSE)
}
rownames(d) <- d$VARIANT_ID

# pedigree
ped <- read.pedigree(opt$ped)
samples <- intersect(ped$ID,colnames(d))

#depth
depth <- opt$depth
if (file_ext(depth)=='rds') {
   depth <- readRDS(depth)
} else if (file_ext(depth)=='csv') {
   depth <- read.csv(depth,check.names=FALSE)
}
rownames(depth) <- depth$VARIANT_ID

#
if (!is.null(opt$min.geno.depth)) {
    min.geno.depth <- as.numeric(opt$min.geno.depth)
    for (i in 1:nrow(d)) {
        d[i,names(which(depth[i,samples]<min.geno.depth))] <- NA
    }
}

#
if (!is.null(opt$avg.geno.depth)) {
    avg.deno.depth <- opt$avg.geno.depth
    d$avg.geno.depth <- rowMeans(d[,grep('depth\\.',colnames(d))])
    message('samples')
    err.cat(length(samples <- grep('geno\\.',colnames(d), value=TRUE)))
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
}

#
output <- opt$output
if (file_ext(output)=='rds') {
   saveRDS(d, file=output)
} else if (file_ext(output)=='csv') {
   write.csv(d, file=output)
}


