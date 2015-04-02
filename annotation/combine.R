#!/usr/bin/env Rscript

# Filtering of variants based on annotation 
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tools))
suppressPackageStartupMessages(library(xtable))

option_list <- list(
    make_option(c('--chr'), help='')
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)
chr <- opt$chr
ann <- read.csv(sprintf('VEP_%s-annotations.csv',chr))
geno <- read.csv(sprintf('VEP_%s-genotypes.csv',chr))
print(dim(geno <- geno[,-which('VARIANT_ID'==names(geno))]))
names(geno) <- paste('geno',names(geno),sep='.')
geno.depth <- read.csv(sprintf('VEP_%s-genotypes_depth.csv',chr))
print(dim(geno.depth <- geno.depth[,-which('VARIANT_ID'==names(geno.depth))]))
names(geno.depth) <- paste('depth',names(geno.depth),sep='.')

d <- cbind(ann,geno,geno.depth)

write.csv(d, quote=FALSE, row.names=FALSE, file=sprintf('VEP_%s.csv',chr))


