#!/usr/bin/env Rscript

# Filtering of variants based on annotation 
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tools))
suppressPackageStartupMessages(library(xtable))

option_list <- list(
    make_option(c('--chr'), default=NULL, help='')
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)
chr <- opt$chr

if (!is.null(chr)) {
    # combine annotation, genotypes and genotypes_depth per chromosome
    ann <- read.csv(sprintf('VEP_%s-annotations.csv',chr),check.names=FALSE)
    geno <- read.csv(sprintf('VEP_%s-genotypes.csv',chr),check.names=FALSE)
    print(dim(geno <- geno[,-which('VARIANT_ID'==names(geno))]))
    names(geno) <- paste('geno',names(geno),sep='.')
    geno.depth <- read.csv(sprintf('VEP_%s-genotypes_depth.csv',chr),check.names=FALSE)
    print(dim(geno.depth <- geno.depth[,-which('VARIANT_ID'==names(geno.depth))]))
    names(geno.depth) <- paste('depth',names(geno.depth),sep='.')
    d <- cbind(ann,geno,geno.depth)
    write.csv(d, quote=FALSE, row.names=FALSE, file=sprintf('VEP_%s.csv',opt$chr))
} else {
    # works on all chromosomes to create one giant file
    # first combine all chromosomes into one file for annotations, genotypes and genotypes_depth
    for (file in c('annotations','genotypes','genotypes_depth')) {
        print(dim(d <- do.call('rbind', lapply(c(seq(1,22),'X'), function(chr) {
            print(chr)
            print(dim(x <- read.csv(sprintf('VEP_%s-%s.csv',chr,file))),check.names=FALSE)
            print(dput(colnames(x)))
            return(x)
        }
        )) ))
        write.csv(d,file=sprintf('VEP-%s.csv',file),row.names=FALSE,quote=FALSE)
    }
    # then combine annotations, genotypes nd genotypes_depth
    ann <- read.csv('VEP-annotations.csv',check.names=FALSE)
    geno <- read.csv('VEP-genotypes.csv',check.names=FALSE)
    print(dim(geno <- geno[,-which('VARIANT_ID'==colnames(geno))]))
    colnames(geno) <- paste('geno',colnames(geno),sep='.')
    geno.depth <- read.csv('VEP-genotypes_depth.csv',check.names=FALSE)
    print(dim(geno.depth <- geno.depth[,-which('VARIANT_ID'==colnames(geno.depth))]))
    colnames(geno.depth) <- paste('depth',colnames(geno.depth),sep='.')
    d <- cbind(ann,geno,geno.depth)
    write.csv(d, quote=FALSE, row.names=FALSE, file='VEP.csv')
}



