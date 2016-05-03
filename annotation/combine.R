#!/usr/bin/env Rscript

# Filtering of variants based on annotation 
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
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
    print('annotations')
    print(dim(ann <- as.data.frame(fread(sprintf('VEP_%s-annotations.csv',chr),header=TRUE))))
    print('genotypes')
    print(dim(geno <- as.data.frame(fread(sprintf('VEP_%s-genotypes.csv',chr),header=TRUE))))
    #print(dim(geno <- geno[,-which('VARIANT_ID'==names(geno))]))
    names(geno) <- paste('geno',names(geno),sep='.')
    print('genotypes depth')
    print(dim(geno.depth <- as.data.frame(fread(sprintf('VEP_%s-genotypes_depth.csv',chr),header=TRUE))))
    #print(dim(geno.depth <- geno.depth[,-which('VARIANT_ID'==names(geno.depth))]))
    names(geno.depth) <- paste('depth',names(geno.depth),sep='.')
    # ignore the first column which is variant id
    d <- cbind(ann,geno[,-1],geno.depth[,-1])
    write.csv(d, quote=FALSE, row.names=FALSE, file=sprintf('VEP_%s.csv',chr))
} else {
    # works on all chromosomes to create one giant file
    # first combine all chromosomes into one file for annotations, genotypes and genotypes_depth
    for (file in c('annotations','genotypes','genotypes_depth')) {
        print(dim(d <- do.call('rbind', lapply(c(seq(1,22),'X'), function(chr) {
            print(chr)
            print(dim(x <- as.data.frame(fread(sprintf('VEP_%s-%s.csv',chr,file),header=TRUE))))
            print(dput(colnames(x)))
            return(x)
        }
        )) ))
        write.csv(d,file=sprintf('VEP-%s.csv',file),row.names=FALSE,quote=FALSE)
    }
    # then combine annotations, genotypes and genotypes_depth
    ann <- as.data.frame(fread('VEP-annotations.csv',header=TRUE))
    geno <- as.data.frame(fread('VEP-genotypes.csv',header=TRUE))
    print(dim(geno <- geno[,-which('VARIANT_ID'==colnames(geno))]))
    colnames(geno) <- paste('geno',colnames(geno),sep='.')
    geno.depth <- as.data.frame(fread('VEP-genotypes_depth.csv',header=TRUE))
    print(dim(geno.depth <- geno.depth[,-which('VARIANT_ID'==colnames(geno.depth))]))
    colnames(geno.depth) <- paste('depth',colnames(geno.depth),sep='.')
    d <- cbind(ann,geno,geno.depth)
    write.csv(d, quote=FALSE, row.names=FALSE, file='VEP.csv')
}



