#!/usr/bin/env Rscript

err.cat <- function(x)  cat(x, '\n', file=stderr())

suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(optparse, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))
suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(data.table, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))
source('/cluster/project8/IBDAJE/scripts/Families/ped/ped-functions.R')
pedigree <- read.pedigree('pedigree_details.csv')


### 
err.cat(dim( d <- as.data.frame(fread('cat /dev/stdin')) ))
#
samples <- gsub('geno\\.','',grep('geno',colnames(d),value=TRUE))
#

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

#group: list of individuals of interest, group.name: subfamily name
calculate <- function(group, group.name) { 
    #With geno. prefix
    geno.group <- paste("geno",group,sep=".")
    #Genotype columns
    group.d.cols <- which(names(d) %in% geno.group)
    #Genotypes
    group.d.geno <- d[,group.d.cols]
    #Number of WT, HET, HOM, MISS
    group.wt <- apply(X=group.d.geno,1,function(X){length(which(X==0))})
    group.het <- apply(X=group.d.geno,1,function(X){length(which(X==1))})
    group.hom <- apply(X=group.d.geno,1,function(X){length(which(X==2))})
    group.miss <- apply(X=group.d.geno,1,function(X){length(which(is.na(X)))})
    #Allele frequency
    group.af <- (group.het+group.hom*2)/(2*(group.wt+group.het+group.hom))
    #Mutant frequency (percentage either HET or WT over total)
    group.mf <- (group.het+group.hom)/(group.wt+group.het+group.hom)
    #As data frame
    group.out <- as.data.frame(cbind(group.wt,group.het,group.hom,group.miss,group.af,group.mf))
    names(group.out) = paste(group.name,c("WT","HET","HOM","MISS","AF","MF"),sep="_")
    return(group.out)
}


ca <- calculate(cases, 'cases')
co <- calculate(controls, 'controls')
d <- cbind(d,ca,co)

write.csv( d , quote=FALSE, file='', row.names=FALSE)





