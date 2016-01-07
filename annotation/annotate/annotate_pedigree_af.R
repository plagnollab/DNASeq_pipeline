#!/usr/bin/env Rscript

err.cat <- function(x)  cat(x, '\n', file=stderr())

suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(optparse, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))
suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(data.table, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))
source('/cluster/project8/IBDAJE/scripts/Families/ped/ped-functions.R')

ped <- read.pedigree('pedigree_details.csv')

err.cat('number of sufamilies:')
err.cat(length(subfamilies <- na.omit(unique(ped$Subfamily))))

#Load filtered variants with genotypes
#d <- read.csv("go-csq-af-filtered.csv")
d <- as.data.frame(fread('VEP_22.csv'),check.names=FALSE)
#d <- as.data.frame(fread('cat /dev/stdin'),check.names=FALSE)

err.cat(samples <- gsub('geno\\.','',grep('geno',colnames(d),value=TRUE)))
#colnames(d) <- gsub('geno\\.','',colnames(d))

#Function to subset individuals and do calculations
#group: list of individuals of interest, group.name: subfamily name
calculate <- function(group, group.name){
    #With geno. prefix
    #geno.group <- paste("geno",group,sep=".")
    geno.group <- group
    #Genotype columns
    group.d.cols <- which(gsub('geno\\.','',colnames(d)) %in% geno.group)
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
    names(group.out) <- paste(group.name,c("WT","HET","HOM","MISS","AF","MF"),sep="_")
    return(group.out)

}

##Group A0 only
#Names of individuals of interest
for (subfamily in subfamilies) {
    err.cat(subfamily)
    group.id <- as.character(ped$ped[which(ped$Subfamily==subfamily&ped$Affection==1)])
    err.cat(name <- sprintf("%s%s_controls",unique(ped$FamilyAlias),subfamily))
    d[,name] <- calculate(group.id,name)
    group.id <- as.character(ped$ped[which(ped$Subfamily==subfamily&ped$Affection==2)])
    err.cat(name <- sprintf("%s%s_cases",unique(ped$FamilyAlias),subfamily))
    d[,name] <- calculate(group.id,name)
}

write.csv(d, file='', quote=FALSE, row.names=FALSE)


