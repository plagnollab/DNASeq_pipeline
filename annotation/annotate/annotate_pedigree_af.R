#!/usr/bin/env Rscript

err.cat <- function(x)  cat(x, '\n', file=stderr())

suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(optparse, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))
suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(data.table, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))
library(kinship2)
#source('/cluster/project8/IBDAJE/scripts/Families/ped/ped-functions.R')
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
calculate <- function(individuals, group.name){
    #With geno. prefix
    #geno.group <- paste("geno",group,sep=".")
    geno.group <- individuals
    #Genotype columns
    #group.d.cols <- which(gsub('geno\\.','',colnames(d)) %in% geno.group)
    err.cat(group.d.cols <- intersect(paste('geno',individuals,sep='.'),colnames(d)))
    if (!length(group.d.cols)) {
        group.wt <- group.het <- group.hom <- group.miss <- group.af <- group.mf <- rep(NA,nrow(d))
        group.out <- as.data.frame(cbind(group.wt,group.het,group.hom,group.miss,group.af,group.mf))
        names(group.out) <- paste(group.name,c("WT","HET","HOM","MISS","AF","MF"),sep="_")
        return(group.out)
    }
    #Genotypes
    group.d.geno <- as.matrix(d[,group.d.cols],nrow=nrow(d),ncol=ncol(d))
    #Number of WT, HET, HOM, MISS
    group.wt <- apply(group.d.geno,1,function(X){length(which(X==0))})
    group.het <- apply(group.d.geno,1,function(X){length(which(X==1))})
    group.hom <- apply(group.d.geno,1,function(X){length(which(X==2))})
    group.miss <- apply(group.d.geno,1,function(X){length(which(is.na(X)))})
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
    err.cat(individuals <- intersect(as.character(ped$ID[which(ped$Subfamily==subfamily&ped$Affection==1)]),samples))
    err.cat(name <- sprintf("%s%s_controls",unique(ped$FamilyAlias),subfamily))
    controls <- calculate(individuals,name)
    d[,names(controls)] <- controls
    err.cat(individuals <- intersect(as.character(ped$ID[which(ped$Subfamily==subfamily&ped$Affection==2)]),samples))
    err.cat(name <- sprintf("%s%s_cases",unique(ped$FamilyAlias),subfamily))
    cases <- calculate(individuals,name)
    d[,names(cases)] <- cases
}

write.csv(d, file='', quote=FALSE, row.names=FALSE)


