#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tools))
suppressPackageStartupMessages(library(xtable))
suppressPackageStartupMessages(library(knitr))

option_list <- list(
    make_option(c('--dir'), help=''),
    make_option(c('--chr'), help='')
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)
dir <- opt$dir
chromosome <- opt$chr

#set working dir
setwd(dir)
dir <- '.'

dim(ann <- read.csv(file.path(dir,sprintf('VEP_%s-annotations.csv',chromosome))))
dim(geno <- read.csv(file.path(dir,sprintf('VEP_%s-genotypes.csv',chromosome))))
dim(genodepth <- read.csv(file.path(dir,sprintf('VEP_%s-genotypes_depth.csv',chromosome))))

#feature type
ann$Feature_type[which(is.na(ann$Feature_type))] <- ''

# samples
length(samples <- colnames(geno)[-1])
genodepth$avg.depth <- rowMeans(genodepth[,-1])
positions <- as.numeric(sapply((strsplit(ann$VARIANT_ID,'_')),'[[',2))

# 1Mb sliding window
w <- 10^6
N <- (ceiling(max(positions)/w))

# genotype depth per window
# number of variants per window
var.count <- vector()
# proportion of variants with annotation MAF per window
#maf.count <- vector()
# average depth
avg.depth <- vector()
#
avg.geno <- matrix(0,N,length(samples))
#
for (i in 1:N) {
    r <- which(i*w < positions & positions < (i+1)*w)
    var.count[[i]] <- length(which(!is.na(ann[r,'GMAF'])))
    #maf.count[[i]] <- length(which()) / length(r)
    avg.depth[[i]] <- mean(genodepth[r,'avg.depth'])
    for (j in 1:length(samples)) {
        avg.geno[i,j] <- mean(geno[r,samples[[j]]],na.rm=T)
    }
}

variants <- sapply( 1:N, function(i) tabulate( 1+as.numeric(!is.na(ann[r <- which(i*w < positions & positions < (i+1)*w),'GMAF'])), nbins=2 ) )
rowSums(variants)
colSums(variants)


ESP <- c('ESP_EA', 'ESP_AA')
EXAC <- c('EXAC_AFR', 'EXAC_AMR', 'EXAC_Adj', 'EXAC_EAS', 'EXAC_FIN', 'EXAC_NFE', 'EXAC_OTH', 'EXAC_SAS')
UCL <- c('UCLEX')
ONEKG <- c('ONEKG_EUR', 'ONEKG_AFR', 'ONEKG_AMR', 'ONEKG_ASN')
AF <- c('GMAF',ESP,EXAC, UCL,ONEKG)

initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)

# PDF output: not great for large tables
#output <- file.path(dir,sprintf('VEP_%s.tex',chromosome))
# pdf figures are too large, png woud be better but does not work on cluster nodes
# also tables split across pages
#Sweave( file=file.path(script.basename,"statsVEP.Rnw"), output=output, debug=TRUE)
#texi2pdf( output, clean=TRUE )

# HTML output
output <- file.path(dir,sprintf('VEP_%s.html',chromosome))
text <- knit_expand(file=file.path(script.basename,'statsVEP.Rhtml'),chromosome=chromosome)
knit(text=text,output=output)


