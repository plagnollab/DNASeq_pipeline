#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tools))
suppressPackageStartupMessages(library(xtable))
suppressPackageStartupMessages(library(knitr))

option_list <- list(
    make_option(c('--generate'), action='store_true', default=FALSE, help='generate the gene_tissues.csv file'),
    make_option(c('--annotate'), default='', help='file to annotate')
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

BASE.DIR <- '/goon2/scratch2/vyp-scratch2/reference_datasets/GeneEnrichmentProfiler'


###
if (opt$generate) {
setwd(BASE.DIR) 
mapping <- read.table('U133A.YBenitaMapping.txt')
colnames(mapping) <- c('ProbeID', 'gene') 
tissues <- read.table('tissues_keep.txt',sep='\t')
colnames(tissues) <- c('tissue','keep') 
enrichment <- read.table('EnrichmentValues.Primary-U133A.txt',sep='\t',header=TRUE, as.is=TRUE, check.names=FALSE)
keep.tissues <- tissues[which(tissues$keep==1),'tissue']
enrichment <- enrichment[, c('ProbeID', keep.tissues) ]
d <- merge(mapping,enrichment)
probeid.gene <- d[,c('ProbeID','gene')]
d <- d[,-(1:2)] 
d[d<=0] <- FALSE
d[d>0] <- TRUE 
#collapse by gene name 
X <- do.call('rbind', by(d, probeid.gene[,'gene'], function(x) colSums(x) )) 
X <- X[-which(rowSums(X)==0),] 
Y <- apply(X,1,function(x) paste(colnames(X)[which(x==1)],collapse=';')) 
Y <- data.frame(gene=names(Y), tissues=Y)
write.csv(Y, file='gene_tissues.csv',quote=FALSE, row.names=FALSE)
}

###
if (!is.null(opt$annotate)) {
annotate.file <- opt$annotate
gene.tissues <- read.csv(file.path(BASE.DIR,'gene_tissues.csv'))
colnames(gene.tissues) <- c('SYMBOL','tissues')
d <- read.csv(annotate.file)
d <- merge(d, gene.tissues)
write.csv(d, file='', quote=FALSE, row.names=FALSE)
}




