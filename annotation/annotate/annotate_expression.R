#!/usr/bin/env Rscript
err.cat <- function(x)  cat(x, '\n', file=stderr())

### 
err.cat(dim( d <- read.csv(file('stdin'),check.names=FALSE) ))


# it is to use grep here I think because ENSG ids are fixed width (15 char)
# so no risk of partial matching

message('*** BARCODE COMMENT ***')
barcode <- read.table('/cluster/project8/IBDAJE/barcode.txt', header=TRUE)
barcode$Gene <- barcode$ENSG
d$barcode.comment <- ''
for (j in 1:nrow(barcode)) {
    b <- barcode[j,]
    d[grep(b['ENSG'],d$Gene),'barcode.comment'] <- b['COMMENT']

}


message('*** XAVIER TISSUES  ***')
xavier <- read.csv('/goon2/scratch2/vyp-scratch2/reference_datasets/GeneEnrichmentProfiler/gene_tissues_ensg.csv')
d$xavier.tissues <- ''
for (j in 1:nrow(xavier)) {
    b <- xavier[j,]
    d[grep(b['gene'],d$Gene),'xavier.tissues'] <- b['tissues']

}


message('*** IBD_GENE_LIST ***')
ibd.gene.list <- read.csv('/cluster/project8/IBDAJE/variants/IBD_GENE_LIST.csv',header=FALSE)
d$ibd.gene.list <- FALSE
for (j in 1:nrow(ibd.gene.list)) {
    gene <- ibd.gene.list[j,1]
    d[grep(gene,d$Gene),'ibd.gene.list'] <- TRUE
}

write.csv( d , quote=FALSE, file='', row.names=FALSE)



