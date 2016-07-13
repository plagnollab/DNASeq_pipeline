#!/usr/bin/env Rscript
err.cat <- function(x)  cat(x, '\n', file=stderr())

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tools))
suppressPackageStartupMessages(library(xtable)) 
suppressPackageStartupMessages(library(snpStats))
suppressPackageStartupMessages(library(GenomicRanges))

### Series of filters suggested by Adam.
# Filtering of variants based on annotation

option_list <- list(
    make_option(c('--root'), help='path to root'),
    make_option(c('--chromosome'),  help='')
) 
option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

print(root <- opt$root)
print(chromosome <- opt$chromosome)

#######
sample.names.file <- paste0(root, '_snpStats/colname_chr', chromosome, '.tab')
sample.names <- scan(sample.names.file, what = character())

####### Now we read the data one chromosome at a time
matrix.calls.snpStats.all <- NULL
first <- TRUE
message('Preparing data for chromosome ', chromosome)
calls.file <- paste(root, '_snpStats/calls_chr', chromosome, '.tab', sep = '')
depth.file <- paste(root, '_snpStats/depth_chr', chromosome, '.tab', sep = '')
rowname.file <- paste(root, '_snpStats/rowname_chr', chromosome, '.tab', sep = '')
annotations.file <- paste(root, '_snpStats/annotations_chr', chromosome, '.csv', sep = '')

if (file.exists(rowname.file) && file.exists(calls.file)) {
  variant.names <- scan(rowname.file, what = character())
  n.calls <- length(variant.names)
  
  annotations.snpStats <- as.data.frame(fread(annotations.file, stringsAsFactors = FALSE))
  annotations.snpStats$signature <- paste(annotations.snpStats$Chr, annotations.snpStats$Start, annotations.snpStats$Ref, annotations.snpStats$Obs, sep = '_')
  annotations.snpStats$clean.signature <- variant.names
  row.names(annotations.snpStats) <- annotations.snpStats$clean.signature
  
  matrix.calls <- matrix(scan(calls.file, what = integer()), 
                         dimnames = list( variant.names, sample.names), byrow = TRUE,
                         ncol = length(sample.names),
                         nrow = n.calls)
  
  matrix.calls.snpStats <- new('SnpMatrix', t(matrix.calls)+1)
  
################# And now the depth matrix
  matrix.depth <- matrix(as.raw( pmin(100, scan(depth.file, na.strings = ".", what = integer()))),
                         dimnames = list( variant.names, sample.names), byrow = TRUE,
                         ncol = length(sample.names),
                         nrow = n.calls)
  
  save(list = c('matrix.depth', 'annotations.snpStats', 'matrix.calls.snpStats'), file = paste(root, '_snpStats/chr', chromosome, '_snpStats.RData', sep = ''))
  
  print(table(annotations.snpStats$FILTER))
}
