library(snpStats)
library(GenomicRanges)

#tstamp <- 'April2014'  ##813 cases, 436889 variants
#tstamp <- Sys.getenv("currentUCLex")
#if (tstamp == "") tstamp <- 'May2014'

getArgs <- function() {
  myargs.list <- strsplit(grep("=",gsub("--","",commandArgs()),value=TRUE),"=")
  myargs <- lapply(myargs.list,function(x) x[2] )
  names(myargs) <- lapply(myargs.list,function(x) x[1])
  return (myargs)
}


myArgs <- getArgs()


root <- '/scratch2/vyp-scratch2/vincent/GATK/mainset_October2014/mainset_October2014'
if ('root' %in% names(myArgs)) root <- myArgs[[ 'root' ]]
if ('chromosome' %in% names(myArgs)) chromosome <- as.character(myArgs[[ 'chromosome' ]])


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
  
  annotations.snpStats <- read.csv(annotations.file, stringsAsFactors = FALSE)
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
