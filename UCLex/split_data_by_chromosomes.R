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


#######

sample.names.file <- paste(root, '_colname.tab', sep = '')
sample.names <- scan(sample.names.file, what = character())


onekgpositions <- read.table('data/filteredPurcell_final.012.pos', header = FALSE, col.names = c('CHROM', 'POS'))
OneKG_ranges <- GRanges(seqnames = Rle(onekgpositions$CHROM), IRanges(onekgpositions$POS, onekgpositions$POS))

####### Now we read the data one chromosome at a time
matrix.calls.snpStats.all <- NULL
first <- TRUE

for (chrom in c(as.character(1:22), 'X')) {
#for (chrom in c(as.character(c(1)))) {
  message('Preparing data for chromosome ', chrom)

  calls.file <- paste(root, '_by_chr/calls_chr', chrom, '.tab', sep = '')
  depth.file <- paste(root, '_by_chr/depth_chr', chrom, '.tab', sep = '')
  rowname.file <- paste(root, '_by_chr/rowname_chr', chrom, '.tab', sep = '')
  annotations.file <- paste(root, '_by_chr/annotations_chr', chrom, '.csv', sep = '')

  if (file.exists(rowname.file) && file.exists(calls.file)) {
    variant.names <- scan(rowname.file, what = character())
    n.calls <- length(variant.names)
    
    annotations.snpStats <- read.csv(annotations.file, stringsAsFactors = FALSE)
    annotations.snpStats$signature <- paste(annotations.snpStats$Chr, annotations.snpStats$Start, annotations.snpStats$Ref, annotations.snpStats$Obs, sep = '_')
    annotations.snpStats$clean.signature <- variant.names
    row.names(annotations.snpStats) <- annotations.snpStats$clean.signature

#### this is for the PCA work, just keep common variants
    annotations_ranges <- GRanges(seqnames = Rle(annotations.snpStats$Chr), IRanges(annotations.snpStats$Start, annotations.snpStats$End))
    overlap <- data.frame(as.matrix(findOverlaps(annotations_ranges, OneKG_ranges)))
    
    matrix.calls <- matrix(scan(calls.file, what = integer()), 
                           dimnames = list( variant.names, sample.names), byrow = TRUE,
                           ncol = length(sample.names),
                           nrow = n.calls)
    
    matrix.calls.snpStats <- new('SnpMatrix', t(matrix.calls)+1)
    #print(table(as(matrix.calls.snpStats[, '1_897325_G_C'], 'numeric')[,1])); stop()
    
################# And now the depth matrix
    matrix.depth <- matrix(as.raw( pmin(100, scan(depth.file, na.strings = ".", what = integer()))),
                           dimnames = list( variant.names, sample.names), byrow = TRUE,
                           ncol = length(sample.names),
                           nrow = n.calls)
    
    save(list = c('matrix.depth', 'annotations.snpStats', 'matrix.calls.snpStats'), file = paste(root, '_by_chr/chr', chrom, '_snpStats.RData', sep = ''))
    
    print(table(annotations.snpStats$FILTER))
    good.pos <- unique(overlap$queryHits)
    
    if (chrom != 'X') { ### no SNP on chrom X so don't merge in this case
      if (first) {
        matrix.calls.snpStats.all <- matrix.calls.snpStats[, good.pos ]
        annotations <- annotations.snpStats[ good.pos,]
        first <- FALSE
      } else {
        if (length(good.pos) > 0) {
          matrix.calls.snpStats.all <- cbind( matrix.calls.snpStats.all, matrix.calls.snpStats[, good.pos])
          annotations <- rbind.data.frame(annotations, annotations.snpStats[ good.pos, ])
        }
      }
    }
    rm(matrix.calls)
    gc()
  }
}

annotations <- annotations[ dimnames(matrix.calls.snpStats.all)[[2]],  ]
save(list = c('annotations', 'matrix.calls.snpStats.all'), file = paste(root, '_forPCA_calls_snpStats.RData', sep = ''))



### correct an odd bug
#my.tab <- table(annotations$signature)
#bad.rows <- data$signature %in% names(which(my.tab == 2))
#data <- data[ ! bad.rows, ]

