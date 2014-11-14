library(DEXSeq)
library(parallel)

#source("/cluster/project2/vyp/vincent/Software/pipeline/RNASeq/dexseq_tools.R")

getArgs <- function() {
  myargs.list <- strsplit(grep("=",gsub("--","",commandArgs()),value=TRUE),"=")
  myargs <- lapply(myargs.list,function(x) x[2] )
  names(myargs) <- lapply(myargs.list,function(x) x[1])
  return (myargs)
}

plotDispEsts <- function( cds, ymin, linecol="#ff000080", xlab = "mean of normalized counts", ylab = "dispersion", log = "xy", cex = 0.45, ... ) {
      px = rowMeans( counts( cds, normalized=TRUE ) )
      sel = (px>0)
      px = px[sel]
      
      py = fData(cds)$dispBeforeSharing[sel]
      if(missing(ymin)) ymin = 10^floor(log10(min(py[py>0], na.rm=TRUE))-0.1)
      
      plot(px, pmax(py, ymin), xlab=xlab, ylab=ylab,
           log=log, pch=ifelse(py<ymin, 6, 16), cex=cex, ... )
      xg = 10^seq( -.5, 5, length.out=100 )
      fun = function(x) { cds@dispFitCoefs[1] + cds@dispFitCoefs[2] / x }
      lines( xg, fun(xg), col=linecol, lwd=4)
    }

plotMA <- function(x, ylim,col = ifelse(x$padj>=0.1, "gray32", "red3"), linecol = "#ff000080", xlab = "mean of normalized counts", ylab = expression(log[2]~fold~change), log = "x", cex=0.45, ...) {
  if (!(is.data.frame(x) && all(c("baseMean", "log2FoldChange") %in% colnames(x))))
    stop("'x' must be a data frame with columns named 'baseMean', 'log2FoldChange'.")
  
  x = subset(x, baseMean!=0)
  py = x$log2FoldChange
  if(missing(ylim))
    ylim = c(-1,1) * quantile(abs(py[is.finite(py)]), probs=0.99) * 1.1
  plot(x$baseMean, pmax(ylim[1], pmin(ylim[2], py)),
       log=log, pch=ifelse(py<ylim[1], 6, ifelse(py>ylim[2], 2, 16)),
       cex=cex, col=col, xlab=xlab, ylab=ylab, ylim=ylim, ...)
  abline(h=0, lwd=4, col=linecol)
}

fitDispersionFunction <- function (ecs) {  ##function that includes a fix by VP
  stopifnot(is(ecs, "ExonCountSet"))
  if (all(is.na(fData(ecs)$dispBeforeSharing))) {
    stop("no CR dispersion estimations found, please first call estimateDispersions function")
  }
  means <- colMeans(t(counts(ecs))/sizeFactors(ecs))
  disps <- fData(ecs)$dispBeforeSharing
  coefs <- c(0.1, 1)
  iter <- 0
  while (TRUE) {
    residuals <- disps/(coefs[1] + coefs[2]/means)
    good <- which((residuals > 1e-04) & (residuals < 15) )
    mm <- model.matrix(disps[good] ~ I(1/means[good]))
    fit <- try(statmod:::glmgam.fit(mm, disps[good], coef.start = coefs),
               silent = TRUE)
    if (inherits(fit, "try-error")) {
      stop("Failed to fit the dispersion function\n")
    }
    oldcoefs <- coefs
    coefs <- coefficients(fit)
    #message(coefs)
    if (coefs[1] < 0) {
      coefs[1] <- 0.001
      warning("Negative intercept value in the dispersion function, it will be set to 0. Check fit diagnostics plot section from the vignette.")
      #print(coefs)
      #break
    }
    if (sum(log(coefs/oldcoefs)^2) < 0.005)
      break
    iter <- iter + 1
    if (iter > 10) {
      warning("Dispersion fit did not converge.")
      break
    }
  }
  ecs@dispFitCoefs <- coefs
  fData(ecs)$dispFitted <- ecs@dispFitCoefs[1] + ecs@dispFitCoefs[2]/colMeans(t(counts(ecs))/sizeFactors(ecs))
  fData(ecs)$dispersion <- pmin(pmax(fData(ecs)$dispBeforeSharing,
                                     fData(ecs)$dispFitted, na.rm = TRUE), 1e+08)
  return(ecs)
}


getArgs <- function() {
  myargs.list <- strsplit(grep("=",gsub("--","",commandArgs()),value=TRUE),"=")
  myargs <- lapply(myargs.list,function(x) x[2] )
  names(myargs) <- lapply(myargs.list,function(x) x[1])
  return (myargs)
}



########################## read arguments
dexseq.compute <- TRUE

annotation.file <- '/cluster/project8/vyp/vincent/Software/pipeline/RNASeq/bundle/mouse/biomart/biomart_annotations_mouse.tab'
iFolder <- '/scratch3/vyp-scratch2/Bochukova_RNASeq/processed'
support.frame <- 'support/Bochukova.tab'
code <- 'Bochukova'


myArgs <- getArgs()
if ('support.frame' %in% names(myArgs)) support.frame <- myArgs[['support.frame']]
if ('code' %in% names(myArgs)) code <- myArgs[['code']]
if ('iFolder' %in% names(myArgs)) iFolder <- myArgs[['iFolder']]
if ('annotation.file' %in% names(myArgs)) annotation.file <- myArgs[['annotation.file']]





###check input files and data frame
message('Now reading ', support.frame)
support <- read.table(support.frame, header = TRUE, stringsAsFactors = FALSE)
my.ids <- support$sample
list.conditions <- grep(names(support), pattern = '^condition.*', value  = TRUE)
annotation <- read.table(annotation.file, header = TRUE, sep = '\t', na.string = c('NA', ''))

files <- paste(iFolder, '/', my.ids, '/dexseq/', my.ids, '_dexseq_counts.txt', sep = '')
if (sum(!file.exists(files)) > 0) stop('Some input files are missing')




### dexseq output folders
dexseq.folder <- paste(iFolder, '/dexseq', sep = '')
dexseq.counts <- paste(dexseq.folder, '/dexseq_counts_', code, '.RData', sep = '')  ##this contains the key data
if (!file.exists(dexseq.folder)) dir.create(dexseq.folder)



for (condition in list.conditions) {

  message('Condition ', condition)
  support.loc <- support

  ##handle the type variable
  support.loc$condition <- factor(support[, condition])
  support.loc <-  subset(support.loc, !is.na(condition) )

  ##handle the type variable
  type.loc <- gsub(x = condition, pattern = 'condition', replacement = 'type')
  if ( (! type.loc %in% names(support.loc)) & ('type' %in% names(support.loc))) {type.loc <- 'type'}  ##if only "type" is present, use it
  if (type.loc %in% names(support)) {
    support.loc$type <- factor(support.loc[, type.loc])
    support.loc <- subset(support.loc, !is.na(type))
  }
  
  loc.code <-  paste(unique(support.loc$condition), collapse = '_')
  message('Support data frame', loc.code)
  print(support.loc)

  
  ################### create the appropriate folders
  loc.dexseq.folder <- paste(iFolder, '/dexseq/', loc.code, sep = '')
  dexseq.figs <- paste(loc.dexseq.folder, '/figs', sep = '')
  dexseq.data <- paste(loc.dexseq.folder, '/dexseq_', code, '_', loc.code, '.RData', sep = '')  ##this will contain the final output of dexseq
  
  for (folder in c(loc.dexseq.folder, dexseq.figs)) {
    if (! file.exists(folder)) dir.create(folder)
  }

  
  if (dexseq.compute) {  
    load(dexseq.counts)  ##object mycounts is key
    #DexSeqExons <- subset(DexSeqExons, c(rep(TRUE, 300), rep(FALSE, nrow(counts(DexSeqExons)) - 300))) 
    

    use.covariate <- FALSE
    if ('type' %in% names(support.loc)) {
      if (length(unique(as.character(support.loc$type))) > 1) {
        use.covariate <- TRUE
      }
    }

    if (use.covariate) {
      formuladispersion <- ~ sample + (condition + type) * exon
      formula0 <-  ~ sample + type * exon + condition
      formula1 <-  ~ sample + type * exon + condition * exon
      my.design <- support.loc[, c('type', 'condition')]
      my.design$type <- factor(my.design$type) ## probably not needed
      my.design$condition <- factor(my.design$condition)  ## probably not needed
      my.design.loc <- my.design  ##just to print basically
    } else {
      formuladispersion <-  ~ sample + condition * exon
      formula0 <-  ~ sample + exon + condition
      formula1 <-  ~ sample + exon + condition * exon
      my.design <- factor(support.loc[, c('condition')])
      my.design.loc <- support.loc[, c('condition'), drop = FALSE]  ##just to print basically
    }

    row.names(my.design.loc) <- support.loc$sample
    write.table(x = my.design.loc, file = paste(loc.dexseq.folder, '/design.tab', sep = ''), row.names = TRUE, quote = FALSE, sep = '\t')


    message('Updating the dexseq object')
    DexSeqExons.loc <- newExonCountSet(countData=counts(DexSeqExons)[, !is.na(support[, condition])],
                                       design= my.design,
                                       geneIDs=geneIDs(DexSeqExons),
                                       exonIDs=exonIDs(DexSeqExons),
                                       transcripts = DexSeqExons@featureData@data$transcripts,
                                       exonIntervals = DexSeqExons@featureData@data[, c('chr', 'start', 'end', 'strand')])
    
    message('Starting the computations')
    DexSeqExons.loc <- estimateSizeFactors(DexSeqExons.loc)

    
    message('Here is the part that takes a lot of time')
    DexSeqExons.loc <- DEXSeq::estimateDispersions(object=DexSeqExons.loc,
                                                    nCores=8,
                                                    minCount=20,
                                                    formula = as.formula(formuladispersion))
    message('Done with estimateDispersions')
    
    DexSeqExons.loc <- fitDispersionFunction(DexSeqExons.loc)
    #DexSeqExons.loc <- testForDEU(DexSeqExons.loc, nCores=8, formula0 = as.formula(formula0), formula1 = as.formula(formula1))
    DexSeqExons.loc <- testForDEU(DexSeqExons.loc, nCores=8)
    DexSeqExons.loc <- estimatelog2FoldChanges(DexSeqExons.loc)

######################### output basic table
    res <- DEUresultTable(DexSeqExons.loc)
    
    names(res)[1:6] <- c("EnsemblID", "exonID", "dispersion", "pvalue", "padjust", "meanBase")
    res <- merge(res, annotation[, c('EnsemblID', 'external_gene_id', 'chromosome_name', 'start_position', 'end_position', 'strand')], by = 'EnsemblID')
    res <- res[ order(res$pvalue),]  
    
    write.csv(x = res,
              file=paste(loc.dexseq.folder, "/", code, "_", loc.code, "_SignificantExons.csv", sep = ''),
              row.names = FALSE)
    
    message('Saving results in ', dexseq.data)
    save(list = c('res', 'DexSeqExons.loc'), file = dexseq.data)
  } else {
    load(dexseq.data)
  }


########################## Now plot a subset
  file.remove(list.files(dexseq.figs, pattern = 'DEXSeq*', full.names = TRUE)) ##remove the old plots

  n.sig <- sum(res$padjust < 0.01, na.rm = TRUE)
  if (n.sig > 50) {
    resSigs <- subset(res, padjust<0.01)
  } else resSigs <- res[1:50,]


  for (i in 1:nrow(resSigs)) {
    gene <- as.character(resSigs$EnsemblID[ i ])
    gene.pretty <- as.character(resSigs$external_gene_id[ i ])
    
    message(i, ' ', gene, ' ', gene.pretty)
    output.pdf <- paste(dexseq.figs, '/DEXSeq-', gene.pretty, '.pdf', sep = '')
    pdf(output.pdf, width = 8, height = 4.9)
    try(plotDEXSeq(DexSeqExons.loc,
                   geneID = gene,  ##I suspect it has to be gene, otherwise it crashes
                   cex.axis = 1.2,
                   cex=1.3,
                   lwd=2,
                   legend=TRUE,
                   displayTranscripts = TRUE,
                   names = TRUE,
                   main = gene.pretty)
        )
    #stop()                                   #color=c("red", "blue", "darkgreen"))
    dev.off()
    print(output.pdf)
  }
  
  
#################################### plot some graphs
  
  pdf(file=paste(dexseq.figs, '/DEXSeq-MeanVsDispPoints.pdf', sep = ''))
  plotDispEsts (DexSeqExons.loc)
  dev.off()
  
    
  pdf(file=paste(dexseq.figs, '/DEXSeq-MeanVsDispCircles.png', sep = ''))
  try(plotMA(with(res,
                  data.frame(baseMean = res[,6],log2FoldChange = res[,7], padj = res[,5]),
                  ylim=c(-4,4), cex=0.8)))
  dev.off()

  message('Done with ', condition)
  rm('DexSeqExons.loc')
  gc()

}

warnings()


#head(modelFrameForGene(DexSeqExons, "Tardbp"))
#testGeneForDEU(DexSeqExons, "Tardbp")


sessionInfo()
