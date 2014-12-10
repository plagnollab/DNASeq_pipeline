library(DESeq)


getArgs <- function() {
  myargs.list <- strsplit(grep("=",gsub("--","",commandArgs()),value=TRUE),"=")
  myargs <- lapply(myargs.list,function(x) x[2] )
  names(myargs) <- lapply(myargs.list,function(x) x[1])
  return (myargs)
}



########################## read arguments
deseq.compute <- TRUE
extra.plots <- TRUE
keep.dups <- FALSE

#annotation.file <- '/cluster/project8/vyp/vincent/Software/pipeline/RNASeq/bundle/mouse/biomart/biomart_annotations_mouse.tab'
#iFolder <- '/scratch2/vyp-scratch2/Bochukova_RNASeq/processed'
#support.frame <- 'support/Bochukova.tab'
#code <- 'Bochukova'

annotation.file <- '/cluster/project8/vyp/vincent/Software/pipeline/RNASeq/bundle/mouse/biomart/biomart_annotations_mouse.tab'
iFolder <- '/scratch2/vyp-scratch2/IoN_RNASeq/Natalia/processed'
support.frame <- 'support/Natalia.tab'
code <- 'Natalia'





myArgs <- getArgs()
if ('support.frame' %in% names(myArgs)) support.frame <- myArgs[['support.frame']]
if ('code' %in% names(myArgs)) code <- myArgs[['code']]
if ('iFolder' %in% names(myArgs)) iFolder <- myArgs[['iFolder']]
if ('annotation.file' %in% names(myArgs)) annotation.file <- myArgs[['annotation.file']]
if ('keep.dups' %in% names(myArgs)) keep.dups <- as.logical(myArgs[['keep.dups']])




###check input files and data frame
message('Now reading ', support.frame)
support <- read.table(support.frame, header = TRUE, stringsAsFactors = FALSE)
list.conditions <- grep(names(support), pattern = '^condition.*', value  = TRUE)
annotation <- read.table(annotation.file, header = TRUE, sep = '\t', na.string = c('', 'NA'))

### deseq output folders and files
deseq.folder <- paste(iFolder, '/deseq', sep = '')
for (folder in c(deseq.folder)) {
  if (! file.exists(folder)) dir.create(folder)
}


########## load the count data
if (keep.dups) deseq.counts <- paste(deseq.folder, '/deseq_counts_', code, '_keep_dups.RData', sep = '')  ##this contains the key data
if (!keep.dups) deseq.counts <- paste(deseq.folder, '/deseq_counts_', code, '.RData', sep = '')  ##this contains the key data
load(deseq.counts)



genes.on.XY <- as.character(subset(annotation, chromosome_name %in% c('X' ,'Y'), 'EnsemblID', drop = TRUE))
message('Prior to removing chr XY probes: ', nrow(genes.counts))
genes.counts <- genes.counts[ ! dimnames(genes.counts)[[1]] %in% genes.on.XY, ]                 
message('After removing chr XY probes: ', nrow(genes.counts))
#genes.counts <- genes.counts[1:1000,]

###loop over all proposed conditions
for (condition in list.conditions) {

  genes.counts.loc <- genes.counts[, !is.na(support[, condition]) ]

  support.loc <-  support[  !is.na(support[, condition]), ]
  support.loc$condition <- factor(support.loc[, condition])
  

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
  
################### create the appropriate folders and specify output file
  loc.deseq.folder <- paste(iFolder, '/deseq/', loc.code, sep = '')
  deseq.figs <- paste(loc.deseq.folder, '/figs', sep = '')

  for (folder in c(loc.deseq.folder, deseq.figs)) {
    if (! file.exists(folder)) dir.create(folder)
  }
  
  if (keep.dups) output.file <- paste(loc.deseq.folder, '/deseq_', code, '_differential_expression_keep_dups.tab', sep = '')
  if (!keep.dups) output.file <- paste(loc.deseq.folder, '/deseq_', code, '_differential_expression.tab', sep = '')

######### Now add a PCA for the subset of individuals being considered
  my.sd <- apply(genes.counts.loc, MAR = 1, FUN = sd)
  mat.for.pca <- t(genes.counts.loc[my.sd > median(my.sd), ])
  pca.data <- prcomp(mat.for.pca, scale = TRUE)  

  output.pca <- paste(deseq.figs, '/', loc.code, '_pca.pdf', sep = '')
  pdf(output.pca)

  for (i in 1:2) {
    message('PCA ', i)
    if (length(unique(condition)) <= 4) {
      col <- c('black', 'red', 'green', 'blue')[ as.numeric(factor(support.loc[, condition ])) ]
    } else {
      col <- as.numeric(factor(support.loc[, condition ]))
    }

    
    plot(x = pca.data$x[,2*i-1],
         y = pca.data$x[,2*i],
         xlab = ifelse (i == 1, 'PC1', 'PC3'),
         ylab = ifelse (i == 1, 'PC2', 'PC4'),
         col = col,
         pch = '+')
    
    text(x = pca.data$x[,2*i -1],
         y = pca.data$x[,2*i],
         labels = as.character(support.loc$sample),
         pos = 3)

    my.levels <- levels(factor(support.loc[, condition]))
    legend(col = 1:length(my.levels),
           legend = my.levels,
           pch = '+',
           x = 'bottomright')
     
  }
  print(output.pca)
  dev.off()

  
###################
  method <- 'pooled'
  use.type <- FALSE
  if ('type' %in% names(support.loc)) {
    if ( length(unique(support.loc$type)) > 1) use.type <- TRUE
  }
  
  if (use.type) {
    formula1 <- count ~ type + condition
    formula0 <- count ~ type
    
    support.loc$type <- factor(support.loc$type)
    design.deseq <- support.loc[, c('condition', 'type')]
    CDS <- newCountDataSet(genes.counts.loc, conditions = design.deseq)    

    if ( (max(table(design.deseq$condition)) > 1) && (max(table(design.deseq$condition, design.deseq$type)) == 1) ) {method <- 'pooled-CR'}
    if (max(table(design.deseq$condition)) == 1) {method <- 'blind'}
    
  } else {
    formula1 <- count ~ condition
    formula0 <- count ~ 1
    design.deseq <- support.loc[, c('condition'), drop = FALSE]
    CDS <- newCountDataSet(genes.counts.loc, condition = support.loc$condition)    
  }

  
  ### output the design information
  row.names(design.deseq) <- support.loc$sample
  write.table(x = design.deseq, file = paste(loc.deseq.folder, '/design.tab', sep = ''), row.names = TRUE, quote = FALSE, sep = '\t')
  message('Estimation method is ', method, ' for ', loc.code)
  
#################### compute the P-values
  CDS <- estimateSizeFactors(CDS)
  CDS <- estimateDispersions(CDS, method = method)

####
  output.pdf <- paste(deseq.figs, '/', code, '_', loc.code, '_disp_estimates.pdf', sep = '')
  pdf(output.pdf)
  plotDispEsts( CDS )
  dev.off()
  print(output.pdf)
####
  fit0 <- fitNbinomGLMs( CDS, formula0 )
  fit1 <- fitNbinomGLMs( CDS, formula1 )

################################
  deseq.pval <- fit1
  deseq.pval$EnsemblID <- row.names( deseq.pval)

  ###################### add the read count by condition
  for (cond in levels(support.loc$condition)) {
    deseq.pval[, paste('average.count.', cond, sep = '')] <- as.numeric(apply(MAR =1, genes.counts.loc[, support.loc$condition == cond, drop = FALSE], FUN = mean))
  }

  #################
  deseq.pval$basic.pval <- signif(nbinomGLMTest( fit1, fit0 ), 4)
  deseq.pval$adj.pval <- signif(p.adjust( deseq.pval$basic.pval, method="BH" ), 4)
  deseq.pval <- merge(deseq.pval, annotation, by = 'EnsemblID', all.x = TRUE)
  deseq.pval <- deseq.pval[ order(deseq.pval$basic.pval),]


  ######### add a zscore
  name.conditions <- grep (names(deseq.pval), pattern = '^condition.*', value = TRUE)  ## need that for the z-score
  name.counts <- grep (names(deseq.pval), pattern = '^average.count.*', value = TRUE)

  if (length(name.conditions) == 1) { ## I think this is not required but just in case
    basic.pval.floor <- ifelse (  deseq.pval$basic.pval < 10^-15, 10^-15, deseq.pval$basic.pval )
    deseq.pval$zscore <- qnorm(p = 1 - basic.pval.floor/2)
    deseq.pval$zscore<- ifelse ( deseq.pval[, name.conditions] > 0, abs(deseq.pval$zscore), - abs(deseq.pval$zscore))
  } else {
    message('Error in z-score computation')
    print(name.conditions)
  }
  
#################### now print things
  my.labels <- c('EnsemblID', 'external_gene_id', 'basic.pval', 'adj.pval', 'zscore', name.conditions, name.counts, 'chromosome_name', 'start_position', 'end_position', 'strand')
  write.table(x = deseq.pval[, subset(my.labels, my.labels %in% names(deseq.pval))], file = output.file, row.names = FALSE, quote = FALSE, sep = '\t')
  print(output.file)


#### more graphs


#### heatmap to start
  if (extra.plots) {
    if (length(unique(support.loc[, condition])) <= 4) {
      genes.diff <- deseq.pval$EnsemblID[1:100]
      genes.counts.loc.heatmap <- genes.counts.loc[ as.character(genes.diff), ]
      for (col in 1:ncol(genes.counts.loc.heatmap)) genes.counts.loc.heatmap[,col] <- genes.counts.loc.heatmap[,col]/sum(genes.counts.loc.heatmap[,col])*1000000

      if (!keep.dups) {
        dimnames(genes.counts.loc.heatmap)[[2]] <- gsub(pattern = '_dexseq_counts.txt', replacement = '', basename(  dimnames(genes.counts.loc.heatmap)[[2]] ))
        output.heatmap <- paste(loc.deseq.folder, '/deseq_', code, '_heatmap.pdf', sep = '')
      }
      if ( keep.dups) {
        dimnames(genes.counts.loc.heatmap)[[2]] <- gsub(pattern = '_dexseq_counts_keep_dups.txt', replacement = '', basename(  dimnames(genes.counts.loc.heatmap)[[2]] ))
        output.heatmap <- paste(loc.deseq.folder, '/deseq_', code, '_heatmap_keep_dups.pdf', sep = '')
      }

      dimnames(genes.counts.loc.heatmap)[[1]] <- as.character(deseq.pval$external_gene_id[ match(dimnames(genes.counts.loc.heatmap)[[1]], table = deseq.pval$EnsemblID) ])
      
      pdf(output.heatmap)
      heatmap(x = genes.counts.loc.heatmap,
              ColSideColors = c('black', 'red', 'blue', 'green')[as.numeric(factor(support.loc[, condition]))],
              margins = c(10, 10),
              cexRow = 0.5)
      dev.off()
      print(output.heatmap)
    }
  }
########
  
  if ( length(levels(factor(support.loc[, condition]))) == 2) {
    output.pdf <- paste(deseq.figs, '/', code, '_', loc.code, '_MAplot.pdf', sep = '')
    pdf(output.pdf)
    for.plot <- data.frame(baseMean = deseq.pval[[ '(Intercept)']],
                           log2FoldChange = deseq.pval[, name.conditions],
                           padj = deseq.pval$adj.pval)
    plotMA(for.plot, ylim = c(-3, 3))
    dev.off()
    print(output.pdf)
  }
}


warnings()


sessionInfo()
