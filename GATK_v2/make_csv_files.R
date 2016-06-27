##create a snpStats object with annotations (including external control set) that has only the variant found in at least one case

###TODO: also save individual per chromosome files with case, controls and basic annotations-> useful for SKAT test
###TODO: create equivalent files with read depth information

library(snpStats)
source('process_multiVCF.R')
options(stringsAsFactors = FALSE)

getArgs <- function() {
  myargs.list <- strsplit(grep("=",gsub("--","",commandArgs()),value=TRUE),"=")
  myargs <- lapply(myargs.list,function(x) x[2] )
  names(myargs) <- lapply(myargs.list,function(x) x[1])
  return (myargs)
}


case.control.analysis <- function(choice.cases = NULL, output = 'processed/support/case_control', known.genes = c(), SKAT = FALSE, fix.names = NULL, annotations.file = '/cluster/project8/vyp/vincent/Software/RNASeq_pipeline/bundle/human/biomart/biomart_extra_annotations_human.tab') {
  message('Now running the case control script') 
  options(stringsAsFactors = FALSE) 
  library(snpStats)
  if (SKAT) library(SKAT) 
  final.frame <- NULL
  first <- TRUE 
  output.file <- paste(output, '_gene_based_summary.csv', sep = '')
  output.single <- paste(output, '_single_variant_summary.csv', sep = '') 
  gene.annotations <- read.table(file = annotations.file, stringsAsFactors = FALSE, na.string = c('', 'NA'), sep = '\t', header = TRUE, comment.char = '|', quote= '') 
  final.table <- data.frame()
  final.single <- data.frame() 
  for (chr in as.character(1:22)) { 
    input.file <- paste('data/case_control_chr', chr, '.RData', sep = '') 
    if (file.exists(input.file)) {
      load(input.file) 
      if (!is.null(fix.names)) {
        message('Fixing the names first in the case control module for chromosome ', chr)
        if (! 'New_correct_name' %in% names(fix.names)) stop('Missing New_correct_name column in the fix name table')
        if (! 'Old_name' %in% names(fix.names)) stop('Missing Old_name column in the fix name table') 
        fix.names$New_correct_name <- paste(fix.names$New_correct_name, 'z', sep = '')
        my.names <- dimnames(case.control.snpStats)[[1]]
        for (i in 1:nrow(fix.names)) {my.names <- ifelse ( my.names == fix.names$Old_name[ i ], fix.names$New_correct_name[ i ], my.names)}
        my.names <- gsub(my.names, pattern = 'z$', replacement = '')
        dimnames(case.control.snpStats)[[1]] <- my.names 
        message('Now done fixing the names')
      }
##########
      #annotations$my.pass <- (annotations$FILTER == 'PASS' | grepl(pattern = 'INDEL', annotations$FILTER))  ###important filter here, I remove the FAIL filters
      annotations$my.pass <- TRUE 
      annotations$cat1 <- annotations$somewhat.rare & annotations$my.pass & (annotations$lof | annotations$non.syn | annotations$splicing)  
      annotations$cat2 <- annotations$rare & annotations$my.pass & annotations$lof 
##################### find the case control labels
      if (first) {
        my.controls <- dimnames(case.control.snpStats)[[1]][ case.control.labels == 0 ]
        my.cases <- dimnames(case.control.snpStats)[[1]] [ case.control.labels == 1 ]
        if (!is.null(choice.cases)) my.cases <- subset(my.cases, my.cases %in% choice.cases)
        first <- FALSE
      } 
###### single variant analysis
      case.control.snpStats <- case.control.snpStats[c(my.cases, my.controls), ]
      my.data <- data.frame ( cc = c(rep(1, length(my.cases)), rep(0, length(my.controls))))
      row.names(my.data) <- dimnames(case.control.snpStats)[[1]] 
      my.pvals <- snp.rhs.tests ( snp.data = case.control.snpStats, formula = as.formula('cc ~ 1'), data = my.data, family = "binomial")
      annotations$single.variant.p.value <- p.value(my.pvals)
      interesting.variants <- subset(annotations, single.variant.p.value < 0.01) 
      for.fisher <- as.matrix(interesting.variants[, c('non.missing.controls', 'non.ref.calls.controls', 'non.missing.cases', 'non.ref.calls.cases') ])
      interesting.variants$fisher.pvalue <- as.numeric(apply(for.fisher, MAR = 1, FUN = function(x) {fisher.test( matrix( data = c(2*x[1] - x[2], x[2], 2*x[3] - x[4], x[4]), nrow = 2, ncol = 2))$p.value}))
      interesting.variants <- subset(interesting.variants, fisher.pvalue < 0.01) 
      my.names <- c('fisher.pvalue', 'dbSNP137', 'HUGO', 'ExonicFunc', 'Description', 'AAChange', 'clean.signature', 'freq.controls', 'freq.cases', 'rare', 'somewhat.rare', 'lof', 'non.syn', 'splicing', 'ESP6500si_ALL', 'X1000g2012apr_ALL', 'Chr', 'Start', 'End', 'Ref', 'Obs', 'FILTER', 'non.missing.controls', 'non.ref.calls.controls', 'non.missing.cases', 'non.ref.calls.cases')
      interesting.variants <- interesting.variants[, my.names ]
      final.single <- rbind.data.frame( final.single, interesting.variants ) 
########### Now the gene based tests
      gene.list <- unique(annotations$ensemblID)
      chr.table <- data.frame (ensemblID = gene.list,
                               HUGO = NA,
                               chr = chr,
                               start.gene = NA,
                               end.gene = NA,
                               cases.biallelic.n.somewhat.rare.func = NA,
                               controls.biallelic.n.somewhat.rare.func = NA,
                               cases.n.somewhat.rare.func = NA,
                               controls.n.somewhat.rare.func = NA,
                               cases.n.hom.somewhat.rare.func = NA,
                               controls.n.hom.somewhat.rare.func = NA,
                               cases.n.rare.lof = NA,
                               controls.n.rare.lof = NA,
                               cases.n.hom.rare.lof = NA,
                               controls.n.hom.rare.lof = NA,
                               case.control.pval.biallelic.func = NA,
                               case.control.pval.func = NA,
                               case.control.pval.lof = NA,
                               sample.size.cases = NA,
                               sample.size.controls = NA) 
      chr.table$HUGO <- gene.annotations$external_gene_id[ match(chr.table$ensemblID, table = gene.annotations$EnsemblID) ]
      chr.table$start.gene <- gene.annotations$start_position [ match(chr.table$ensemblID, table = gene.annotations$EnsemblID) ]
      chr.table$end.gene <- gene.annotations$end_position [ match(chr.table$ensemblID, table = gene.annotations$EnsemblID) ]
      if (length(known.genes) > 0) {chr.table$known.gene <- chr.table$HUGO %in% known.genes} 
################# prepare the case control data in a nice friendly format
      control.data <- as(case.control.snpStats[ my.controls, , drop = FALSE], 'numeric')
      case.data <- as(case.control.snpStats[ my.cases, , drop = FALSE], 'numeric')
      if (SKAT) {
        chr.table$SKAT.pvalue <- NA
        binary.phenotype <- c(rep(1, times = length(my.cases)), rep(0, times = length(my.controls)) )
        obj <- SKAT_Null_Model(as.numeric(binary.phenotype) ~ 1, out_type="D") 
        case.control.data <- as(case.control.snpStats, 'numeric')
      } 
############# Now loop over the genes
      n.genes <- length(gene.list)
      for (i in 1:n.genes) {
        gene <- gene.list[ i ]
        symbol <- chr.table$HUGO[ i ]
        if (i %% 100 == 0) {message('Looking at gene ', symbol, ' on chromosome ', chr, ', nb ', i, ' out of ', n.genes)} 
######### Determine the sample size first
        gene.variants <- which( annotations$ensemblID == gene & annotations$my.pass )
        if (length(gene.variants) > 0) {
          case.data.loc <- case.data[, gene.variants, drop = FALSE]
          control.data.loc <- control.data[, gene.variants, drop = FALSE]
          sample.size.controls <- round(sum(!is.na(control.data.loc)) /ncol(control.data.loc))
          sample.size.cases <- round(sum(!is.na(case.data.loc)) /ncol(case.data.loc))
          ##if ( sample.size.cases > 30) browser() 
          chr.table$sample.size.controls[ i ] <- sample.size.controls
          chr.table$sample.size.cases[ i ] <- sample.size.cases  ##VP 
        ######### Now the SKAT test
          if (SKAT) {
            case.control.data.loc <- case.control.data[, gene.variants, drop = FALSE]
            chr.table$SKAT.pvalue[ i ]<- SKAT(case.control.data.loc, obj)$p.value
          } 
          proba.mut <- sample.size.cases  / ( sample.size.cases +  sample.size.controls)
        } else {proba.mut <- NA} 
########## now the non syn, splicing and LOF variants first
        good.variants <- which (annotations$ensemblID == gene & annotations$cat1) 
        if (length(good.variants) == 0) {
          chr.table$cases.n.somewhat.rare.func[ i ] <- 0
          chr.table$controls.n.somewhat.rare.func[ i ] <- 0
          chr.table$case.control.pval.func[ i ] <- 1
          chr.table$cases.n.hom.somewhat.rare.func[ i ] <- 0
          chr.table$controls.n.hom.somewhat.rare.func[ i ] <- 0
        } else {
          case.data.loc <- case.data[, good.variants, drop = FALSE]
          control.data.loc <- control.data[, good.variants, drop = FALSE] 
          chr.table$cases.biallelic.n.somewhat.rare.func[ i ] <-  sum(apply(case.data.loc, MAR = 1, FUN = sum, na.rm = TRUE) >= 2)
          chr.table$controls.biallelic.n.somewhat.rare.func[ i ] <-  sum(apply(control.data.loc, MAR = 1, FUN = sum, na.rm = TRUE) >= 2) 
          chr.table$cases.n.somewhat.rare.func[ i ] <- sum(case.data.loc, na.rm = TRUE)
          chr.table$controls.n.somewhat.rare.func[ i ] <- sum(control.data.loc, na.rm = TRUE) 
          chr.table$cases.n.hom.somewhat.rare.func[ i ] <- sum(case.data.loc == 2, na.rm = TRUE)
          chr.table$controls.n.hom.somewhat.rare.func[ i ] <- sum(control.data.loc == 2, na.rm = TRUE) 
          if ( (proba.mut > 1/1000) && !is.na(proba.mut)) {
            chr.table$case.control.pval.func[ i ] <- pbinom(q = chr.table$cases.n.somewhat.rare.func[ i ] - 1,
                                                            size = chr.table$cases.n.somewhat.rare.func[ i ] + chr.table$controls.n.somewhat.rare.func[ i ],
                                                            prob = proba.mut,
                                                            lower.tail = FALSE) 
            chr.table$case.control.pval.biallelic.func[ i ] <- pbinom(q = chr.table$cases.biallelic.n.somewhat.rare.func[ i ] - 1, size = chr.table$cases.biallelic.n.somewhat.rare.func[ i ] + chr.table$controls.biallelic.n.somewhat.rare.func[ i ], prob = proba.mut, lower.tail = FALSE) 
          } else  {chr.table$case.control.pval.func[ i ] <- 1} 
        } 
################
        good.variants <- which (annotations$ensemblID == gene & annotations$cat2) 
        if (length(good.variants) == 0) {
          chr.table$cases.n.rare.lof[ i ] <- 0
          chr.table$controls.n.rare.lof[ i ] <- 0
          chr.table$case.control.pval.lof[ i ] <- 1
          chr.table$cases.n.hom.rare.lof[ i ] <- 0
          chr.table$controls.n.hom.rare.lof[ i ] <- 0
        } else {
          case.data.loc <- case.data[, good.variants, drop = FALSE]
          control.data.loc <- control.data[, good.variants, drop = FALSE] 
          chr.table$cases.n.rare.lof[ i ] <- sum(case.data.loc, na.rm = TRUE)
          chr.table$controls.n.rare.lof[ i ] <- sum(control.data.loc, na.rm = TRUE) 
          chr.table$cases.n.hom.rare.lof[ i ] <- sum(case.data.loc == 2, na.rm = TRUE)
          chr.table$controls.n.hom.rare.lof[ i ] <- sum(control.data.loc == 2, na.rm = TRUE)
        } 
        if ( (proba.mut > 1/100) && !is.na(proba.mut) )  {
          chr.table$case.control.pval.lof[ i ] <- pbinom(q = chr.table$cases.n.rare.lof[ i ] - 1, size = chr.table$cases.n.rare.lof[ i ] + chr.table$controls.n.rare.lof[ i ], prob = proba.mut, lower.tail = FALSE)
        } else {
          chr.table$case.control.pval.lof[ i ] <- 1
        }
      } 
      print(subset(chr.table, case.control.pval.func < 0.01 | case.control.pval.lof < 0.01))
      final.table <- rbind.data.frame ( final.table, chr.table) 
      #write.csv(x = final.table[ order(pmin(final.table$case.control.pval.lof, final.table$case.control.pval.func), decreasing = FALSE), ],
      #          file = output.file, row.names = FALSE)
      ### order the gene based list by recessive P-values
      write.csv(x = final.table[ order(final.table$case.control.pval.biallelic.func, decreasing = FALSE), ], file = output.file, row.names = FALSE) 
      write.csv( x = final.single[ order(final.single$fisher.pvalue, decreasing = FALSE), ], file = output.single, row.names = FALSE) 
    }
  } 
}




annotate.standard.annovar.output <- function(data,   ##this function does NOT include the control data to define the rare/somewhat.rare flags
                                             threshold.rare = 0.002,
                                             threshold.somewhat.rare = 0.005,
                                             bad.genes = c(),
                                             freq.fields = c( 'X1000g2012apr_ALL', 'ESP6500si_ALL'),
                                             choice.transcripts = NULL,
                                             biomart.gene.description = '/cluster/project8/vyp/vincent/Software/RNASeq_pipeline/bundle/human/biomart/biomart_extra_annotations_human.tab') {

  data$Chr <- gsub(pattern = '^chr', replacement = '', data$Chr)
  data$signature <- paste(data$Chr, data$Start, data$Ref, data$Obs, sep = '_')
  data$is.indel <- nchar(as.character(data$Ref)) > 1 | nchar(as.character(data$Obs)) > 1 | data$Ref == '-' | data$Obs == '-'
  data$indel.length <- ifelse (data$is.indel, pmax(1, nchar(as.character(data$Ref)), nchar(as.character(data$Obs))), 0) 
  if (sum(grepl(pattern = '^ENS', data$Gene)) > 0) {
    message('Replacing the Ensembl gene names with usual HUGO names')
    geneMapping <- read.table(file = biomart.gene.description, stringsAsFactors = FALSE, na.string = c('', 'NA'), sep = '\t', header = TRUE, comment.char = '|', quote= '')
    splicing.info <- ifelse( grepl(data$Gene, pattern = '\\(.*\\)$'), gsub(pattern = '.*\\(', replacement = '\\(', x = data$Gene), '') 
    gene.only <- gsub(pattern = '\\(.*', replacement = '', x = data$Gene)
    gene.only <- gsub(pattern = ';.*', replacement = '', gene.only)
    data$ensemblID <- gsub(pattern = ',.*', replacement = '', gene.only)
    data$ensemblID.bis <- ifelse (data$ensemblID != gene.only, gsub(pattern = '.*,', replacement = '', gene.only), NA) 
    data$Description <- geneMapping$description[ match(data$ensemblID, table = geneMapping$EnsemblID) ] 
    HUGO1 <- geneMapping$external_gene_id [ match(data$ensemblID, table = geneMapping$EnsemblID) ]
    HUGO2 <- geneMapping$external_gene_id [ match(data$ensemblID.bis, table = geneMapping$EnsemblID) ]
    HUGO <- ifelse (!is.na(HUGO2), paste(HUGO1, HUGO2, sep = ','), HUGO1)
    data$HUGO.no.splice.info <- HUGO
    data$HUGO <- paste(HUGO, splicing.info, sep = '')
    message('Done replacing gene names and adding description')
  } 
  ###########
  message('Defining the main rare/somewhat rare and functional flags')
  data$dup.region <- !is.na(data$SegDup) & data$SegDup >= 0.98 & ! data$Gene %in% c('GBA')
  data$somewhat.rare <- !data$dup.region & ! data$Gene %in% bad.genes & (data$cg69 < 0.1 | is.na(data$cg69)) 
  data$rare <- data$somewhat.rare 
  data$novel <- data$rare & is.na(data$cg69) 
  for (field in freq.fields) {
    message('Field ', field)
    data$somewhat.rare <- data$somewhat.rare & (data[, field] <= threshold.somewhat.rare | is.na(data[, field]) )
    data$rare <- data$rare & (data[, field] <= threshold.rare | is.na(data[, field]) )
    data$novel <- data$novel &  (data[, field] == 0 | is.na(data[, field]) )
  } 
  ######## now the functional flags
  message('Adding functional flags')
  data$exonic.splicing <- grepl(pattern = 'splicing', data$Func) & grepl(pattern = 'exonic', data$Func)
  data$splicing <- grepl(pattern = 'splicing', data$Func) &  (! grepl(pattern = 'exonic', data$Func))
  data$core.splicing <- data$splicing & grepl(pattern = '\\+1|\\-1|\\+2|\\-2', data$Gene) ##splicingVP 
  ##
  lof <- c('frameshift deletion', 'frameshift insertion', 'frameshift substitution', 'stopgain SNV', 'stoploss SNV')
  non.syn <-   c('nonsynonymous SNV', 'nonframeshift substitution', 'nonframeshift deletion', 'nonframeshift insertion')
  data$lof <- data$ExonicFunc %in% lof | data$core.splicing
  data$non.syn <- data$ExonicFunc %in% non.syn 
  ##################
  message('Now working on excluding the bad transcripts')
  data$remove.bad.transcripts <- FALSE
  if (!is.null(choice.transcripts)) { 
    my.genes <- unique(as.character(choice.transcripts$EnsemblID))
    for (gene in my.genes) { 
      transcripts <- as.character(subset(choice.transcripts, EnsemblID == gene, 'Transcript', drop = TRUE))
      bad <- (data$ensemblID == gene)  ## we can only be bad with the bad gene
      for (transcript in transcripts) {bad <- bad & ! grepl(pattern = transcript, data$AAChange) & ! grepl(pattern = transcript, data$Gene)}
      data$remove.bad.transcripts <- data$remove.bad.transcripts | bad
    }
  }
  print(table(data$remove.bad.transcripts))
  ####################
  message("Done with annotations")
  return(data)
}


process.multiVCF <- function(calls,
                             depth,
                             annotations,
                             my.cases,
                             case.control = TRUE,
                             case.control.SKAT = FALSE,
                             known.genes = c(),
                             hom.mapping = FALSE,
                             oFolder = 'processed',
                             explained.cases = NULL,
                             depth.threshold.conf.homs = 6,
                             print.individual.files = TRUE,
                             nb.IDs.to.show = 5,
                             fix.names = NULL,
                             also.fix.names.in.case.list = TRUE,
                             choice.transcripts = NULL, ##data frame, need columns "EnsemblID" and "Transcript"
                             extra.gene.annotations.file = NULL,
                             PCA.pop.structure = '/cluster/project8/vyp/exome_sequencing_multisamples/mainset/mainset_January2014/mainset_January2014_PCA.RData') {
  library(snpStats)
  require(VPlib)
  require(snpStats)
  ######### first thing first, we fix the names if there are names to fix
  if (!is.null(fix.names)) {
    #browser()
    message('Fixing the sample names first')
    if (! 'New_correct_name' %in% names(fix.names)) stop('Missing New_correct_name column in the fix name table')
    if (! 'Old_name' %in% names(fix.names)) stop('Missing Old_name column in the fix name table')
    fix.names$New_correct_name <- paste(fix.names$New_correct_name, 'z', sep = '')
    my.names <- dimnames(calls)[[1]]
    for (i in 1:nrow(fix.names)) {my.names <- ifelse ( my.names == fix.names$Old_name[ i ], fix.names$New_correct_name[ i ], my.names)}
    my.names <- gsub(my.names, pattern = 'z$', replacement = '')
    dimnames(calls)[[1]] <- my.names
    dimnames(depth)[[2]] <- my.names ### and update the depth matrix as well
    message('Done fixing the sample names')
    if (also.fix.names.in.case.list) {
      for (i in 1:nrow(fix.names)) { my.cases <- ifelse (my.cases == fix.names$Old_name[ i ], fix.names$New_correct_name[ i ], my.cases)}
      my.cases <- gsub(my.cases, pattern = 'z$', replacement = '')
      #if (sum (! my.cases %in% dimnames(calls)[[1]]) > 0) browser()
    }
  }
  ######################
  all.variants.folder <- paste(oFolder, '/all_variants', sep= '')
  hom.variants.folder <- paste(oFolder, '/homozygous_variants', sep= '')
  chrX.folder <- paste(oFolder, '/chrX', sep= '')
  hom.mapping.folder <- paste(oFolder, '/hom_mapping', sep= '')
  all.rare.variants.folder <- paste(oFolder, '/rare_variants', sep= '')
  known.genes.folder <- paste(oFolder, '/known_genes', sep= '')
  compound.hets.folder <- paste(oFolder, '/compound_hets', sep= '')
  supportFolder <- paste(oFolder, '/support', sep= '')
for (folder in c(oFolder, chrX.folder, all.variants.folder, hom.variants.folder, all.rare.variants.folder, known.genes.folder, hom.mapping.folder, compound.hets.folder, supportFolder)) {
    if (!file.exists(folder)) dir.create(folder) else {
      old.files <- list.files(folder, full.names = TRUE, include.dir = FALSE)
      message('Will remove old files')
      print(old.files)
      print(file.remove( old.files) )
    }
  }
  good.ids <- my.cases %in% dimnames(calls)[[1]]
  correct.ids <- subset(my.cases, good.ids)
  if (sum(good.ids) < dim(calls)[1]) {
    calls <- calls[ correct.ids, , drop = FALSE]
    depth <- depth[ , correct.ids , drop = FALSE]
    my.sum <- col.summary(calls)
    poly.in.cases <- (my.sum$RAF == 1 | my.sum$MAF > 0) & !is.na(my.sum$MAF) & my.sum$Calls > 0  ##this seems overly complicated. I think I really only need RAF > 0 ?
    depth <- depth[ poly.in.cases, , drop = FALSE]
    calls <- calls[, poly.in.cases, drop = FALSE ]
    annotations <- annotations[ poly.in.cases, ]
    my.sum <- my.sum[ poly.in.cases, ]
   } else {
    my.sum <- col.summary(calls)
  }
######### Update the frequency information in cases
  annotations$non.missing.cases <- my.sum$Calls
  annotations$freq.cases <- my.sum$RAF
  annotations$ncarriers.cases <- round(my.sum$Calls*(my.sum$P.AB + my.sum$P.BB))
  annotations$non.ref.calls.cases <- round( my.sum$Calls * my.sum$RAF*2  )
  annotations$missing.rate.cases <- 1 - my.sum$Call.rate
################### define the basic output folders, check what IDs are present
  if (is.null(explained.cases)) explained.cases <- rep(FALSE, length(my.cases))
  explained.cases <- subset(explained.cases, good.ids)
  nsamples <- length(correct.ids)
  message(nsamples, ' present out of ', length(my.cases), ' specified')
  print(my.cases)
  support.list <- data.frame(samples = correct.ids, explained = explained.cases)
  write.csv(x = support.list, row.names = FALSE, file = paste(supportFolder, '/list_cases.csv', sep = ''))
################ Now plot the PCA matrix to see where the samples come from
  if (!is.null(PCA.pop.structure)) {
    if (file.exists(PCA.pop.structure)) {
      load(PCA.pop.structure)
      output.pdf <- paste(supportFolder, '/PCA_plot.pdf', sep = '')
      output.tab <- paste(supportFolder, '/PCA_data.csv', sep = '')
#### Now in cases some individuals are missing from the table
      in.PCA.table <- subset(correct.ids, correct.ids %in% row.names(pcs))
      message('Output PCA plot in ', output.pdf)
      pdf(output.pdf)
      plot (x = pcs[,1],
            y = pcs[,2],
            xlab = 'PC1',
            ylab = 'PC2',
            col = 'grey',
            pch = '+')
      if (length( in.PCA.table ) > 1) {
        PCA.cases <- pcs[ in.PCA.table, ]
        points(x = PCA.cases[,1],
               y = PCA.cases[,2],
               col = 'red',
               pch = 20)
        write.csv(x = PCA.cases, file = output.tab)  ##now give the location of the cases on the PCA plot
      }
      
      ### add 1KG labels
      uniq <- unique(population$pop)
      for (i in 1:length(uniq)){
        match <-  which(population$pop == uniq[i])
        match.select <- sample(match, 1)
        text(pcs[match.select,1], pcs[match.select,2], population$pop[match.select ] ,cex=0.7, pos=4, col="black")
      }
      
      dev.off()
      
    }
  }

  

################# build the matrix of calls

####### get the sample names
  excess.message <- paste('More than', nb.IDs.to.show)
  cases.names <- dimnames(calls)[[1]]
  calls.num <- as(calls, 'numeric')

  annotations$Samples <- sapply(1:ncol(calls),
                         FUN = function(x) { if (annotations$ncarriers.cases[x] > nb.IDs.to.show) {return(excess.message)}
                         else  {paste(cases.names[ which(calls.num[,x] > 0) ], collapse = ';')}},
                         simplify = TRUE)


############## Now the case control step
  if (case.control) {
    if (sum(explained.cases) > 0) {  #### If some cases are explained we run this thing twice
      res <- case.control.analysis(choice.cases =  subset(correct.ids, !explained.cases),
                                   fix.names = fix.names,
                                   output = paste(supportFolder, '/case_control_non_explained_cases', sep = ''),
                                   known.genes = known.genes, SKAT = case.control.SKAT )
    }
        
    res1 <- case.control.analysis(choice.cases = correct.ids,
                                  fix.names = fix.names,
                                  output = paste(supportFolder, '/case_control_all_cases', sep = ''),
                                  known.genes = known.genes, SKAT = case.control.SKAT )


  }


  
################# Now start the proper filtering
  annotations$potential.comp.het <- FALSE
  annotations$n.cases.conf.homs <- 0
  summary.frame <- data.frame(ids = correct.ids,
                              n.exonic.calls = NA,
                              percent.homozyg = NA,
                              frac.het.on.X = NA,
                              n.somewhat.rare.exonic.calls = NA,
                              n.rare.exonic.calls = NA,
                              ngenes.comp.het.lof = NA,
                              ngenes.comp.het.func = NA,
                              n.func.rare.calls = NA,
                              n.func.rare.hom.calls = NA,
                              n.lof.rare.calls = NA,
                              n.lof.rare.hom.calls = NA)

  
  for (sample in 1:length(correct.ids)) {
  #for (sample in 1:1) {
    id <- correct.ids[ sample ]
    message('Sample ', id)

    depth.loc <- as(depth[, id ], 'numeric')
    calls.loc <- calls.num[id, ]
    annotations$n.cases.conf.homs <- annotations$n.cases.conf.homs + ((depth.loc >= depth.threshold.conf.homs) &  calls.loc == 2 & !is.na(calls.loc) )
    
    annotations[, id] <- paste( ifelse ( !is.na(calls.loc), c('0|0', '0|1', '1|1')[ calls.loc +1 ], "NA"), depth.loc, sep = ':')

    ###### from this point we restrict to poly in that specific case
    selected <- calls.loc > 0 & !is.na(calls.loc)
    annotations.loc <- annotations[ selected, ]
    calls.loc <- calls.loc[ selected ]
    depth.loc <- depth.loc[ selected ]
    good.hom <- (depth.loc >= depth.threshold.conf.homs) & (calls.loc == 2) & !is.na(calls.loc)

    
    
    summary.frame$n.func.rare.calls[ sample ] <- sum(annotations.loc$rare & (annotations.loc$splicing | annotations.loc$lof | annotations.loc$non.syn) & !annotations.loc$remove.bad.transcripts , na.rm = TRUE)
    summary.frame$n.func.rare.hom.calls[ sample ] <- sum(annotations.loc$rare & (annotations.loc$splicing | annotations.loc$lof | annotations.loc$non.syn) & !annotations.loc$remove.bad.transcripts  & good.hom, na.rm = TRUE)
    summary.frame$n.lof.rare.calls[ sample ] <- sum(annotations.loc$rare & annotations.loc$lof, na.rm = TRUE)
    summary.frame$n.lof.rare.hom.calls[ sample ] <- sum(annotations.loc$rare & annotations.loc$lof & good.hom, na.rm = TRUE)
    
    summary.frame$n.exonic.calls[ sample ] <- sum(annotations.loc$Func == 'exonic', na.rm = TRUE)
    summary.frame$n.rare.exonic.calls[ sample ] <- sum(annotations.loc$Func == 'exonic' & annotations.loc$rare, na.rm = TRUE)
    summary.frame$n.somewhat.rare.exonic.calls[ sample ] <- sum(annotations.loc$Func == 'exonic' & annotations.loc$somewhat.rare, na.rm = TRUE)
    summary.frame$frac.het.on.X[ sample ] <- sum( (annotations.loc$Chr == 'X') & calls.loc == 1, na.rm = TRUE)/ sum(annotations.loc$Chr == 'X' & calls.loc %in% c(1, 2), na.rm = TRUE)

######### preferred choice to print the labels
    my.block <- c('Samples', 'Func', 'ExonicFunc', 'HUGO', 'Description', 'non.ref.calls.cases', 'ncarriers.cases', 'missing.rate.cases', 'freq.controls', 'non.missing.controls', 'non.ref.calls.controls',
                  'non.missing.external.controls', 'freq.external.controls', 'AAChange', id, 'Depth', 'is.indel', 'QUAL')
    my.names2 <-  c(subset(my.block, my.block %in% names(annotations.loc)), subset(names(annotations.loc), ! (names(annotations.loc) %in% c(correct.ids, my.block) ) ))
    my.names2 <- subset(my.names2, ! my.names2 %in% c('Otherinfo', 'Gene.Start..bp.', 'Gene.End..bp.', 'MIM.Gene.Description', 'n.cases.conf.homs', 'HetNames', 'HomNames'))
##### output the full list of calls for this sample
    if (print.individual.files) {
      output.all <- paste(all.variants.folder, '/', id, '.csv', sep = '')
      write.csv(x = annotations.loc[, my.names2], file = output.all, row.names = FALSE)
      message('Outputting all variants in ', output.all, ', ncalls: ', nrow(annotations.loc))
    }

    
############# Now the compound hets
    tab.genes <- table(subset(annotations.loc$HUGO.no.splice.info, annotations.loc$somewhat.rare & (annotations.loc$splicing | annotations.loc$non.syn | annotations.loc$lof) & !annotations.loc$remove.bad.transcripts ))
    potential.comp.het.genes <- subset(tab.genes, tab.genes >= 2)
    comp.het.frame <- subset( annotations.loc, HUGO.no.splice.info %in% names(potential.comp.het.genes) & somewhat.rare & (splicing | non.syn | lof) & !annotations.loc$remove.bad.transcripts )
    summary.frame$ngenes.comp.het.func[ sample ] <- length( potential.comp.het.genes )
    
    tab.genes.lof <- table(subset(annotations.loc$HUGO.no.splice.info, annotations.loc$somewhat.rare & annotations.loc$lof & !annotations.loc$remove.bad.transcripts))
    potential.comp.het.genes.lof <- subset(tab.genes.lof, tab.genes.lof >= 2)
    summary.frame$ngenes.comp.het.lof[ sample ] <- length( potential.comp.het.genes.lof)

    output.comp.hets <- paste(compound.hets.folder, '/', id, '.csv', sep = '')
    write.csv(x = comp.het.frame[, my.names2], file = output.comp.hets, row.names = FALSE)
    message('Outputting all compound het functional variants in ', output.comp.hets, ', ncalls: ', nrow(comp.het.frame))
  
    annotations$potential.comp.het <- annotations$potential.comp.het | annotations$signature %in% comp.het.frame$signature
    
    
######## homozygosity mapping
    if (hom.mapping) {
      output.pdf <- paste(hom.mapping.folder, '/', id, '.pdf', sep = '')
      pdf(output.pdf, width = 8, height = 14)
      
      hzig.frame <- data.frame(Het = (calls.loc == 1), positions = annotations.loc$Start, chromosome = annotations.loc$Chr, depth = depth.loc, SegDup = annotations.loc$SegDup, is.indel = annotations.loc$is.indel)
      hzig.frame <- hzig.frame[ is.na(hzig.frame$SegDup) & ! hzig.frame$is.indel, ]
      
      hzig.frame <- subset(hzig.frame, chromosome %in% as.character(1:22))
      hzig.frame <- subset(hzig.frame, ! ( !Het  & depth < 7) ) ##remove low depth homozygous calls
      hzig.frame <- hzig.frame[ order(as.numeric(hzig.frame$chromosome), hzig.frame$positions), ]
      my.homozyg <- homozyg.mapping.v2 (hzig.frame$Het, positions = hzig.frame$positions, chromosome = hzig.frame$chromosome, plot = TRUE)
      dev.off()
      message('Homozygosity mapping plot in ', output.pdf)
      regions <- my.homozyg[[2]]
      
      if (nrow(regions) >= 1) {
        summary.frame$percent.homozyg[ sample ] <- sum (regions$end - regions$start) / (3*10^9)
      } else {summary.frame$percent.homozyg[ sample ] <- 0}
        
      if (nrow(regions) >= 1) {
        hom.mapping.candidates <- data.frame()
        for (region in 1:nrow(regions)) {
          hmap.loc <- subset(annotations.loc, somewhat.rare & calls.loc > 0 & (non.syn | splicing | lof) & !remove.bad.transcripts & (Start > regions$start[region]) & (Start < regions$end[region]) & (Chr == regions$chromosome[region]))
          hom.mapping.candidates <- rbind.data.frame(hom.mapping.candidates, hmap.loc)
        }
        
        hom.mapping.candidates <- hom.mapping.candidates[, my.names2]
        output.hmap <- paste(hom.mapping.folder, '/', id, '.csv', sep = '')
        write.csv(x = hom.mapping.candidates[, my.names2], file = output.hmap, row.names = FALSE)
        message('Outputting hom mapping variants in ', output.hmap, ', ncalls: ', nrow(hom.mapping.candidates))
      
        output.regions <- paste(hom.mapping.folder, '/', id, '_regions.csv', sep = '')
        write.csv(x = regions, file = output.regions, row.names = FALSE)
        message('Outputting hom mapping regions in ', output.regions, ', nregions: ', nrow(regions))      
      }
    }


##### now the somewhat.rare homs- note that I exclude chrX calls for now in that list
    rare.homs <- subset(annotations.loc, (! Chr %in% c('X', 'Y'))  & somewhat.rare & calls.loc == 2 & (non.syn | splicing | lof) & !remove.bad.transcripts & (depth.loc >= depth.threshold.conf.homs))
    rare.homs <- rare.homs[, my.names2]

    output.homs <- paste(hom.variants.folder, '/', id, '.csv', sep = '')
    write.csv(x = rare.homs, file = output.homs, row.names = FALSE)
    message('Outputting all rare homozygous functional variants in ', output.homs, ', ncalls: ', nrow(rare.homs))
    
    
##### now the somewhat.rare hets
    rare.hets <- subset(annotations.loc, somewhat.rare & calls.loc >= 1 & (non.syn | splicing | lof) & !remove.bad.transcripts)
    rare.hets <- rare.hets[, my.names2]

    output.hets <- paste(all.rare.variants.folder, '/', id, '.csv', sep = '')
    write.csv(x = rare.hets, file = output.hets, row.names = FALSE)
    message('Outputting all rare heterozygous functional variants in ', output.hets, ', ncalls: ', nrow(rare.hets))

######## now the X linked stuff
    chrX <- subset( annotations.loc, Chr == 'X' & somewhat.rare & calls.loc >= 1 & (non.syn | splicing | lof) & !remove.bad.transcripts)
    chrX <- chrX[, my.names2]
    output.X <- paste(chrX.folder, '/', id, '.csv', sep = '')
    write.csv(x = chrX, file = output.X, row.names = FALSE)
    message('Outputting all rare X linked variants ', output.X, ', ncalls: ', nrow(chrX))
    
##### now the variants in known genes, keep the wrong transcripts in this folder for now
    known.genes.calls <- subset(annotations.loc, somewhat.rare & calls.loc >= 1 & (non.syn | splicing | lof) 
                                & (Gene %in% known.genes | gsub(pattern = '\\(.*', replacement = '', annotations.loc$HUGO) %in% known.genes | HUGO %in% known.genes))
    known.genes.calls <- known.genes.calls[, my.names2]
    output.known <- paste(known.genes.folder, '/', id, '.csv', sep = '')
    write.csv(x = known.genes.calls, file = output.known, row.names = FALSE)
    message('Outputting all rare variants in known genes ', output.known, ', ncalls: ', nrow(known.genes.calls))
  }



  
############################################################ preferred choice to print the labels
  my.block <- c('Samples', 'Func', 'ExonicFunc', 'HUGO', 'Description', 'non.missing.cases', 'non.ref.calls.cases', 'ncarriers.cases', 'maf.cases', 'n.cases.conf.homs', 'freq.controls', 'non.missing.controls', 'non.ref.calls.controls', 'non.missing.external.controls', 'freq.external.controls', 'pval.cc.single', 'AAChange', 'is.indel', 'QUAL')
  my.names2 <-  c(subset(my.block, my.block %in% names(annotations)),
                  subset(names(annotations), ! (names(annotations) %in% c(correct.ids, my.block, c('HomNames', 'HetNames')) ) ),
                  correct.ids)

  #### output all sorts of overall tables
  output.all.variants <- paste(oFolder, '/all_somewhat_rare_variants_cases.csv', sep = '')
  all.variants <- subset(annotations, non.missing.cases > 0 & somewhat.rare)
  write.csv(x = all.variants[, my.names2], file = output.all.variants, row.names = FALSE)

  ##known genes
  known.genes.calls <- subset(annotations, (non.syn | splicing | lof) &
                              (HUGO %in% known.genes | Gene %in% known.genes | gsub(pattern = '\\(.*', replacement = '', annotations$HUGO) %in% known.genes ) & somewhat.rare)
  output.known.genes <- paste(oFolder, '/known_genes.csv', sep = '')
  write.csv(x = known.genes.calls[, my.names2], file = output.known.genes, row.names = FALSE)

  ##compound hets
  comp.hets <- subset(annotations, potential.comp.het)
  output.comp.hets<- paste(oFolder, '/comp_het_candidates.csv', sep = '')
  write.csv(x = comp.hets[, my.names2], file = output.comp.hets, row.names = FALSE)

  ##hom candidates
  is.hom.in.at.least.one <- apply(calls.num, MAR = 2, FUN = max, na.rm = TRUE)
  hom.candidates <- subset(annotations, is.hom.in.at.least.one == 2 & somewhat.rare & (lof | splicing | non.syn))
  output.homs <- paste(oFolder, '/hom_candidates.csv', sep = '')
  write.csv(x = hom.candidates[, my.names2], file = output.homs, row.names = FALSE)

  ##signatures
  signatures <- annotations[, c('freq.cases', 'non.missing.cases', 'signature')]
  write.csv(x = signatures, file = paste(oFolder, '/signatures.csv', sep = ''), row.names = FALSE)

  ###
  write.csv(x = summary.frame,
            file = paste(supportFolder, '/summary_frame.csv', sep = ''),
            row.names = FALSE)

  stamp.file <- paste(oFolder, '/time_stamp.txt', sep = '')
  cat(date(), file = stamp.file)
}

#### default arguments
Prion.setup <- FALSE
missing.depth.threshold <- 0
root <- '/scratch2/vyp-scratch2/vincent/GATK/mainset_October2014/mainset_October2014'

myArgs <- getArgs()

if ('Prion.setup' %in% names(myArgs)) Prion.setup <- as.logical(myArgs[[ 'Prion.setup' ]])
if ('minDepth' %in% names(myArgs)) missing.depth.threshold <- as.numeric(myArgs[[ 'minDepth' ]])
if ('root' %in% names(myArgs)) root <- as.character(myArgs[[ 'root' ]])

threshold.somewhat.rare <- 0.005
threshold.rare <- 0.001

case.names <- scan('data/caseKeywords.tab', what = character())
control.names <- scan('data/controlKeywords.tab', what = character())

first <- TRUE

for (chr in  c(as.character(1:22), 'X')) {
#for (chr in  c(as.character(16:16))) {
  message('Chromosome ', chr)
  #input.data <- paste(root, '_by_chr/chr', chr, '_snpStats.RData', sep = '')
  input.data <- paste0(root, '_snpStats/chr', chr, '_snpStats.RData')
  if (!file.exists(input.data)) {
    stop("File ", input.data, " does not exist")
  }
  if (file.exists(input.data)) {
    message('Will now load ', input.data)
    load(input.data)
    #print(table(as(matrix.calls.snpStats[, '1_897325_G_C'], 'numeric')[,1])); stop()
    if (first) {
      sample.frame <- data.frame(sampleIDs = dimnames(matrix.calls.snpStats)[[1]])
      sample.frame$case <- 0
      sample.frame$control <- 0
      sample.frame$external.control <- 0
      for (case in case.names) {
          sample.frame$case <- ifelse ( grepl(pattern = case, sample.frame$sampleIDs), 1, sample.frame$case )
      }
      for (control in control.names) {
          sample.frame$control <- ifelse ( grepl(pattern = control, sample.frame$sampleIDs), 1, sample.frame$control )
      }
      if (max(sample.frame$case + sample.frame$control ) == 2) stop('A sample is a case and a control at the same time')
################# Now we consider the Prion unit setup
      if (Prion.setup) {
        sample.frame$external.control <- ifelse (  sample.frame$case + sample.frame$control == 0, 1, 0)
      } else {
        exclude.rows <- which ( sample.frame$control + sample.frame$case == 0 )  ### that's the rows that are NOTHING AT ALL
        sample.frame.excluded.samples <- sample.frame[ exclude.rows, ]
        case.control.frame <- sample.frame[ -exclude.rows, ]  ##remove the rows that are neither from cases or controls
        write.table(x = sample.frame.excluded.samples, file = 'data/samples_not_used.tab', sep = '\t', row.names = FALSE, quote = FALSE)
        control.rows <- which(sample.frame$control == 1)
        external.controls.rows <- control.rows[ which( rep(c(0, 0, 0, 1), times = floor(length(control.rows)/4)) == 1) ]
        sample.frame$external.control[ external.controls.rows ] <- 1
        sample.frame$control [ external.controls.rows ] <- 0
      }
########### Now we work on the external set, need to define it
      cases <- subset(sample.frame, case == 1, 'sampleIDs', drop = TRUE)
      controls <- subset(sample.frame, control == 1, 'sampleIDs', drop = TRUE)
      external.controls <- subset(sample.frame, external.control == 1, 'sampleIDs', drop = TRUE)
      case.control.labels <- c(rep(1, length(cases)), rep(0, length(controls)))
    }
############# apply missingness threshold
    ##some quick checks first of all
    print(table(dimnames(matrix.depth)[[2]] == dimnames(matrix.calls.snpStats)[[1]]))
    print(table(dimnames(matrix.depth)[[1]] == dimnames(matrix.calls.snpStats)[[2]]))
    if (missing.depth.threshold > 0) {  ### apply missingness if requested
      for (i in 1:nrow(matrix.calls.snpStats)) {
        message('Fixing sample ', i, ' min depth ', missing.depth.threshold)
        matrix.calls.snpStats[ i, which(as.numeric(matrix.depth[,i]) < missing.depth.threshold) ] <- as.raw(0)
      }
    }
    ######### Now we remove the variants we don't need
    external.controls.snpStats <- matrix.calls.snpStats[ external.controls, ]
    case.control.snpStats <- matrix.calls.snpStats[ c(cases, controls), ]
    matrix.depth <- matrix.depth[ , c(cases, controls)]  ##reduce the depth matrix to the right set of cases and controls
    case.control.variant.summary <- col.summary( case.control.snpStats )
    poly.in.case.control <- case.control.variant.summary$RAF > 0 & case.control.variant.summary$Calls > 0
    case.control.snpStats <- case.control.snpStats[, poly.in.case.control]
    external.controls.snpStats <- external.controls.snpStats[, poly.in.case.control]
    annotations <- annotations.snpStats[poly.in.case.control, ]
    matrix.depth <- matrix.depth[ poly.in.case.control, ]
###########
    annotations <- annotate.standard.annovar.output(annotations)
  ################
    external.summary <- col.summary( external.controls.snpStats )
    annotations$non.ref.calls.external.controls <- round( external.summary$Calls * external.summary$RAF*2  )
    annotations$freq.external.controls <- external.summary$RAF
    annotations$non.missing.external.controls <- external.summary$Calls
    controls.snpStats <- case.control.snpStats[ controls, ]
    control.summary <- col.summary( controls.snpStats )
    annotations$non.ref.calls.controls <- round( control.summary$Calls * control.summary$RAF*2  )
    annotations$freq.controls <- control.summary$RAF
    annotations$non.missing.controls <- control.summary$Calls
    cases.snpStats <- case.control.snpStats[ cases, ]
    case.summary <- col.summary( cases.snpStats )
    annotations$non.ref.calls.cases <- round( case.summary$Calls * case.summary$RAF*2  )
    annotations$freq.cases <- case.summary$RAF
    annotations$non.missing.cases <- case.summary$Calls
############# Now use the external control data to refine the flags
    message('Number of rare variants before control filter: ', sum(annotations$rare))
    annotations$somewhat.rare <- annotations$somewhat.rare &
    ( (annotations$non.ref.calls.external.controls <= 2) |  (annotations$freq.external.controls <= threshold.somewhat.rare) | is.na(annotations$freq.external.controls) )
    annotations$rare <- annotations$rare & ( (annotations$non.ref.calls.external.controls <= 1) |  (annotations$freq.external.controls <= threshold.rare) | is.na(annotations$freq.external.controls) )
    annotations$novel <- annotations$novel  &  (annotations$freq.external.controls == 0 | is.na(annotations$freq.external.controls)) & ( is.na(annotations$freq.controls) | annotations$freq.controls == 0 )
    message('Number of rare variants after control filter: ', sum(annotations$rare))
##############  Now we need to save the full case control file
    output.file <- paste('data/case_control_chr', chr, '.RData', sep = '')
    save(list = c('case.control.labels', 'annotations', 'case.control.snpStats'), file = output.file)
########## restrict the calls to the polymorphic data in cases and only keep info about cases
    case.summary <- col.summary(case.control.snpStats[  cases, ])
    poly.in.cases <- case.summary$RAF > 0 & case.summary$Calls > 0
    annotations <- annotations[ poly.in.cases, ]
    matrix.depth <- matrix.depth[poly.in.cases, cases]
    case.control.snpStats <- case.control.snpStats[cases, poly.in.cases ]
    message('Number of polymorphic variants in cases for chromosome ', chr, ': ', sum(poly.in.cases, na.rm = TRUE) )
    if (first) {
      combined.snpStats <- case.control.snpStats
      combined.annotations <- annotations
      combined.matrix.depth <- matrix.depth
      first <- FALSE
    } else {
      combined.annotations <- rbind.data.frame( combined.annotations, annotations)    
      combined.snpStats <- cbind( combined.snpStats, case.control.snpStats)
      combined.matrix.depth <- rbind ( combined.matrix.depth, matrix.depth)
      if (nrow(combined.annotations) != ncol(combined.snpStats)) {stop('Non matching numbers')}
    }
    gc()
  }
}

message('Now saving the data in data/poly_in_cases_snpStats.RData')
save(list = c('combined.matrix.depth', 'combined.snpStats', 'combined.annotations'), file = 'data/poly_in_cases_snpStats.RData')
