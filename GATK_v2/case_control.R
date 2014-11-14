
case.control.analysis <- function(choice.cases = NULL, output = 'processed/support/case_control', known.genes = c(), SKAT = FALSE, fix.names = NULL,
                                  annotations.file = '/cluster/project8/vyp/vincent/Software/pipeline/RNASeq/bundle/human/biomart/biomart_extra_annotations_human.tab') {
  message('Now running the case control script')
  
  options(stringsAsFactors = FALSE)
  
  library(snpStats)
  if (SKAT) library(SKAT)
  
  final.frame <- NULL
  first <- TRUE

  output.file <- paste(output, '_gene_based_summary.csv', sep = '')
  output.single <- paste(output, '_single_variant_summary.csv', sep = '')
  
  gene.annotations <- read.table(file = annotations.file, stringsAsFactors = FALSE, na.string = c('', 'NA'), sep = '\t', header = TRUE,
                                 comment.char = '|', quote= '')
  
  final.table <- data.frame()
  final.single <- data.frame()
  
  for (chr in as.character(1:22)) {
                                        #for (chr in as.character(16)) {
    
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
      
      my.names <- c('fisher.pvalue', 'dbSNP137', 'HUGO', 'ExonicFunc', 'Description', 'AAChange', 'clean.signature', 'freq.controls', 'freq.cases', 'rare', 'somewhat.rare', 'lof', 'non.syn', 'splicing', 'ESP6500si_ALL',
                    'X1000g2012apr_ALL', 'Chr', 'Start', 'End', 'Ref', 'Obs', 'FILTER', 'non.missing.controls', 'non.ref.calls.controls', 'non.missing.cases', 'non.ref.calls.cases')
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
            
            chr.table$case.control.pval.biallelic.func[ i ] <- pbinom(q = chr.table$cases.biallelic.n.somewhat.rare.func[ i ] - 1,
                                                                      size = chr.table$cases.biallelic.n.somewhat.rare.func[ i ] + chr.table$controls.biallelic.n.somewhat.rare.func[ i ],
                                                                      prob = proba.mut,
                                                                      lower.tail = FALSE)
            
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
          chr.table$case.control.pval.lof[ i ] <- pbinom(q = chr.table$cases.n.rare.lof[ i ] - 1,
                                                         size = chr.table$cases.n.rare.lof[ i ] + chr.table$controls.n.rare.lof[ i ],
                                                         prob = proba.mut,
                                                         lower.tail = FALSE)
        } else {
          chr.table$case.control.pval.lof[ i ] <- 1
        }
      }
      
      print(subset(chr.table, case.control.pval.func < 0.01 | case.control.pval.lof < 0.01))
      final.table <- rbind.data.frame ( final.table, chr.table)
      
      #write.csv(x = final.table[ order(pmin(final.table$case.control.pval.lof, final.table$case.control.pval.func), decreasing = FALSE), ],
      #          file = output.file, row.names = FALSE)

      ### order the gene based list by recessive P-values
      write.csv(x = final.table[ order(final.table$case.control.pval.biallelic.func, decreasing = FALSE), ],
                file = output.file, row.names = FALSE)
      
      write.csv( x = final.single[ order(final.single$fisher.pvalue, decreasing = FALSE), ], file = output.single, row.names = FALSE)
      
    }
  }
    
}


#test <- case.control.analysis (choice.cases = NULL, oFolder = 'processed/support', known.genes = c())
