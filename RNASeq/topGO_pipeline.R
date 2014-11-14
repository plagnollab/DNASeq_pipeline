## We want to look for gene set enrichment

## topGO required data:
## 
## - a) List of gene identifiers and optionally the gene-wise
##   scores. The score can be the t-test statistic (or the p-value)
##   for differential expression, correlation with a phenotype, or any
##   other relevant score.
##
## - b) Mapping between gene identifiers and GO terms.
##
## - c) GO hierarchical structure. This structure is obtained from the
##   GO.db package. At the moment topGO supports only the ontology
##   definition provided by GO.db.


### TODO LIST ###
## Use Adj p-values (need to catch empty objects)
## Use p-vals + DE info
## Replicate pathwayGO (long!)
## try different params like node.size



## we better keep the data in data frames as strings
options(stringsAsFactors = FALSE)

library(RCurl)
library(biomaRt)


## Load topGO library


## library("graph", lib.loc="~/software/R/x86_64-unknown-linux-gnu-library/2.15/")
## library("SparseM", lib.loc="~/software/R/x86_64-unknown-linux-gnu-library/2.15/")
##library("topGO", lib.loc="~/software/R/x86_64-unknown-linux-gnu-library/2.15/")
library("topGO")


## support.frame='~/Zanda_Uveitis_RNASeq/manu/data/Zanda_Uveitis_ConTh17_ConTh0.tab'
## code='Zanda_uveitis'
## mart='ensembl'
## db='hsapiens_gene_ensembl'
## iFolder='/SAN/biomed/biomed14/vyp-scratch/Zanda_Uveitis_RNASeq/processed'

############# functions ##################
topDiffGenes <- function(allScore) {
  return(allScore < 0.05)
}


colMap <- function(x) {
  .col <- rep(rev(heat.colors(length(unique(x)))), time = table(x))
  return(.col[match(1:length(x), order(x))])
}


## Sanity check on KS test
sanityCheck_KS <- function(loc.code, res, results, GOdata, resultsKS, resultsKS.elim, pdf=TRUE, csv=TRUE){

  pValue.classic <- score(resultKS)
  pValue.elim <- score(resultKS.elim)[names(pValue.classic)]
  gstat <- termStat(GOdata, names(pValue.classic))
  gSize <- gstat$Annotated / max(gstat$Annotated) * 4

  gCol <- colMap(gstat$Significant)

  if(pdf== TRUE){
    pdf.file <- paste(loc.topGO.folder, '/topGO_', loc.code, '_KSsanityCheck.pdf', sep="")
  ## sanity check classic v elim method
    pdf(pdf.file)
    plot(pValue.classic, pValue.elim, xlab = "KS test p-value classic", ylab = "KS test p-value elim",
         pch = 19, cex = gSize, col = gCol)
    dev.off()
  }
  ## select outliers 
  sel.go <- names(pValue.classic)[pValue.elim < pValue.classic]
  if(length(sel.go) > 0){
    outliers <- cbind(termStat(GOdata, sel.go),
                      elim = pValue.elim[sel.go],
                      classic = pValue.classic[sel.go])
    outliers$go_id <- row.names(outliers)
  
    outliers <- merge(outliers, results)
    outliers <- merge(outliers, res[, c('EnsemblID', 'GeneName')], by.x="ensembl_gene_id", by.y="EnsemblID")
    outliers <- outliers[order(outliers$name_1006, outliers$go_id), ]
    print(outliers)
    ## save outliers
    if(csv == TRUE){
      outliers.file <- paste(loc.topGO.folder, '/topGO_KSsanityCheck_outliers.csv', sep="")
      write.table(outliers, outliers.file, sep=",", quote=FALSE, row.names=FALSE)
    }
  }

}


plotGenes <- function(goodGO, test="Fisher", res, results, loc.topGO.folder, topGO.figs, loc.code){

  final <- data.frame()

  ## for each go term
  for (i in 1:length(goodGO)) {
    message(goodGO[ i ])

    ## Genes associated with GO name
    res$match <- res$EnsemblID %in% subset(results, name_1006 == goodGO[ i ], 'ensembl_gene_id', drop = TRUE)
    head(res)
    message('Num gene matches: ', sum(as.integer(res$match)))

    all.size  <- prod(dim(table(res$match, res$P.Value < 0.05)))
    up.size   <- prod(dim(table(res$match, res$P.Value < 0.05 & sign(res$Condition) == 1)))
    down.size <- prod(dim(table(res$match, res$P.Value < 0.05 & sign(res$Condition) == -1)))

    my.p.both <- my.p.up <- my.p.down <- NA
 
    my.p <- NULL
    ## P-values, all, down and up regulated
    if(all.size == 4){
      my.p.both <- fisher.test(table(res$match, res$P.Value < 0.05))$p.value
      my.p <- c(my.p, my.p.both)

      if(up.size == 4){
        my.p.up <- fisher.test(table(res$match, res$P.Value < 0.05 & sign(res$Condition) == 1))$p.value
        my.p <- c(my.p, my.p.up)
      }else{
        message("No up-regulated genes are significant for this pathway")
      }
      if(down.size == 4){
        my.p.down <- fisher.test(table(res$match, res$P.Value < 0.05 & sign(res$Condition) == -1))$p.value
        my.p <- c(my.p, my.p.down)
      } else{
        message("No down-regulated genes are significant for this pathway")
      }
    }else{
      message("No genes are significant for this pathway")
    }

    ## Most significant p-value
    my.p <- min(my.p)
    print(my.p)
    

    ##my.p <- min(c(my.p.both, my.p.up, my.p.down))
    
    ## if at least one significant p-value... 
    if (!is.na(my.p) && my.p < pval.thr && sum(res$match & res$P.Value < 0.05) > 1 ) {
      message( 'Good: ', goodGO[ i ], ' ', my.p)
      final <- rbind.data.frame(final, list(GO = goodGO[ i ], P.both = my.p.both, P.up = my.p.up, P.down = my.p.down))
      final$GO <- as.character(final$GO)
      clean.term <-  gsub(pattern = ' ', replacement = '_', goodGO[i])
      clean.term <- gsub(pattern = ',', replacement = '_', clean.term)
      term <- goodGO[i]

      interesting <- subset(res, P.Value < 0.05 & EnsemblID %in% subset(results, name_1006 == goodGO[ i ], 'ensembl_gene_id', drop = TRUE))

      print.clean.term <- gsub(pattern = '\\/', replacement = '_', clean.term)
      output.pdf <- paste(topGO.figs,'/topGO_',  loc.code, '_pathway_signif_', test, '_', print.clean.term, '.pdf', sep = '')
      ##output.pdf <- gsub(pattern = '\\/', replacement = '_', output.pdf)
      pdf(output.pdf, width = 12, height = 4)
      plot (x = NA, y = NA, ylim = c(-1, 1), xlim = c(-4, 4),
            xlab = 'log2 fold change',
            ylab = '',
            yaxt = 'n',
            sub = paste(nrow(interesting), 'significant at P < 0.05 out of ', sum(res$match), ', Fisher P =', signif(my.p, 3)))
    
      points (x = subset(res$Condition, res$match & res$P.Value >= 0.05), y = rep(0, sum(res$match & res$P.Value >= 0.05 )), col = 'black', pch = '+')
      points (x = interesting$Condition, y = rep(0, nrow(interesting)), col = 'red', pch = '+')
      
    my.y <- c(-1, -0.8, -0.6, -0.4, -0.2, 0.2, 0.4, 0.6, 0.8, 1) [ (1:nrow(interesting)) %% 10  + 1]
      text(x = interesting$Condition,
           y = my.y,
           labels = as.character(interesting$GeneName),
           cex = 0.9)
      
      title(paste(loc.code, ', GO search term:', term))
      dev.off()
      
    }
  }
  
  if(nrow(final) > 0){
    final <- final[ order(apply(FUN = min, final[, 2:4], MAR = 1)), ]

    final$GO<- gsub(pattern = ',', replacement = '_', final$GO)
    final$P.both <- signif(final$P.both, 3)
    final$P.up <- signif(final$P.up, 3)
    final$P.down <- signif(final$P.down, 3)

    write.table(final, sep = ' &  ', row.names = FALSE, eol = '\\\\\n', quote = FALSE)
    write.csv(final, file = paste(loc.topGO.folder, '/GOterms_with_significant_genes_', test, '.csv', sep = ''), row.names = FALSE)
  }else{
    message("No significant genes")
  }
  
}


## loc.code=loc.code
## iFolder=iFolder
## results=results
## res=res
## GOdata=GOdata
## result.test=resultFisher
## test.type="Fisher"
## plot.nodes=5
## top.nodes=20

analyseTest <- function(loc.code, iFolder, results, res, GOdata, result.test, test.type="Fisher", plot.nodes=5, top.nodes=20){
  
  loc.topGO.folder <- paste(iFolder, '/deseq/', loc.code, '/topGO',sep = '')
  topGO.figs <- paste(loc.topGO.folder, '/figs', sep = '')


 
  ## Plot graph
  ##  pdf.plot <- paste(loc.topGO.folder, "/topGO_", loc.code, "_graph_top_", test.type, sep="")
  ##   printGraph(GOdata,
  ##              result.test, 
  ##              firstSigNodes = plot.nodes,
  ##              fn.prefix = pdf.plot,
  ##              useInfo = "all",
  ##              pdfSW=TRUE
  ##              )
  
  ## Create table top nodes
  my.res.table <- GenTable(GOdata, classic = result.test, topNodes = top.nodes, numChar=100, ranksOf="classic")

  ## mapping with go terms
  mapping <- unique(results[, c('go_id', 'name_1006')])
  names(mapping) <- c("GO.ID", "GO_term")
  mapping <- subset(mapping, GO.ID %in% my.res.table$GO.ID)
  
  ## Merge with mapping to do some formatting
  my.res.table <- merge(my.res.table, mapping, all.x=TRUE)
  my.res.table$GO_term[is.na(my.res.table$GO_term)] <- my.res.table$Term[is.na(my.res.table$GO_term)]
  my.res.table$GO_term <-  gsub(pattern = ' ', replacement = '_', my.res.table$GO_term)

  output.table <- my.res.table[, c('GO.ID', 'GO_term', 'Annotated', 'Significant', 'Expected', 'classic')]

output.table$GO_term<- gsub(pattern = ',', replacement = '_', output.table$GO_term)

  names(output.table) <-  gsub("\\.", "_",  names(output.table))
  sigpath.file <- paste(loc.topGO.folder, '/topGOterms_', test.type, ".csv", sep="")
  write.table(output.table, sigpath.file, sep=",", quote=FALSE, row.names=FALSE)

  goodGO <- unique(results[results$go_id %in% my.res.table$GO.ID, 'name_1006'])
  head(goodGO)
  plotGenes(goodGO, test=test.type, res, results, loc.topGO.folder, topGO.figs, loc.code)




  
}


##############################
getArgs <- function() {
  myargs.list <- strsplit(grep("=",gsub("--","",commandArgs()),value=TRUE),"=")
  myargs <- lapply(myargs.list,function(x) x[2] )
  names(myargs) <- lapply(myargs.list,function(x) x[1])
  return (myargs)
}

iFolder <- "/SAN/biomed/biomed14/vyp-scratch/Zanda_Uveitis_RNASeq/processed"
##iFolder <- "~/uveitis_scratch/processed"

### input/output directory to be replaced with script arguments
##input.file <- "~/uveitis_scratch/processed/deseq/report/ConTh17_ConTh0/deseq_Zanda_uveitis_differential_expression.tab"

support.frame <- "/cluster/project4/vyp/Zanda_Uveitis_RNASeq/manu/data/Zanda_Uveitis_ConTh17_ConTh0.tab"
output.path <- "~/uveitis_scratch/try_topGO/test_output"
db <- "hsapiens_gene_ensembl"
mart    <- "ensembl"
num.genes <- 300
pval.thr <- 0.00104
use.adj.pval <- FALSE
node.size = 10
code <- "Zanda_uveitis"

myArgs <- getArgs()
if ('support.frame' %in% names(myArgs)) support.frame <- myArgs[['support.frame']]
if ('mart' %in% names(myArgs)) mart <- myArgs[['mart']]
if ('db' %in% names(myArgs)) db  <- myArgs[['db']]
if ('code' %in% names(myArgs)) code <- myArgs[['code']]
if ('iFolder' %in% names(myArgs)) iFolder <- myArgs[['iFolder']]
##if ('use.adj.pval' %in% names(myArgs)) use.adj.pval <- myArgs[['use.adj.pval']]


###check input files and data frame
message('Now reading ', support.frame)
support <- read.table(support.frame, header = TRUE, stringsAsFactors = FALSE)

list.conditions <- grep(names(support), pattern = '^condition.*', value  = TRUE)


loc.topGO.folder <- ''
topGO.figs <- ''

###loop over all proposed conditions
for (condition in list.conditions) {

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
  message('Support data frame: ', loc.code)
  print(support.loc)
  
################### create appropriate  output folders

  loc.topGO.folder <- paste(iFolder, '/deseq/', loc.code, '/topGO',sep = '')
  topGO.figs <- paste(loc.topGO.folder, '/figs', sep = '')

  for (folder in c(loc.topGO.folder, topGO.figs)) {
    if (! file.exists(folder)) dir.create(folder)
  }

}

de.genes <- paste(iFolder, '/deseq/', loc.code, '/deseq_', code, '_differential_expression.tab' , sep = '')  ##this contains the key data
  

## Step 1) Read DB from biomart or from input file
##################################################

GO.file <- paste(loc.topGO.folder,  '/topGO_', loc.code, '.RData', sep = '')
if (!file.exists(GO.file)) {
  
  tab <- read.delim(de.genes, stringsAsFactors = FALSE)
  all.ids <- unique(tab$EnsemblID)

  message("Loading ", db, " dataset from ", mart)
  my.mart <- useMart(biomart = mart, dataset = db )

  results <- getBM(attributes = c('ensembl_gene_id', 'go_id', 'name_1006'), 
                   filters = "ensembl_gene_id", 
                   values = all.ids, mart = my.mart)

  ##Remove entries with empty GO terms
  results <- results[apply(results, 1, function(x) nchar(x[2]) > 0 & nchar(x[3]) > 0), ]
  
  ## Convert data frame to list
  ensembl2goID <- by(results$go_id, results$ensembl_gene_id, function(x) as.character(x))
  ##ensembl2goID <- split(results$go_id, as.factor(results$ensembl_gene_id), drop=TRUE)

  ##ensembl2goName <- split(results$name_1006, as.factor(results$ensembl_gene_id), drop=TRUE)
  ensembl2goName <- by(results$name_1006, results$ensembl_gene_id, function(x) as.character(x))

  save(list = c('results', 'ensembl2goID', 'ensembl2goName'), file = GO.file)
  
} else {
  load(GO.file)
}






## Step 2) Load gene data
#########################

message("Read table with differentially expressed genes")

res <- read.delim(de.genes, stringsAsFactors = FALSE)
res <- res[, c(1, 3, 4, 5, 2)]

##condition <- gsub("condition", "condition_", names(res)[3]) ## can prob be removed
names(res) <- c('EnsemblID', 'P.Value', 'Adj.P.Value', 'Condition', 'GeneName')

## This line is looking for the smallest pvalue, but it also remove NAs!
res <- do.call(rbind.data.frame, by (res, IND = res$EnsemblID, FUN = function(x) {x[ which.min(x$P.Value), ]}))


### Step 3) GSEA
################


## 3a) Subset relevant signal
#########################
pvalues <- NULL
if(use.adj.pval == TRUE){
  pvalues <- res$Adj.P.Value
  names(pvalues) <- res$EnsemblID
}else{
  pvalues <- res$P.Value
  names(pvalues) <- res$EnsemblID
}


message("N=", sum(topDiffGenes(pvalues)), " genes with ", ifelse(use.adj.pval, "adjusted", "basic"), " p-values < 0.05")

## 3b) Create object
####################
GOdata <- new("topGOdata", 
                    description = loc.code,  ontology = "BP",
                    allGenes = pvalues, geneSel = topDiffGenes,
                    nodeSize = node.size,
                    annot = annFUN.gene2GO, gene2GO=ensembl2goID)

## 3c) Run tests
################
resultFisher  <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultKS      <- runTest(GOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")

## 3d) Bind tests together
##########################
all.tests.file <- paste(loc.topGO.folder, '/topGO_top20terms_allTests.csv', sep="")
allRes <- GenTable(GOdata, classicFisher = resultFisher,
                   classicKS = resultKS,  elimKS = resultKS.elim,
                   orderBy = "classicFisher", ranksOf = "elimKS", topNodes = 20)


allRes$Term<- gsub(pattern = ',', replacement = '_', allRes$Term)
print(allRes)
write.table(allRes, all.tests.file, quote=FALSE, row.names=FALSE, sep=",")

## 3e) run sanity check on KS test
##################################
sanityCheck_KS(loc.code, res, results, GOdata, resultsKS, resultsKS.elim, pdf=TRUE, csv=TRUE)


## 3f) analyse results
######################
analyseTest(loc.code, iFolder, results, res, GOdata, resultFisher, test.type="Fisher", plot.nodes=5, top.nodes=20)
analyseTest(loc.code, iFolder, results, res, GOdata, resultKS, test.type="KS", plot.nodes=5, top.nodes=20)
analyseTest(loc.code, iFolder, results, res, GOdata, resultKS.elim, test.type="KSelim", plot.nodes=5, top.nodes=20) 

warnings()


sessionInfo()


