library(RCurl)
library(biomaRt)

getArgs <- function() {
  myargs.list <- strsplit(grep("=",gsub("--","",commandArgs()),value=TRUE),"=")
  myargs <- lapply(myargs.list,function(x) x[2] )
  names(myargs) <- lapply(myargs.list,function(x) x[1])
  return (myargs)
}

##iFolder <- "/SAN/biomed/biomed14/vyp-scratch/Zanda_Uveitis_RNASeq/processed"
iFolder <- "~/uveitis_scratch/processed"


## support.frame <- "~/Zanda_Uveitis_ConTh17_ConTh0.tab"
## output.path <- "~/uveitis_scratch/try_topGO/test_output"
db        <- "hsapiens_gene_ensembl"
mart      <- "ensembl"
num.genes <- 300
pval.thr  <- 0.00104
code <- "Zanda_uveitis"
myArgs <- getArgs()
if ('support.frame' %in% names(myArgs)) support.frame <- myArgs[['support.frame']]
if ('mart' %in% names(myArgs)) mart <- myArgs[['mart']]
if ('db' %in% names(myArgs)) db <- myArgs[['db']]
if ('code' %in% names(myArgs)) code <- myArgs[['code']]
if ('iFolder' %in% names(myArgs)) iFolder <- myArgs[['iFolder']]


###check input files and data frame
message('Now reading ', support.frame)
support <- read.table(support.frame, header = TRUE, stringsAsFactors = FALSE)

list.conditions <- grep(names(support), pattern = '^condition.*', value  = TRUE)


### deseq output folders and files

loc.pathwayGO.folder <- ''
pathwayGO.figs <- ''



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
  
################### create the appropriate folders
  

  loc.pathwayGO.folder <- paste(iFolder, '/deseq/', loc.code, '/pathwayGO',sep = '')
  pathwayGO.figs <- paste(loc.pathwayGO.folder, '/figs', sep = '')

  for (folder in c(loc.pathwayGO.folder, pathwayGO.figs)) {
    if (! file.exists(folder)) dir.create(folder)
  }
  

}
de.genes <- paste(iFolder, '/deseq/', loc.code, '/deseq_', code, '_differential_expression.tab' , sep = '')  ##this contains the key data
  

## Step 1) Read DB from biomart
###############################
message("Loading ", db, " dataset from ", mart)
my.mart <- useMart(biomart = mart, dataset = db )

GO.file <- paste(loc.pathwayGO.folder,  '/GO_', loc.code, '.RData', sep = '')
if (!file.exists(GO.file)) {
  
  tab <- read.delim(de.genes, stringsAsFactors = FALSE)
  all.ids <- unique(tab$EnsemblID)
  
  results <- getBM(attributes = c('ensembl_gene_id', 'go_id', 'name_1006'), 
                   filters = "ensembl_gene_id", 
                   values = all.ids, mart = my.mart)

  ##Remove entries with empty GO terms
  results <- results[apply(results, 1, function(x) nchar(x[2]) > 0 & nchar(x[3]) > 0), ]
  
  save(list = 'results', file = GO.file)
  
} else {
  load(GO.file)
}


## Step 2) Read table with differentially expressed genes
message("Read table with differentially expressed genes")

res <- read.delim(de.genes, stringsAsFactors = FALSE)
res <- res[, c(1, 3, 5, 2)]

##condition <- gsub("condition", "condition_", names(res)[3]) ## can prob be removed
names(res) <- c('EnsemblID', 'P.Value', 'Condition', 'GeneName')

## This line is looking for the smallest pvalue, but it also remove NAs!
res <- do.call(rbind.data.frame, by (res, IND = res$EnsemblID, FUN = function(x) {x[ which.min(x$P.Value), ]}))



## Subset first 300 go terms associated with first 300 refseq entries in res
goodGO <- unique(subset(results, ensembl_gene_id %in% res$EnsemblID[1:num.genes], 'name_1006', drop = TRUE))


final <- data.frame()


head(goodGO)

## for each go term
for (i in 1:length(goodGO)) {
  message(goodGO[ i ])

  ## Genes associated with GO name
  res$match <- res$EnsemblID %in% subset(results, name_1006 == goodGO[ i ], 'ensembl_gene_id', drop = TRUE)

  ## P-values, all, down and up regulated
  my.p.both <- fisher.test(table(res$match, res$P.Value < 0.05))$p.value
  my.p.up <- fisher.test(table(res$match, res$P.Value < 0.05 & sign(res$Condition) == 1))$p.value
  my.p.down <- fisher.test(table(res$match, res$P.Value < 0.05 & sign(res$Condition) == -1))$p.value
  

  ## Most significant p-value
  my.p <- min(c(my.p.both, my.p.up, my.p.down))

  ## if at least one significant p-value... 
  if (my.p < pval.thr && sum(res$match & res$P.Value < 0.05) > 1 ) {
    message( 'Good: ', goodGO[ i ], ' ', my.p)
    final <- rbind.data.frame(final, list(GO = goodGO[ i ], P.both = my.p.both, P.up = my.p.up, P.down = my.p.down))
    final$GO <- as.character(final$GO)
    clean.term <-  gsub(pattern = ' ', replacement = '_', goodGO[i])
    clean.term <- gsub(pattern = ',', replacement = '_', clean.term)    
    term <- goodGO[i]

    
    interesting <- subset(res, P.Value < 0.05 & EnsemblID %in% subset(results, name_1006 == goodGO[ i ], 'ensembl_gene_id', drop = TRUE))
    
    print.clean.term <- gsub(pattern = '\\/', replacement = '_', clean.term)
    output.pdf <- paste(pathwayGO.figs,'/pathway_signif_', loc.code, '_', print.clean.term, '.pdf', sep = '')
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

if(nrow(final)>0){
    
  final <- final[ order(apply(FUN = min, final[, 2:4], MAR = 1)), ]
  final$GO<- gsub(pattern = ',', replacement = '_', final$GO)
  final$P.both <- signif(final$P.both, 3)
  final$P.up <- signif(final$P.up, 3)
  final$P.down <- signif(final$P.down, 3)

  write.table(final, sep = ' &  ', row.names = FALSE, eol = '\\\\\n', quote = FALSE)
  write.csv(final, file = paste(loc.pathwayGO.folder, '/GOterms_', loc.code, '.csv', sep = ''), row.names = FALSE)
}else{
  message("Sorry no pathway detected")
}


warnings()


sessionInfo()
