



computer <- 'CS'

if (computer == 'CS') {
  README <- '/cluster/project8/vyp/vincent/Software/pipeline/filter_rare_variants/README'
  biomart.gene.description <- '/cluster/project8/vyp/vincent/data/biomart/support/biomart_geneNames_hsapiens_gene_ensembl.tab'
  source('/cluster/project8/vyp/vincent/Software/pipeline/GATK_v2/case_control.R')
}


annotate.standard.annovar.output <- function(data,   ##this function does NOT include the control data to define the rare/somewhat.rare flags
                                             threshold.rare = 0.002,
                                             threshold.somewhat.rare = 0.005,
                                             bad.genes = c(),
                                             freq.fields = c( 'X1000g2012apr_ALL', 'ESP6500si_ALL'),
                                             choice.transcripts = NULL,
                                             biomart.gene.description = '/cluster/project8/vyp/vincent/Software/pipeline/RNASeq/bundle/human/biomart/biomart_extra_annotations_human.tab') {
  #biomart.gene.description = '/cluster/project8/vyp/vincent/data/biomart/support/biomart_geneNames_hsapiens_gene_ensembl.tab') {

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
  het.variants.folder <- paste(oFolder, '/heterozygous_variants', sep= '')
  known.genes.folder <- paste(oFolder, '/known_genes', sep= '')
  compound.hets.folder <- paste(oFolder, '/compound_hets', sep= '')
  supportFolder <- paste(oFolder, '/support', sep= '')

  for (folder in c(oFolder, chrX.folder, all.variants.folder, hom.variants.folder, het.variants.folder, known.genes.folder, hom.mapping.folder, compound.hets.folder, supportFolder)) {
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
    ##if (id == 'Vulliamy_Sample_3365') browser()
    rare.hets <- subset(annotations.loc, somewhat.rare & calls.loc >= 1 & (non.syn | splicing | lof) & !remove.bad.transcripts)
    rare.hets <- rare.hets[, my.names2]

    output.hets <- paste(het.variants.folder, '/', id, '.csv', sep = '')
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
