##create a snpStats object with annotations (including external control set) that has only the variant found in at least one case

###TODO: also save individual per chromosome files with case, controls and basic annotations-> useful for SKAT test
###TODO: create equivalent files with read depth information

library(snpStats)
source('/cluster/project8/vyp/vincent/Software/pipeline/GATK_v2/process_multiVCF.R')
options(stringsAsFactors = FALSE)


getArgs <- function() {
  myargs.list <- strsplit(grep("=",gsub("--","",commandArgs()),value=TRUE),"=")
  myargs <- lapply(myargs.list,function(x) x[2] )
  names(myargs) <- lapply(myargs.list,function(x) x[1])
  return (myargs)
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
  input.data <- paste(root, '_by_chr/chr', chr, '_snpStats.RData', sep = '')

  if (file.exists(input.data)) {
    message('Will now load ', input.data)
    load(input.data)
    #print(table(as(matrix.calls.snpStats[, '1_897325_G_C'], 'numeric')[,1])); stop()
    
    if (first) {
      sample.frame <- data.frame(sampleIDs = dimnames(matrix.calls.snpStats)[[1]])
      sample.frame$case <- 0
      sample.frame$control <- 0
      sample.frame$external.control <- 0
      
      for (case in case.names) {sample.frame$case <- ifelse ( grepl(pattern = case, sample.frame$sampleIDs), 1, sample.frame$case )}
      for (control in control.names) {sample.frame$control <- ifelse ( grepl(pattern = control, sample.frame$sampleIDs), 1, sample.frame$control )}
      
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
