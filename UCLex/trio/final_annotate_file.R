source('/cluster/project8/vyp/vincent/Software/pipeline/GATK_v2/process_multiVCF.R')


getArgs <- function () {
  myargs.list <- strsplit(grep("=", gsub("--", "", commandArgs()),
                               value = TRUE), "=")
  myargs <- lapply(myargs.list, function(x) x[2])
  names(myargs) <- lapply(myargs.list, function(x) x[1])
  return(myargs)
}


input.file <- 'processed/trio1_56_57_58/trio1_56_57_58_db.genome_summary.csv'
output.file <- 'default'


myArgs <- getArgs()

if ('input.file' %in% names(myArgs)) input.file <- myArgs[[ 'input.file' ]]
if ('output.file' %in% names(myArgs)) output.file <- myArgs[[ 'output.file' ]]


my.names <- c('Func','Gene','ExonicFunc','AAChange','Conserved','SegDup','ESP6500si_ALL','1000g2012apr_ALL','dbSNP137','AVSIFT','LJB_PhyloP','LJB_PhyloP_Pred','LJB_SIFT','LJB_SIFT_Pred','LJB_PolyPhen2','LJB_PolyPhen2_Pred','LJB_LRT','LJB_LRT_Pred','LJB_MutationTaster','LJB_MutationTaster_Pred','LJB_GERP++','cg69','Omim','Chr','Start','End','Ref','Obs','dunno','pos2','rsid2', 'Ref2' ,'Obs2', 'QUAL', 'junk1', 'junk2', 'INFO', 'Sample1', 'Sample2', 'Sample3')


data <- read.csv( input.file, col.names = my.names, na.strings = c('', 'NA'))
data <- data[, ! names(data) %in% c('pos2', 'rsid2', 'Ref2', 'Obs2', 'junk1', 'junk2', 'dunno') ]

data <- annotate.standard.annovar.output (data)
print(names(data))


data <- subset(data, somewhat.rare)


block1 <- c('HUGO', 'Description', 'Func', 'ExonicFunc', 'AAChange', 'ESP6500si_ALL', '1000g2012apr_ALL', 'dbSNP137', 'Sample1', 'Sample2', 'Sample3')
my.names <- c( subset( block1, block1 %in% names(data)), subset( names(data), ! (names(data) %in% block1)))

print(my.names)


write.csv(x = data[, my.names], file = output.file, row.names = FALSE)
