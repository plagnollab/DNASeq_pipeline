# Process the output of GATK PhaseByTransmission.
# Rscript filter_de_novo.file.R --input.file <input> --output.file <output>

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(xtable))

option_list <- list(
    make_option(c('--input.file'), default='processed/trio1_56_57_58/trio1_56_57_58_noMendel.tab', help=''),
    make_option(c('--output.file'), default='processed/trio1_56_57_58/trio1_56_57_58_noMendel_clean.tab', help='')
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)


input.file <- opt$input.file
output.file <- opt$output.file


data <- read.table(input.file, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
data <- subset(data,  FATHER_GT != "./." & MOTHER_GT != "./." & CHILD_GT != "./.")
data <- subset(data, MOTHER_GT == FATHER_GT)


data$MOTHER_AD.1 <- as.numeric(gsub(data$MOTHER_AD, pattern = ',.*', replacement = ''))
data$MOTHER_AD.2 <- as.numeric(gsub(data$MOTHER_AD, pattern = '.*,', replacement = ''))
data$MOTHER.ratio <- pmin ( data$MOTHER_AD.1 / (data$MOTHER_AD.1 + data$MOTHER_AD.2), data$MOTHER_AD.2 / (data$MOTHER_AD.1 + data$MOTHER_AD.2))

data$FATHER_AD.1 <- as.numeric(gsub(data$FATHER_AD, pattern = ',.*', replacement = ''))
data$FATHER_AD.2 <- as.numeric(gsub(data$FATHER_AD, pattern = '.*,', replacement = ''))
data$FATHER.ratio <- pmin ( data$FATHER_AD.1 / (data$FATHER_AD.1 + data$FATHER_AD.2), data$FATHER_AD.2 / (data$FATHER_AD.1 + data$FATHER_AD.2))

data$CHILD_AD.1 <- as.numeric(gsub(data$CHILD_AD, pattern = ',.*', replacement = ''))
data$CHILD_AD.2 <- as.numeric(gsub(data$CHILD_AD, pattern = '.*,', replacement = ''))
data$CHILD.ratio <- pmin ( data$CHILD_AD.1 / (data$CHILD_AD.1 + data$CHILD_AD.2), data$CHILD_AD.2 / (data$CHILD_AD.1 + data$CHILD_AD.2))


data <- subset(data, MOTHER_AD.1  + MOTHER_AD.2 > 6)
data <- subset(data, FATHER_AD.1  + FATHER_AD.2 > 6)
data <- subset(data, MOTHER.ratio < 0.05 & FATHER.ratio < 0.05 & CHILD.ratio > 0.2 & CHILD_AD.2 > 5)


write.table( x = data, file = output.file, sep = '\t', row.names = FALSE, quote = FALSE, col.names = TRUE)
write.table( x = data[, c('CHROM', 'POS') ], file = paste(output.file, 'positions', sep = '.'), sep = '\t', row.names = FALSE, quote = FALSE, col.names = FALSE)
message('Output in ', output.file)


