

# Defining the MAF levels

MAF levels are defined based on reference sets. The first definition occurs in the crunch_data.R file, defined as follows:

```R

threshold.somewhat.rare <- 0.005
threshold.rare <- 0.001

############# Now use the external control data to refine the flags
annotations$somewhat.rare <- annotations$somewhat.rare & ( (annotations$non.ref.calls.external.controls <= 2) |  (annotations$freq.external.controls <= threshold.somewhat.rare) | is.na(annotations$freq.external.controls) )

annotations$rare <- annotations$rare & ( (annotations$non.ref.calls.external.controls <= 1) |  (annotations$freq.external.controls <= threshold.rare) | is.na(annotations$freq.external.controls) )
 
annotations$novel <- annotations$novel  &  (annotations$freq.external.controls == 0 | is.na(annotations$freq.external.controls)) & ( is.na(annotations$freq.controls) | annotations$freq.controls == 0 )

```

Following this step, these filters are precised using external datasets in the function annotate.standard.annovar.output.

```R
annotate.standard.annovar.output <- function(data,   ##this function does NOT include the control data to define the rare/somewhat.rare flags
                                             threshold.rare = 0.002,
                                             threshold.somewhat.rare = 0.005,
                                             bad.genes = c(),
                                             freq.fields = c( 'X1000g2012apr_ALL', 'ESP6500si_ALL'),
                                             choice.transcripts = NULL,
                                             biomart.gene.description = '/cluster/project8/vyp/vincent/Software/RNASeq_pipeline/bundle/human/biomart/biomart_extra_annotations_hum\
an.tab')

```

This is where the filtering actually occurs:

```R
for (field in freq.fields) {
    message('Field ', field)
    data$somewhat.rare <- data$somewhat.rare & (data[, field] <= threshold.somewhat.rare | is.na(data[, field]) )
    data$rare <- data$rare & (data[, field] <= threshold.rare | is.na(data[, field]) )
    data$novel <- data$novel &  (data[, field] == 0 | is.na(data[, field]) )
  }
```

# The rare variants category

Note that this category in fact filters on the somewhat rare filter, so at MAF < 0.5%, plus a functional filtering (non synonymous, loss-of-function or splice altering).

```R
all.rare.variants.folder <- paste(oFolder, '/rare_variants', sep= '')
rare.hets <- subset(annotations.loc, somewhat.rare & calls.loc >= 1 & (non.syn | splicing | lof) & !remove.bad.transcripts)
```