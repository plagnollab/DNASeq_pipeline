#!/usr/bin/env Rscript
err.cat <- function(x)  cat(x, '\n', file=stderr())

### Series of filters suggested by Adam.
message('*** AF FILTERING ***')
d <- read.csv(file('stdin'))

# Filtering of variants based on annotation
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tools))
suppressPackageStartupMessages(library(xtable))


option_list <- list(
    make_option(c('--exac.thresh'), default=0.01, help='pop freq threshold'),
    make_option(c('--onekg.thresh'), default=0.05, help='pop freq threshold'),
    make_option(c('--esp.thresh'), default=0.05, help='pop freq threshold'),
    make_option(c('--ajcontrols.thresh'), default=NULL, type='numeric', help='pop freq threshold'),
    make_option(c('--uclex.thresh'), default=NULL, type='numeric', help='pop freq threshold')
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

#
af.filter <- function(xx,xx.thresh) {
    xx.filter <- apply(xx, 1, function(x) {
        pop.af <- x[1]
        pop.description <- x[2]
        #message(sprintf('%s: rare (less than %.2f pct or more than %.2f pct) in %s',pop.af,xx.thresh*100,100-xx.thresh*100,pop.description))
        #err.cat(table(xx.filter <- (d[,pop.af] < xx.thresh & (d$HOM>=1|d$HET>=1)) | (d[,pop.af] > (1-xx.thresh) & (d$WT>=1 | d$HET >= 1)), useNA='always' ))
        # do not keep flipped variants
        message(sprintf('%s: rare (less than %.2f pct) in %s',pop.af,xx.thresh*100,pop.description))
        err.cat(table(xx.filter <- (d[,pop.af] < xx.thresh & (d$HOM>=1|d$HET>=1)), useNA='always'))
        return(xx.filter)
    })
    colnames(xx.filter) <- xx$pop
    i <- rowSums(xx.filter) 
    d <- d[which(is.na(i) | i==ncol(xx.filter)),]
    return(d)
}

#AJcontrols
if (!is.null(opt$ajcontrols.thresh)) {
    ajcontrols.thresh <- as.numeric(opt$ajcontrols.thresh)
    message('*** AJcontrols FILTERING ***')
    ajcontrols <- data.frame(
    pop=c('AJcontrols'),
    description=c('ajcontrols')
    )
    d <- af.filter(ajcontrols,ajcontrols.thresh)
}

#BroadAJcontrols
if (!is.null(opt$ajcontrols.thresh)) {
    ajcontrols.thresh <- as.numeric(opt$ajcontrols.thresh)
    message('*** BroadAJcontrols FILTERING ***')
    ajcontrols <- data.frame(
    pop=c('BroadAJcontrols'),
    description=c('broadajcontrols')
    )
    d <- af.filter(ajcontrols,ajcontrols.thresh)
}

#UCLEX
if (!is.null(opt$uclex.thresh)) {
    uclex.thresh <- opt$uclex.thresh
    message('*** UCLEX FILTERING ***')
    uclex <- data.frame(
    pop=c('UCLEX'),
    description=c('uclex')
    )
    d <- af.filter(uclex,uclex.thresh)
}


#EXAC
if (!is.null(opt$exac.thresh)) {
    exac.thresh <- opt$exac.thresh
    message('*** EXAC FILTERING ***')
    exac <- data.frame(
    pop=c('EXAC_AFR', 'EXAC_AMR', 'EXAC_Adj', 'EXAC_EAS', 'EXAC_FIN', 'EXAC_NFE', 'EXAC_OTH', 'EXAC_SAS'),
    description=c('african', 'american', 'adj', 'east asian', 'finnish', 'nfe', 'others', 'sas')
    )
    d <- af.filter(exac,exac.thresh)
}


#ONEKG
if (!is.null(opt$onekg.thresh)) {
    onekg.thresh <- opt$onekg.thresh
    message('*** ONEKG FILTERING ***')
    onekg <- data.frame(
    pop=c('ONEKG_EUR','ONEKG_AFR','ONEKG_AMR','ONEKG_ASN'),
    description=c('european', 'african', 'american', 'asian')
    )
    d <- af.filter(onekg,onekg.thresh)
}


#ESP
if (!is.null(opt$esp.thresh)) {
    esp.thresh <- opt$esp.thresh
    message('*** ESP FILTERING ***')
    esp <- data.frame(
    pop=c('ESP_EA','ESP_AA'),
    description=c('european-african', 'african-american')
    )
    d <- af.filter(esp,esp.thresh)
}

#c('GMAF','AFR_MAF','AMR_MAF','ASN_MAF','EUR_MAF','AA_MAF','EA_MAF')

write.csv(d, file='', quote=FALSE, row.names=FALSE)


