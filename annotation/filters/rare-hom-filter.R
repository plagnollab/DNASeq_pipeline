#!/usr/bin/env Rscript

# Filtering of variants based on annotation
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tools))
suppressPackageStartupMessages(library(xtable))


option_list <- list(
    make_option(c('--sample'), default=NULL, help='sample'),
    make_option(c('--exac.thresh'), default=0.0001, help='pop freq threshold'),
    make_option(c('--onekg.thresh'), default=0.0001, help='pop freq threshold'),
    make_option(c('--esp.thresh'), default=0.0001, help='pop freq threshold'),
    make_option(c('--ajcontrols.thresh'), default=NULL, type='numeric', help='pop freq threshold'),
    make_option(c('--uclex.thresh'), default=0.0001, type='numeric', help='pop freq threshold'),
    make_option(c('--depth.thresh'), default=10, type='numeric', help='depth threshold'),
    make_option(c('--out'), help='outfile'),
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

sample <- opt$sample

variants <- data.frame()
for (chr in 1:22) {
    f <- sprintf('VEP_%s-genotypes.csv',chr)
    geno <- read.csv(f)
    variants <- rbind( variants, geno )
}


ann.variants <- data.frame()
for (chr in 1:22) {
    f <- sprintf('VEP_%s-annotations.csv',chr)
    (ann <- read.csv(f))
    rownames(ann) <- ann$VARIANT_ID
    ann.variants <- rbind(ann.variants, ann)
}

variants <- cbind(ann.variants[variants$VARIANT_ID,],variants[,-1])

depth.variants <- data.frame()
for (chr in 1:22) {
    f <- sprintf('VEP_%s-genotypes_depth.csv',chr)
    depth <- read.csv(f)
    colnames(depth) <- paste(colnames(depth),'depth',sep='.')
    rownames(depth.variants) <- depth.variants$VARIANT_ID
    depth.variants <- rbind(depth.variants, depth)
}

variants <- cbind(variants,depth.variants[variants$VARIANT_ID,][,-1])

variants <- variants[,unique(colnames(variants))]

#remove these columns:
#DISTANCE CADD_RAW  CADD_PHRED
variants <- variants[,-which(colnames(variants) %in% c('DISTANCE','CADD_RAW','CADD_PHRED'))]

variants <- variants[order(variants$CADD,decreasing=TRUE),]

d <- variants

# filter variants with Mendelian Error
if ('ME' %in% colnames(d)) {
    d <- d[which(!d[,'ME']),]
}

#
af.filter <- function(xx,xx.thresh) {
    xx.filter <- apply(xx, 1, function(x) {
        pop.af <- x[1]
        pop.description <- x[2]
        xx.filter <- ( d[,pop.af] < xx.thresh | is.na(d[,pop.af]) )
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
    ajcontrols <- data.frame(
    pop=c('AJcontrols'),
    description=c('ajcontrols')
    )
    d <- af.filter(ajcontrols,ajcontrols.thresh)
}

#BroadAJcontrols
if (!is.null(opt$ajcontrols.thresh)) {
    ajcontrols.thresh <- as.numeric(opt$ajcontrols.thresh)
    ajcontrols <- data.frame(
    pop=c('BroadAJcontrols'),
    description=c('broadajcontrols')
    )
    d <- af.filter(ajcontrols,ajcontrols.thresh)
}

#UCLEX
if (!is.null(opt$uclex.thresh)) {
    uclex.thresh <- opt$uclex.thresh
    uclex <- data.frame(
    pop=c('UCLEX'),
    description=c('uclex')
    )
    d <- af.filter(uclex,uclex.thresh)
}


#EXAC
if (!is.null(opt$exac.thresh)) {
    exac.thresh <- opt$exac.thresh
    exac <- data.frame(
    pop=c('EXAC_AFR', 'EXAC_AMR', 'EXAC_Adj', 'EXAC_EAS', 'EXAC_FIN', 'EXAC_NFE', 'EXAC_OTH', 'EXAC_SAS'),
    description=c('african', 'american', 'adj', 'east asian', 'finnish', 'nfe', 'others', 'sas')
    )
    d <- af.filter(exac,exac.thresh)
}


#ONEKG
if (!is.null(opt$onekg.thresh)) {
    onekg.thresh <- opt$onekg.thresh
    onekg <- data.frame(
    pop=c('ONEKG_EUR','ONEKG_AFR','ONEKG_AMR','ONEKG_ASN'),
    description=c('european', 'african', 'american', 'asian')
    )
    d <- af.filter(onekg,onekg.thresh)
}


#ESP
if (!is.null(opt$esp.thresh)) {
    esp.thresh <- opt$esp.thresh
    esp <- data.frame(
    pop=c('ESP_EA','ESP_AA'),
    description=c('european-african', 'african-american')
    )
    d <- af.filter(esp,esp.thresh)
}


# all variant should have depth greater than depth.thresh
i <- grep('.depth',colnames(d))
d <- d[which(rowSums(d[,i]>opt$depth.thresh)==length(i)),]

# not seen in any of the other individuals in the batch
i <- setdiff(colnames(geno),c(sample,'VARIANT_ID'))
d <- d[which(rowSums(d[,i]==2)==0),]


variants <- d


write.csv( variants[which(variants[,sample]==2),], file=opt$out, quote=FALSE, row.names=FALSE  )






