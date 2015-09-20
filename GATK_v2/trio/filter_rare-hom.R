# Filtering of variants based on annotation
suppressPackageStartupMessages(library(xtable))

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tools))
suppressPackageStartupMessages(library(xtable))


option_list <- list(
    make_option(c('--chr'), default=NULL, help='chromosome'),
    make_option(c('--ped'), default='pedigree_details.ped', help='pedigree_details.ped'),
    make_option(c('--exac.thresh'), default=0.0001, help='pop freq threshold'),
    make_option(c('--onekg.thresh'), default=0.0001, help='pop freq threshold'),
    make_option(c('--esp.thresh'), default=0.0001, help='pop freq threshold'),
    make_option(c('--uclex.thresh'), default=0.0001, type='numeric', help='pop freq threshold'),
    make_option(c('--depth.thresh'), default=20, type='numeric', help='depth threshold')
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

#
exac.thresh <- opt$exac.thresh
onekg.thresh <- opt$onekg.thresh
esp.thresh <- opt$esp.thresh
uclex.thresh <- opt$uclex.thresh
depth.thresh <- opt$depth.thresh

ped <- read.table(opt$ped,col.names=c('fid','id','dadid','momid','gender','affection'))

# combined file containing annotations, genotypes and genotype.depths
d <- read.csv(sprintf('VEP_%s.csv',opt$chr))
rownames(d) <- d$VARIANT_ID

d$chr <- unlist(lapply(strsplit(d$VARIANT_ID,'_'),`[[`,1))
d$pos <- as.numeric(lapply(strsplit(d$VARIANT_ID,'_'),`[[`,2))

geno.id <- paste('geno',ped[1,'id'],sep='.')
geno.momid <- paste('geno',ped[1,'momid'],sep='.')
geno.dadid <- paste('geno',ped[1,'dadid'],sep='.')

depth.id <- paste('depth',ped[1,'id'],sep='.')
depth.momid <- paste('depth',ped[1,'momid'],sep='.')
depth.dadid <- paste('depth',ped[1,'dadid'],sep='.')

#remove these columns:
#DISTANCE CADD_RAW  CADD_PHRED AF
d <- d[,-which(colnames(d) %in% c('DISTANCE','CADD_RAW','CADD_PHRED','AF'))]
d <- d[order(d$CADD,decreasing=TRUE),]

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

#UCLEX
uclex <- data.frame(
    pop=c('UCLEX'),
    description=c('uclex')
    )
d <- af.filter(uclex,uclex.thresh)

#EXAC
exac <- data.frame(
    pop=c('EXAC_AFR', 'EXAC_AMR', 'EXAC_Adj', 'EXAC_EAS', 'EXAC_FIN', 'EXAC_NFE', 'EXAC_OTH', 'EXAC_SAS'),
    description=c('african', 'american', 'adj', 'east asian', 'finnish', 'nfe', 'others', 'sas')
    )
d <- af.filter(exac,exac.thresh)

#ONEKG
onekg <- data.frame(
    pop=c('ONEKG_EUR','ONEKG_AFR','ONEKG_AMR','ONEKG_ASN'),
    description=c('european', 'african', 'american', 'asian')
    )
d <- af.filter(onekg,onekg.thresh)

#ESP
esp <- data.frame(
    pop=c('ESP_EA','ESP_AA'),
    description=c('european-african', 'african-american')
    )
d <- af.filter(esp,esp.thresh)

# all variant should have depth greater than depth.thresh
d <- d[which(rowSums(d[,c(depth.id,depth.dadid,depth.momid)]>depth.thresh)==3),]

# remove variants which are also hom in parents since they are unaffected
d <- d[which(d[,geno.id]==2&d[,geno.dadid]!=2&d[,geno.momid]!=2),]

# het in the parents: rare recessive hom
i <- which(d[,geno.id]==2 & d[,geno.dadid]==1 & d[,geno.momid]==1)
write.csv(d[i,], file=sprintf('chr%s_rare-recessive-hom.csv',opt$chr), quote=FALSE, row.names=FALSE  )

# de novo hom
write.csv(d[-i,], file=sprintf('chr%s_denovo-hom.csv',opt$chr), quote=FALSE, row.names=FALSE  )





