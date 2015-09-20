
# Loss of heterozygosity (LOH) is a gross chromosomal event that results in loss of the entire gene and the surrounding chromosomal region.
# Most of the chromosomes within somatic cells of individuals are paired, allowing for SNP locations to be potentially heterozygous.
# However, one parental copy of a region can sometimes be lost, which results in the region having just one copy.
# The single copy cannot be heterozygous at SNP locations and therefore the region shows loss of heterozygosity (LOH).
# Loss of heterozygosity due to loss of one parental copy in a region is also called hemizygosity in that region.


ped <- read.table('pedigree_details.ped',col.names=c('fid','id','dadid','momid','gender','affection'))

# combined file containing annotations, genotypes and genotype.depths
d <- read.csv('VEP.csv')
rownames(d) <- d$VARIANT_ID

d$chr <- unlist(lapply(strsplit(d$VARIANT_ID,'_'),`[[`,1))
d$pos <- as.numeric(lapply(strsplit(d$VARIANT_ID,'_'),`[[`,2))

geno.id <- paste('geno',ped[1,'id'],sep='.')
geno.momid <- paste('geno',ped[1,'momid'],sep='.')
geno.dadid <- paste('geno',ped[1,'dadid'],sep='.')

depth.id <- paste('depth',ped[1,'id'],sep='.')
depth.momid <- paste('depth',ped[1,'momid'],sep='.')
depth.dadid <- paste('depth',ped[1,'dadid'],sep='.')

dd <- subset(d,chr=='22')

pdf('~/chr22.pdf')
plot(dd[,c('pos',geno.id)],col='black',pch=20,type='l')
#lines(dd[,c('pos',geno.dadid)],col='red',pch=20)
#lines(dd[,c('pos',gemo.momid)],col='green',pch=20)
dev.off()

dd <- subset(dd,dd[,depth.id]>20)

# Run length method
loh <- rle(dd[,id]==2)
loh <- data.frame(values=loh$values,lengths=loh$lengths)
loh$start <- c(1,cumsum(loh$lengths)[-length(loh$lengths)])
loh$end <- loh$start+loh$lengths-1
loh$bp.start <- dd$pos[loh$start]
loh$bp.end <- dd$pos[loh$end]
loh$bp.length <- loh$bp.end-loh$bp.start
# regions which has the most consecutive hom
loh <- subset(loh,loh$values)
loh[which.max(loh$lengths),]

# Sliding window method (very slow)
window.width <- 1000
hom.rate <- data.frame()
for (i in min(dd$pos):max(dd$pos)) {
 j <- which(i <= dd$pos & dd$pos <= i+window.width)
 if (length(j)==0) next
 hom.count <- prop.table(table(dd[j,geno.id]==2))
 if (!'TRUE' %in% names(hom.count)) next
 hom.rate <- rbind(hom.rate,c(start=i,end=i+window.width,rate=hom.count[['TRUE']]))
}




