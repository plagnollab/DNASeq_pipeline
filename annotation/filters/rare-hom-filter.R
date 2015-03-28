
variants <- data.frame()
for (chr in 1:21) {
    f <- sprintf('VEP_%s-genotypes.csv',chr)
    geno <- read.csv(f)
    variants <- rbind( variants, geno )
}

print( dim(variants) )

ann.variants <- data.frame()
for (chr in 1:21) {
    f <- sprintf('VEP_%s-annotations.csv',chr)
    ann <- read.csv(f)
    rownames(ann) <- ann$VARIANT_ID
    ann.variants <- rbind(ann.variants, ann)
}

print( dim(variants <- cbind(ann.variants[variants$VARIANT_ID,],variants)) )

depth.variants <- data.frame()
for (chr in 1:21) {
    f <- sprintf('VEP_%s-genotypes_depth.csv',chr)
    depth <- read.csv(f)
    colnames(depth) <- paste(colnames(depth),'depth',sep='.')
    rownames(depth.variants) <- depth.variants$VARIANT_ID
    depth.variants <- rbind(depth.variants, depth)
}

print( dim(variants <- cbind(variants,depth.variants[variants$VARIANT_ID,])) )

variants <- variants[,unique(colnames(variants))]

#remove these columns:
#DISTANCE CADD_RAW  CADD_PHRED
variants <- variants[,-which(colnames(variants) %in% c('DISTANCE','CADD_RAW','CADD_PHRED'))]

variants <- variants[order(variants$CADD,decreasing=TRUE),]

variants <- variants[which(variants$EXAC_Adj<.001),]

#no missingness
variants <- variants[which(variants$MISS==0),]

samples <- read.csv('samples.txt',header=FALSE)[,1] 
samples <- gsub('-','.',samples)

#less than half of the samples are wt
variants <- variants[which(variants$WT < (length(samples)/2)),]

print(dim(variants))

message('rare hom')
write.csv( variants[which(rowSums(variants[,samples]==2)==length(samples)),], file='', quote=FALSE, row.names=FALSE  )

message('rare het')
write.csv( variants[which(variants$HET==length(samples)),], file='', quote=FALSE, row.names=FALSE  )
 

message('top 100 CADD score')
write.csv( variants, file='', quote=FALSE, row.names=FALSE )





