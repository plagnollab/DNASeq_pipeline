#!/usr/bin/env Rscript

# first combine all chromosomes into one file for annotations, genotypes and genotypes_depth
for (file in c('annotations','genotypes','genotypes_depth')) {
    print(dim(d <- do.call('rbind', lapply(c(seq(1,22),'X'), function(chr) {
        print(chr)
        print(dim(x <- read.csv(sprintf('VEP_%s-%s.csv',chr,file))))
        print(dput(colnames(x)))
        return(x)
    }
    )) ))
    write.csv(d,file=sprintf('VEP-%s.csv',file),row.names=FALSE,quote=FALSE)
}


# then combine annotations, genotypes nd genotypes_depth
ann <- read.csv('VEP-annotations.csv')
geno <- read.csv('VEP-genotypes.csv')
print(dim(geno <- geno[,-which('VARIANT_ID'==colnames(geno))]))
colnames(geno) <- paste('geno',colnames(geno),sep='.')
geno.depth <- read.csv('VEP-genotypes_depth.csv')
print(dim(geno.depth <- geno.depth[,-which('VARIANT_ID'==colnames(geno.depth))]))
colnames(geno.depth) <- paste('depth',colnames(geno.depth),sep='.')

d <- cbind(ann,geno,geno.depth)

write.csv(d, quote=FALSE, row.names=FALSE, file='VEP.csv')


