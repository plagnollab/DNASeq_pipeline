
#parents
#Levine_Macrogen_1477
#Levine_Macrogen_1481
#siblings
# affected:
# Levine_Macrogen_1476
# unaffected:
# Levine_Macrogen_1479
# Levine_Macrogen_1478 not sequenced
# Levine_Macrogen_1480


unaffected.parents <- c('Levine_Macrogen_1477', 'Levine_Macrogen_1481') 
siblings <- c( 'Levine_Macrogen_1476', 'Levine_Macrogen_1478', 'Levine_Macrogen_1479', 'Levine_Macrogen_1480' )
not.sequenced <- 'Levine_Macrogen_1478'
siblings <- setdiff(siblings, not.sequenced)
aff <- 'Levine_Macrogen_1476' 
#aff <- 'Levine_Macrogen_1479' 
#aff <- 'Levine_Macrogen_1480' 
unaffected.siblings <- setdiff(siblings, aff)


variants <- data.frame()
for (chr in 1:21) {
    f <- sprintf('VEP_%s-genotypes.csv',chr)
    geno <- read.csv(f)
    i <- which(
    (geno[,aff] == 2) &
    #(geno[,aff] != geno[,unaffected.siblings[[1]]]) &
    #( (geno[,unaffected.siblings[[1]]] != geno[,unaffected.siblings[[2]]]) )  &
    (geno[,unaffected.siblings[[1]]] != 2) & (geno[,unaffected.siblings[[2]]] != 2)  & 
    (geno[,unaffected.parents[[1]]]!=2) & (geno[,unaffected.parents[[1]]]!=2) &
    (rowSums(geno[, unaffected.parents]) == 2)
    )
    variants <- rbind( variants, geno[i,] )
}

#variants[,unaffected.parents]
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

variants <- variants[order(variants$CADD,decreasing=TRUE),]

variants <- variants[which(variants$EXAC_Adj<.001),]

print(variants[,c('VARIANT_ID','Consequence','SYMBOL','CADD','EXAC_Adj',aff,unaffected.siblings,unaffected.parents,grep('depth',colnames(variants),value=T))])

#remove these columns:
#DISTANCE CADD_RAW  CADD_PHRED
variants <- variants[,-which(colnames(variants) %in% c('DISTANCE','CADD_RAW','CADD_PHRED'))]


write.csv( variants, file='', quote=FALSE, row.names=FALSE  )



