
v <- do.call('rbind', lapply(c(seq(1,22),'X'), function(chr) read.csv(sprintf('filtered-%s.csv',chr))))
v <- v[,-which(colnames(v) %in% c('CAROL','CLIN_SIG','Condel','SIFT','PolyPhen','SYMBOL_SOURCE','Existing_variation'))]

v <- v[order(v$AF,decreasing=TRUE),]

#ignore column all na
v <- v[,colSums(is.na(v))<nrow(v)]

h <- c('VARIANT_ID','ID','SYMBOL','Allele','Gene','Feature','Consequence','AF','WT','HET','HOM','MISS','R', 'ca.AF','ca.WT','ca.HET','ca.HOM','ca.MISS', 'co.AF','co.WT','co.HET','co.HOM','co.MISS', 'cDNA_position','CDS_position','Protein_position','Amino_acids','Codons','STRAND','HGNC_ID','CANONICAL','CADD','AA_MAF','EA_MAF','EXAC_AFR','EXAC_AMR','EXAC_Adj','EXAC_EAS','EXAC_FIN','EXAC_NFE','EXAC_OTH','EXAC_SAS','UCLEX')

write.csv(v[,h],file='', quote=FALSE, row.names=FALSE)






