
v <- do.call('rbind', lapply(c(seq(1,22),'X'), function(chr) read.csv(sprintf('filtered-%s.csv',chr))))
v <- v[,-which(colnames(v) %in% c('CAROL','CLIN_SIG','Condel','SIFT','PolyPhen','SYMBOL_SOURCE','Existing_variation'))]

v <- v[order(v$AF,decreasing=TRUE),]

#ignore column all na
v <- v[,colSums(is.na(v))<nrow(v)]

write.csv(v,file='')






