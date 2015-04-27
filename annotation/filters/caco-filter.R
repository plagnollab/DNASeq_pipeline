# pedigree
message('dim of pedigree')
err.cat(nrow(pedigree <- read.csv(pedigree)))
sample.affection <- pedigree[ which( pedigree$uclex.sample %in% samples ), c('uclex.sample','Affection')]
# cases
message('cases')
err.cat(cases <- sample.affection[which(sample.affection$Affection==2),'uclex.sample'])
message('number of cases')
err.cat(length(cases))
# controls
message('controls')
err.cat(controls <- sample.affection[which(sample.affection$Affection==1),'uclex.sample'])
message('number of controls')
err.cat(length(controls))


