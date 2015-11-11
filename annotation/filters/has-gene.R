#!/usr/bin/env Rscript
err.cat <- function(x)  cat(x, '\n', file=stderr())

message('*** coding FILTERING ***')
d <- read.csv(file('stdin'))

d <- d[which(!is.na(d$SYMBOL)),]

write.csv(d, file='', quote=FALSE, row.names=FALSE)

