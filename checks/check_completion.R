library(ggplot2)
releases <- c("July2016", "June2016")

#file.type <- "exome_table.csv"
#file.type <- "_recal_filtered2.vcf"
file.type <- ".RData"

complete <- data.frame()
for (release in releases) {

  path <- paste0("/SAN/vyplab/UCLex/mainset_", release, "/mainset_", release, "_snpStats")
  my.files <- paste0(path, "/chr", c(1:22, "X"), "_snpStats", file.type)
  loc.size <- file.info(my.files)[,1, drop = FALSE]
  loc.size$name <- row.names(loc.size)
  loc.size$release <- release
  loc.size$chromosome <- c(1:22, "X")
  
  if (nrow(loc.size) != 23) stop("Missing files ", nrow(loc.size))
  complete <- rbind.data.frame( complete, loc.size )
}

complete$size <- complete$size / 1000
complete$chromosome <- factor(complete$chromosome, levels = c(1:22, "X"))
p <- ggplot(data = complete, aes(x = chromosome, y = size, colour = release))
p <- p + geom_point()
ggsave(p, file = "check.pdf")
