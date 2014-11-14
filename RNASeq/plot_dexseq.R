
stop()
drawPlot <- function (matr, ylimn, ecs, intervals, rango, fitExpToVar, numexons,
                      textAxis, rt, color, colorlines, ...)  {
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, max(matr)))
  print(matr)
  makevstaxis(1/ncol(matr), ylimn, ecs, ...)
  intervals <- (0:nrow(matr))/nrow(matr)
  middle <- apply(cbind(intervals[rango], (intervals[rango +
                                                     1] - ((intervals[rango + 1]) - intervals[rango]) * 0.2)),
                  1, median)
  matr <- rbind(matr, NA)
  j <- seq <- len(ncol(matr))
  segments(intervals[rango], matr[rango, j], intervals[rango +
                                                       1] - ((intervals[rango + 1] - intervals[rango]) * 0.2),
           matr[rango, j], col = color, ...)
  segments(intervals[rango + 1] - ((intervals[rango + 1] -
                                    intervals[rango]) * 0.2), matr[rango, j], intervals[rango +
                                                                                        1], matr[rango + 1, j], col = color, lty = "dotted",
           ...)
  abline(v = middle[rango], lty = "dotted", col = colorlines)
  mtext(textAxis, side = 2, adj = 0.5, line = 1.5, outer = FALSE,
        ...)
  axis(1, at = middle[seq(along = rt)], labels = featureData(ecs)$exonID[rt],
       ...)
}
stop()

plotDEXSeq <- function (ecs, geneID, FDR = 0.1, fitExpToVar = "condition",
                        norCounts = FALSE, expression = TRUE, splicing = FALSE, displayTranscripts = FALSE,
                        names = FALSE, legend = FALSE, color = NULL, color.samples = NULL,
                        ...)
  {
    stopifnot(is(ecs, "ExonCountSet"))
    if (any(is.na(sizeFactors(ecs)))) {
      stop("Please estimate sizeFactors first\n")
    }
    if (!fitExpToVar %in% ecs@designColumns) {
      stop("fitExpToVar parameter is not in the design columns, double check ecs@designColumns")
    }
    if (sum(is.na(featureData(ecs)$dispersion)) == nrow(counts(ecs))) {
      stop("No dispersion parameters found, first call function fitDispersionFunction...\n")
    }
    op <- sum(c(expression, splicing, norCounts))
    if (op == 0) {
      stop("Please indicate what would you like to plot\n")
    }
    rt <- which(featureData(ecs)$geneID == geneID)
    count <- countTableForGene(ecs, geneID, normalized = TRUE)
    if (sum(count) == 0) {
      warning("No read counts falling in this gene, there is nothing to plot.")
      return()
    }
    if (FDR > 1 | FDR < 0) {
      stop("FDR has to be a numeric value between 0 - 1")
    }
    rango <- seq(along = rt)
    intervals <- (0:nrow(count))/nrow(count)
    numcond <- length(unique(design(ecs, drop = FALSE)[[fitExpToVar]]))
    numexons <- nrow(count)
    each <- featureData(ecs)$padjust[rt]
    exoncol <- ifelse(each <= FDR, "#F219ED", "#CCCCCC")
    exoncol[is.na(exoncol)] <- "white"
    colorlines <- ifelse(each <= FDR, "#F219ED60", "#B3B3B360")
    colorlines[is.na(colorlines)] <- "#B3B3B360"
    colorlinesB <- ifelse(each <= FDR, "#9E109B", "#666666")
    colorlinesB[is.na(colorlinesB)] <- "#666666"
    if (!any(is.na(featureData(ecs)$start))) {
      sub <- data.frame(start = featureData(ecs)$start[rt],
                        end = featureData(ecs)$end[rt], chr = featureData(ecs)$chr[rt],
                        strand = featureData(ecs)$strand[rt])
      rel <- (data.frame(sub$start, sub$end)) - min(sub$start)
      rel <- rel/max(rel[, 2])
      if (displayTranscripts == TRUE & !is.null(featureData(ecs)$transcripts)) {
        transcripts <- sapply(featureData(ecs)$transcripts[rt],
                              function(x) {
                                strsplit(x, ";")
                              })
        trans <- Reduce(union, transcripts)
        if (length(trans) > 40) {
          warning("This gene contains more than 40 transcripts annotated, only the first 40 will be plotted\n")
        }
        mat <- seq <- DEXSeq::len(3 + min(length(trans), 40))
        hei <- c(8, 1, 1.5, rep(1.5, min(length(trans), 40)))
      }
      else {
        mat <- 1:3
        hei <- c(5, 1, 1.5)
      }
      if (op > 1) {
        hei <- c(rep(hei[1], op - 1), hei)
        mat <- c(mat, length(mat) + seq(along = op))
      }
      hei <- c(hei, 0.2)
      mat <- c(mat, length(mat) + 1)
      layout(matrix(mat), heights = hei)
      par(mar = c(2, 4, 4, 2))
    }
    else if (op > 1) {
      par(mfrow = c(op, 1))
    }
    if (is.null(color)) {
      color <- rgb(colorRamp(c("#D7191C", "#FFFFBF", "#2B83BA"))(seq(0,
                                                                     1, length.out = numcond)), maxColorValue = 255, alpha = 175)
    }
    names(color) <- sort(levels(design(ecs, drop = FALSE)[[fitExpToVar]]))
    if (expression) {
      es <- DEXSeq:::fitAndArrangeCoefs(ecs, geneID, frm = as.formula(paste("count ~",
                                              fitExpToVar, "* exon")))
      if (is.null(es)) {
        warning(sprintf("glm fit failed for gene %s", geneID))
        return()
      }
      coeff <- as.matrix(t(DEXSeq:::getEffectsForPlotting(es, averageOutExpression = FALSE,
                                                          groupingVar = fitExpToVar)))
      coeff <- exp(coeff)
      ylimn <- c(0, max(coeff, na.rm = TRUE))
      coeff <- DEXSeq:::vst(coeff, ecs)
      DEXSeq:::drawPlot(matr = coeff, ylimn, ecs, intervals, rango,
               textAxis = "Fitted expression", rt = rt, color = rep(color[sort(levels(design(ecs,
                   drop = FALSE)[[fitExpToVar]]))], each = numexons),
               colorlines = colorlines, ...)
    }
    if (splicing) {
      coeff <- as.matrix(t(DEXSeq:::getEffectsForPlotting(fitAndArrangeCoefs(ecs,
                                                                    geneID, frm = as.formula(paste("count ~", fitExpToVar,
                                                                              "* exon"))), averageOutExpression = TRUE, groupingVar = fitExpToVar)))
      coeff <- exp(coeff)
      ylimn <- c(0, max(coeff, na.rm = TRUE))
      coeff <- DEXSeq:::vst(coeff, ecs)
      DEXSeq:::drawPlot(matr = coeff, ylimn, ecs, intervals, rango,
               textAxis = "Fitted splicing", rt = rt, color = rep(color[sort(levels(design(ecs,
                                                        drop = FALSE)[[fitExpToVar]]))], each = numexons),
               colorlines = colorlines, ...)
    }
    if (norCounts) {
      ylimn <- c(0, max(count, na.rm = TRUE))
      count <- DEXSeq:::vst(count, ecs)
      if (is.null(color.samples)) {
        colorcounts <- rep(color[as.character(design(ecs,
                                                     drop = FALSE)[[fitExpToVar]])], each = numexons)
      }
      else {
        colorcounts <- rep(color.samples, each = numexons)
      }
      DEXSeq:::drawPlot(matr = count, ylimn, ecs, intervals, rango,
               textAxis = "Normalized counts", rt = rt, color = colorcounts,
               colorlines = colorlines, ...)
    }
    print(featureData(ecs)$start)
    if (!any(is.na(featureData(ecs)$start))) {
      par(mar = c(0, 4, 0, 2))
      plot.new()
      segments(apply((rbind(rel[rango, 2], rel[rango, 1])),
                     2, median), 0, apply(rbind(intervals[rango], intervals[rango +
                                                                            1] - ((intervals[rango + 1] - intervals[rango]) *
                                                                                  0.2)), 2, median), 1, col = colorlinesB)
      par(mar = c(1.5, 4, 0, 2))
      DEXSeq:::drawGene(min(sub$start), max(sub$end), tr = sub, rango,
               exoncol = exoncol, names, trName = "Gene model",
               cex = 0.8)
      if (!is.null(featureData(ecs)$transcripts)) {
        print("tada2")
        if (displayTranscripts) {
          print("tada")
          for (i in seq <- DEXSeq::len(min(length(trans), 40))) {
            logicexons <- sapply(transcripts, function(x) {
              length(which(x == trans[i]))
            })
            tr <- data.frame(featureData(ecs)$start[rt][logicexons ==
                                                        1], featureData(ecs)$end[rt][logicexons ==
                                                                                     1])
            DEXSeq:::drawGene(min(sub$start), max(sub$end), tr = tr,
                     rango, exoncol = NULL, names, trName = trans[i],
                     cex = 0.8)
          }
        }
      }
      axis(1, at = round(seq(min(sub$start), max(sub$end),
                length.out = 10)), labels = round(seq(min(sub$start),
                                     max(sub$end), length.out = 10)), pos = 0, lwd.ticks = 0.2,
           padj = -0.7, ...)
    }
    if (legend) {
      mtext(paste(geneID, unique(featureData(ecs)$strand[rt])),
            side = 3, adj = 0.25, padj = 1.5, line = 0, outer = TRUE,
            cex = 1.5)
      posforlegend <- seq(0.7, 0.9, length.out = numcond)
      for (i in seq(along = color)) {
        mtext(names(color[i]), side = 3, adj = posforlegend[i],
              padj = 1.5, line = 0, outer = TRUE, col = color[i],
              ...)
      }
    }
    else {
      mtext(paste(geneID, unique(featureData(ecs)$strand[rt])),
            side = 3, adj = 0.5, padj = 1.5, line = 0, outer = TRUE,
            cex = 1.5)
    }
  }
