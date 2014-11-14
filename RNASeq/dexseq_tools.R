
plotDispEsts <- function( cds, ymin, linecol="#ff000080", xlab = "mean of normalized counts", ylab = "dispersion", log = "xy", cex = 0.45, ... ) {
      px = rowMeans( counts( cds, normalized=TRUE ) )
      sel = (px>0)
      px = px[sel]
      
      py = fData(cds)$dispBeforeSharing[sel]
      if(missing(ymin)) ymin = 10^floor(log10(min(py[py>0], na.rm=TRUE))-0.1)
      
      plot(px, pmax(py, ymin), xlab=xlab, ylab=ylab,
           log=log, pch=ifelse(py<ymin, 6, 16), cex=cex, ... )
      xg = 10^seq( -.5, 5, length.out=100 )
      fun = function(x) { cds@dispFitCoefs[1] + cds@dispFitCoefs[2] / x }
      lines( xg, fun(xg), col=linecol, lwd=4)
    }

plotMA <- function(x, ylim,col = ifelse(x$padj>=0.1, "gray32", "red3"), linecol = "#ff000080", xlab = "mean of normalized counts", ylab = expression(log[2]~fold~change), log = "x", cex=0.45, ...) {
  if (!(is.data.frame(x) && all(c("baseMean", "log2FoldChange") %in% colnames(x))))
    stop("'x' must be a data frame with columns named 'baseMean', 'log2FoldChange'.")
  
  x = subset(x, baseMean!=0)
  py = x$log2FoldChange
  if(missing(ylim))
    ylim = c(-1,1) * quantile(abs(py[is.finite(py)]), probs=0.99) * 1.1
  plot(x$baseMean, pmax(ylim[1], pmin(ylim[2], py)),
       log=log, pch=ifelse(py<ylim[1], 6, ifelse(py>ylim[2], 2, 16)),
       cex=cex, col=col, xlab=xlab, ylab=ylab, ylim=ylim, ...)
  abline(h=0, lwd=4, col=linecol)
}


fitDispersionFunction <- function (ecs) {
  stopifnot(is(ecs, "ExonCountSet"))
  if (all(is.na(fData(ecs)$dispBeforeSharing))) {
    stop("no CR dispersion estimations found, please first call estimateDispersions function")
  }
  means <- colMeans(t(counts(ecs))/sizeFactors(ecs))
  disps <- fData(ecs)$dispBeforeSharing
  coefs <- c(0.1, 1)
  iter <- 0
  while (TRUE) {
    residuals <- disps/(coefs[1] + coefs[2]/means)
    good <- which((residuals > 1e-04) & (residuals < 15) )
    mm <- model.matrix(disps[good] ~ I(1/means[good]))
    fit <- try(statmod:::glmgam.fit(mm, disps[good], coef.start = coefs),
               silent = TRUE)
    if (inherits(fit, "try-error")) {
      stop("Failed to fit the dispersion function\n")
    }
    oldcoefs <- coefs
    coefs <- coefficients(fit)

    if (coefs[1] < 0) {
      coefs[1] <- 0.001
      warning("Negative intercept value in the dispersion function, it will be set to 0. Check fit diagnostics plot section from the vignette.")
      #break
    }
    
    if (sum(log(coefs/oldcoefs)^2) < 0.005)
      break
    iter <- iter + 1
    if (iter > 10) {
      warning("Dispersion fit did not converge.")
      break
    }
  }
  ecs@dispFitCoefs <- coefs
  fData(ecs)$dispFitted <- ecs@dispFitCoefs[1] + ecs@dispFitCoefs[2]/colMeans(t(counts(ecs))/sizeFactors(ecs))
  fData(ecs)$dispersion <- pmin(pmax(fData(ecs)$dispBeforeSharing,
                                     fData(ecs)$dispFitted, na.rm = TRUE), 1e+08)
  return(ecs)
}






