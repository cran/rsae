plot.meanssaemodel <-
function(x, y=NULL, sort=NULL, ...){
   # y is part of the generic (later, we will use y to allow plot comparison)
   k <- attr(x, "robustness")
   fe <- x$fixeff
   re <- x$raneff
   means <- x$means
   g <- length(re)
   areaNames <- rownames(re)
   # sorting
   if (!is.null(sort)){
      ord <- switch(sort, ranef=order(re), fixef=order(fe), means=order(means)) 
      re <- re[ord]
      fe <- fe[ord]
      means <- means[ord]
      areaNames <- areaNames[ord]
   }
   ra <- range(means)
   ra[1] <- min(fe) 
   # add an null-line to the plot (for the legend) 
   g <- g + 1
   at <- 1:g
   re <- c(re, NA)
   fe <- c(fe, NA)
   means <- c(means, NA)  
   areaNames <- c(areaNames, "") 
   # prepare the plot
   op <- par(mfcol=c(1,1), mar=c(4,8,2,4))
   plot(means, at, type="b", col=2, lwd=2, axes=FALSE, xlab="predicted mean", ylab="", xlim=ra, main="Predicted means")
   lines(fe, at, type="b", col=1, lwd=2, xlim=ra)
   axis(2, seq(1, g), labels=areaNames, las=1)
   axis(1)
   grid(col="gray65", lty=2, lwd=1)
   box()
   legend("top", pch=c(1,1), lty=c(1,1), col=c(1,2), legend=c("fixeff prediction", "full prediction"), bg="white", ncol=2) 
}

