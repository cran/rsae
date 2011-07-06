print.meanssaemodel <-
function(x, digits=4, ...){
   cat("Robustly Estimated/Predicted Area-Level Means:\n")
   all <- cbind(x$raneff, x$fixeff, x$means)
   colnames(all) <- c("raneff", "fixeff", "predicted mean")
   print.default(format(all, digits = digits), print.gap = 2, quote = FALSE)
   k <- attr(x, "robustness")
   cat(paste("(robusteness tuning constant k = ", k, ")\n", sep=""))
}

