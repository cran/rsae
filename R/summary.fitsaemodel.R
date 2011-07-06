summary.fitsaemodel <-
function (object, full=FALSE, digits=3, ...){
   saemodel <- attr(object, "saemodel")
   n <- saemodel$n
   # retrieve the estimating method
   method <- attr(object, "method") 
   cat("SUMMARY: ESTIMATES OF SAE-MODEL (model type B) \n")
   # check whether the model converged
   converged <- object$converged
   if (converged != 1){
      cat("NOTE: ALGORITHM DID NOT CONVERGE! (see acc and niter, below)\n")
   }
   cat("Method: ", method$type, "\n")
   # branch: robust vs. non-robust methods
   if (length(method) > 1){
      tuning <- method$tuning
      if (length(tuning) == 1){
	 cat(paste("Robustness tuning constant: ",names(tuning), " = ", as.numeric(tuning), "\n"))
      }else{
	 for (i in 1:length(tuning)){
	    cat(tuning[i], "\n")  
	 }
      }
   }
   #----------------------
   # robustness properties
   robustness <- attr(object, "robustness") 
   # branch robust vs non-robust
   if (!is.null(robustness)){
      wgt <- t(robustness$wgt) / n
      colnames(wgt) <- c("fixeff", "residual var", "area raneff var")
      rownames(wgt) <- "sum(wgt)/n"
      cat("---\n")
      cat("Degree of downweighting/winsorization:\n")
      cat("\n")
      print.default(format(t(wgt), digits = digits), print.gap = 2, quote = FALSE)
   }
   #----------------------
   # niter and acc specification
   optim <- attr(object, "optim")
   acc <- optim$acc
   niter <- c(optim$niter, 100) # 100 is the max value defined in estimation of "d"
   together <- cbind(niter, acc)
   colnames(together) <- c("niter", "acc")
   rownames(together) <- c("overall loop", "fixeff", "residual var", "area raneff var")
   cat("---\n")
   cat("User specified number of iterations (niter) and \nnumeric precision (acc):\n")
   cat("\n")
   print.default(format(together, digits = 1), print.gap = 2, quote = FALSE)
   #----------------------
   # used iters
   iters <- optim$usediter
   if (dim(iters)[1] == 1){
      iters <- as.matrix(iters)
   }else{
      #remove the zero entires in iters
      iters <- iters[rowSums(iters) != 0 ,]
   }
   colnames(iters) <- c("fixeff", "residual var", "area raneff var")
   rownames(iters) <- paste(seq(1, dim(iters)[1]))
   cat("---\n")
   cat(paste("Number of runned EE-specific iterations in each \ncall (given the user-defined specs), reported for \neach of the ",dim(iters)[1], "overall iterations separately: \n"))
   cat("\n")
   print.default(format(iters, digits = 1), print.gap = 2, quote = FALSE)
   #----------------------
   # if not converged or if full=TRUE
   if (converged != 1 | full == TRUE){
      tau <- as.matrix(optim$tau)
      if (dim(tau)[1] > 1){
	 tau <- tau[rowSums(tau) != 0 ,]
      }
      colnames(tau) <- c(colnames(saemodel$X), "residual var", "area raneff var") 
      rownames(tau) <- paste(seq(1, dim(tau)[1]))
      cat("Convergence progress (estimates at each overall iteration):\n")
      print.default(format(tau, digits = digits), print.gap = 2, quote = FALSE)
   }
   cat("-EOF-\n")
}

