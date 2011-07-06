robpredict <-
function(fit, k, areadata=NULL){
   if (!inherits(fit, "fitsaemodel")) stop("fit must be of class 'fitsaemodel'")
   if (k <= 0 ) stop("Robustness tuning constant k must be > 0!\n")
   # set k=inf to a numerically valid value
   if (k == Inf){
      k <- 20000
      kappa <- 1
   }else{
      kappa <- .computekappa(k)
   }
   # from the model definitions; used in order to compute the random effects
   model <- attr(fit, "saemodel")
   areaNames <- attr(model, "areaNames")
   x <- model$X
   y <- model$y
   n <- model$n
   p <- model$p
   g <- model$g
   nsize <- model$nsize
   # from the fitted model
   modelk <- attr(fit, "method")$tuning$k
   if (is.null(modelk)){
      modelk <- 20000
   }
   beta <- fit$beta
   v <- fit$theta[1]
   d <- fit$theta[2] / v
   # tau
   tau <- c(beta, v, d)
   #
   predre <- rep(0, g)
   predfe <- rep(0, g)
   tmp <- .Fortran("drsaehubpredict", n=as.integer(n), p=as.integer(p), g=as.integer(g), nsize=as.integer(nsize), k=as.double(k), kappa=as.double(kappa), d=as.double(d), v=as.double(v), beta=as.matrix(beta), yvec=as.matrix(y), xmat=as.matrix(x), predfe=as.matrix(predfe), predre=as.matrix(predre))
   # retrieve the area-level random effects; it is used whether new data is present or not
   raneff <- tmp$predre
   # branch: old vs new data
   if (is.null(areadata)){
      fixeff <- as.matrix(tmp$predfe)
   }else{
      # check whether the new data are proper
      if(!is.matrix(areadata)){
	 areadata <- as.matrix(areadata)
      }
      # check the dimensions
      if (dim(areadata)[1] != g) stop("'areadata' is not of conformable size! \n")
      if (dim(areadata)[2] != p) stop("'areadata' is not of conformable size! \n")
      # compute the fixed-effect spredictions (at the area level)
      fixeff <- areadata %*% beta
   }
   means <- raneff + fixeff
   rownames(fixeff) <- areaNames
   rownames(raneff) <- areaNames
   rownames(means) <- areaNames
   # compute the residuals of the model (i.e. e_ij = y_ij - X_ij*beta - u_i)
   vn <- numeric(n)
   getres <- .Fortran("drsaeresid", n=as.integer(n), p=as.integer(p), g=as.integer(g), nsize=as.integer(nsize), k=as.double(modelk), tau=as.matrix(tau), u=as.matrix(raneff), xmat=as.matrix(x), yvec=as.matrix(y), res=as.matrix(vn), stdres=as.matrix(vn), wgt=as.matrix(vn))
   #
   result <- list(fixeff=fixeff, raneff=raneff, means=means, res=getres$res, stdres=getres$stdres, wgt=getres$wgt)
   attr(result, "robustness") <- k
   attr(result, "fit") <- fit
   class(result) <- "meanssaemodel"
   return(result)
}

