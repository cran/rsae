.computekappa <-
function(k){
   2*(k^2*(1-pnorm(k)) + pnorm(k) - 0.5  - k*dnorm(k))
}

.fitsaemodel.huberm <-
function(method, model, k, control=fitsaemodel.control(...), ...){
   #machine eps
   eps <- .Machine$double.eps
   #switch method
   if (method == "ml"){
      k <- control$maxk
      kappa <- 1
      methodName <- list(type="Maximum likelihood estimation")
   }else{
      #check whether k exists 
      if (missing(k)) stop("Robustness tuning constant is missing! \n")
      #check whether the range of k makes sense
      if (k < 0) stop("Robustness tuning constant k must be > 0.\n")
      kmin <- eps^(1/4) 
      if (k < kmin) stop(paste("Robustness tuning constant k is too small (note kmin = ", kmin, ").\n"))
      #consistency correction term for Huber-Proposal-2 scale estimate
      kappa <- .computekappa(k)
      methodName <- list(type="Huber-type M-estimation", tuning=list(k=k))
   } 
   #data preparation
   x <- model$X
   y <- model$y
   #vector of area size
   nsize <- model$nsize
   #number of observations
   n <- model$n
   #number of areas
   g <- model$g
   #number of fixed effects
   p <- model$p
   #total number of overall iterations
   niter <- control$niter
   #vector of iterations delivered drsaebetaiter and drsaehubdest
   iter <- control$iter
   #matrix recording the number of iterations
   iterrecord <- matrix(0, niter, 3)
   #acc specification (note that we take abs() to ensure that it is not negative)
   allacc <- control$acc[1]
   acc <- control$acc[2:4]
   #sum of Huber downweighting
   sumwgt <- c(0, 0, 0)
   #initialize the full parameter vector
   if(length(control$add) == 0){
      init <- .initmethod(model, control$init)
   }else{
      init <- .initmethod(model, control$init, control$add)
   }
   tau <- init
   #matrix recodring iteration-specific estimates
   taurecord <- matrix(0, niter, (p+2))
   #define epsd (minimal d that is different from zero)
   epsd <- eps^(1/4)
   #call
   tmp <- .Fortran("drsaehub", n=as.integer(n), p=as.integer(p), g=as.integer(g), niter=as.integer(niter), nsize=as.integer(nsize), iter=as.integer(iter), iterrecord=as.matrix(iterrecord), allacc=as.double(allacc), acc=as.matrix(acc), sumwgt=as.matrix(sumwgt), xmat=as.matrix(x), yvec=as.matrix(y), k=as.double(k), kappa=as.double(kappa), epsd=as.double(epsd), tau=as.matrix(tau), taurecord=as.matrix(taurecord), converged=as.integer(0))
   #return values
#HOTFIX!, check for cycling an choose the parameter-vector estimate whose estimate of v is closer to the (robust) init (i.e., either the "lts" or "s" estimate) value. This method is not supported for init=default or ml
   converged <- tmp$converged
   if (control$init > 0 & converged == 0){
      taurecord <- tmp$taurecord
      # take the second diff, take the mean over each parameter vector, and take the last that fullfiled the criterion
      u <- max(which(rowMeans(diff(taurecord, 2)) <= eps^(1/2)))
      # take the difference from init for at u
      uat <- abs(taurecord[u, (p+1)] - init[(p+1)])
      # take the difference from init for one before u
      ubefore <- abs(taurecord[(u-1), (p+1)] - init[(p+1)])
      if (uat <= ubefore){
	 tau <- taurecord[u, ]
	 converged <- 1
      }else{
	 tau <- taurecord[(u-1), ]
	 converged <- 1
      }
   }else{
      tau <- tmp$tau
   }
   #compute the covariance matrix of the fixed effects
   if (converged == 1){
      vcovbeta <- matrix(0, p, p)
      o = .Fortran("drsaehubvariance", n=as.integer(n), p=as.integer(p), g=as.integer(g), nsize=as.integer(nsize), kappa=as.double(kappa), v=as.double(tau[p+1]), d=as.double(tau[p+2]), xmat=as.matrix(x), vcovbeta=as.matrix(vcovbeta))
      vcovbeta <- o$vcovbeta
   }else{
      vcovbeta <- NULL
   }
   res <- list(beta=tau[1:p], theta=c(tau[p+1], tau[p+1]*tau[p+2]), converged=converged, vcovbeta=vcovbeta)
   #additional attributes
   attr(res, "call") <- match.call()
   attr(res, "optim") <- list(acc=c(allacc, acc), niter=c(niter, iter), usediter=tmp$iterrecord, tau=tmp$taurecord)
   if (method == "huberm"){
      attr(res, "robustness") <- list(wgt=tmp$sumwgt)
   }
   attr(res, "init") <- init
   attr(res, "method") <- methodName 
   attr(res, "saemodel") <- model
   class(res) <- "fitsaemodel" 
   return(res)  
}

.initmethod <-
function(model, init, ...){
   n <- model$n
   p <- model$p
   intercept <- model$intercept
   # default
   if (init == 0){
      result <- rep(1, (p+2))
   }
   # check wheter robustbase must be loaded
   if (init > 0){
      checkrobustbase <- require(robustbase)
      if(!checkrobustbase) stop("You cannot use either 'lts' or 's', because the \npackage 'robustbase' is not installed! \n")
   }
   # lts
   if (init == 1){
      x <- as.matrix(model$X)
      #check if it has an intercept
      if (intercept == 1){
	 x <- as.matrix(x[, (2:p)])
	 intercept <- TRUE
      }
      else{
	 cat
	 intercept <- FALSE
      }
      y <- model$y
      tmp <- robustbase::ltsReg(x=x, y=y, intercept=intercept, ...)
      # compute (robust) variance bound of d
      # repare return value (beta, v, d)
      result <- as.numeric(c(tmp$coefficients, tmp$raw.scale^2, 1))
   }
   # lmrob.S
   if (init == 2){
      x <- as.matrix(model$X)
      y <- model$y
      control <- robustbase::lmrob.control(...)
      tmp <- robustbase::lmrob.S(x=x, y=y, control=control)
      result <- as.numeric(c(tmp$coefficients, tmp$scale^2, 1))
   }
   return(result)
}

.mspe <-
function(fit, reps, areameans, fixeff){
   #fitted-model attrs 
   #as scale not variance, since rnorm wants it like that
   theta <- sqrt(fit$theta)
   beta <- fit$beta
   #model attrs
   model <- attr(fit, "saemodel")
   X <- as.matrix(model$X)
   Xbeta <- X %*% beta
   n <- model$n
   g <- model$g
   nsize <- model$nsize
   #
   predicts <- matrix(NA, reps, g)
   # 
   for (j in 1:reps){
      #draw model error, e
      e <- rnorm(n, 0, theta[1])
      #draw raneff, v
      v <- unlist(sapply(nsize, function(u) rep(rnorm(1, 0, theta[2]), u), simplify=TRUE))
      #modify pred (add random effec)
      predrf <- fixeff +  (unique(v))# * as.double(nsize))
      #generate bootstrap samples (and fix it to the model)
      model$y <- Xbeta + e + v 
      #compute the model parameters using ml
      tmp <- fitsaemodel("ml", model)
      #predict 
      predicts[j, ] <- t(robpredict(tmp, areameans, k=20000, reps=NULL)$means) - t(predrf)
   }
   res <- colMeans(predicts^2)
   return(res)
}

