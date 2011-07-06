fitsaemodel.control <-
function(niter=40, iter=c(200, 200), acc=1e-5, init="default", ...){
   # define acc
   if (length(acc) != 4){
      acc = rep(acc, 4)
   }
   if (length(iter) != 2){
      iter = rep(iter[1], 2)
   }
   # implicitly check for postitivity
   acc = abs(acc)
   iter = abs(iter)
   niter = abs(niter[1])
   # define maxk (define ml method)
   maxk = 20000
   # machine eps
   eps <- .Machine$double.eps^(1/4) 
   # make them all positive
   init <- switch(init, "default"=0, "lts"=1, "s"=2) 
   res = list(niter=niter, iter=iter, acc=acc, maxk=maxk, init=init, add=list(...))
   return(res)
}

