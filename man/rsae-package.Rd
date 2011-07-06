\name{rsae-package}
\alias{rsae-package}
\alias{rsae}
\docType{package}
\title{
Robust Small Area Estimation
}
\description{
Computes robust basic unit- and area-level and predicts area-specific means
}
\details{
\tabular{ll}{
Package: \tab rsae\cr
Type: \tab Package\cr
Version: \tab 0.1-3\cr
Date: \tab 2011-07-06\cr
Suggests: \tab robustbase, nlme\cr
License: \tab GPL (>=2) | FreeBSD\cr
LazyLoad: \tab yes\cr
}

\subsection{Implemented methods}{
   \itemize{
      \item maximum likelihood (as reference)
      \item Huber-type M-estimation
      \item \emph{[S-estimation; not released, yet]}
      \item \emph{[Simple and robust, minimal-iterative estimation; not released, yet]}
   }
}

\subsection{How to}{

Data analysis with \code{rsae} involves the following steps: 
\enumerate{
   \item prepare the data/ set up the model for estimation; see \code{\link{saemodel}}
   \item fit the model by various (robust) methods; see \code{\link{fitsaemodel}}
   \item (robustly) predict the random effects and the area means; see \code{\link{robpredict}}
   }
}
}

\author{
Tobias Schoch

Maintainer: Tobias Schoch <tobias.schoch@fhnw.ch>
}
\references{
Rao, J.N.K. (2003) \emph{Small Area Estimation}, New York: John Wiley and Sons.

Richardson, A.M. and A.H. Welsh (1995) \emph{Robust restricted maximum likelihood in mixed linear model}, Biometrics 51, pp. 1429-1439.

Schoch, T. (2011) \emph{The robust basic unit-level small area model. A simple and fast algorithm for large datasets}, in: Proceedings of the New Technologies and Techniques Conference (NTTS), EUROSTAT, Brussels.

Sinha, S.K. and J.N.K. Rao (2009) \emph{Robust small area estimation}, Canadian Journal of Statistics 37, pp. 381-399.
}

