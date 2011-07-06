\name{robpredict}
\alias{robpredict}
\alias{print.meanssaemodel}
\alias{plot.meanssaemodel}
\alias{residuals.meanssaemodel}

\title{
Robust prediction of random effects, fixed effects, and area-specific means
}
\description{
Robust prediction of random effects, fixed effects, and area-specific means. It can predict based on new, directly delivered, areadata. 
}
\usage{
robpredict(fit, areameans=NULL, k=NULL, reps=NULL)

\method{print}{meanssaemodel}(x, digits=4, ...)
\method{plot}{meanssaemodel}(x, y=NULL, type="e", sort=NULL, ...)
\method{residuals}{meanssaemodel}(object, ...)
}

\arguments{
  \item{fit}{
a fitted SAE model; object of class \code{fitsaemodel}
}
   \item{areameans}{
   numeric matrix (typically, with area-level means); the no. of rows must be equal to the no. of areas; the no. of columns must be equal to the no. of fixed-effects coefficients (incl. intercept). By default: \code{areadata=NULL}, i.e., predictions are based on those data that have been used to estimate the model.
}
  \item{k}{
robustness tuning constant (of the Huber psi-function) for robust prediction. Notice that \code{k} does not necessarily be the same as the \code{k} that has been used in \code{fitsaemodel}. By default, \code{k} is equal to the tuning constant used in estimating the model parameters. 
}
  \item{reps}{
no. of bootstrap replicates for mean squared prediction error; default: \code{reps=NULL}
}

   \item{x}{
object of the class \code{"meanssaemodel"}; this argument is only used in the \code{print} method.
}
   \item{digits}{
integer, defining the number of decimal places to be shown in the \code{print} method (default: \code{digits=4})
}
   \item{y}{
has no meaning, yet! (default: \code{y=NULL}; needs to included in the args list, because it is part of plot's generic arg definition)
}
   \item{type}{
character specifying the \code{plot} method; either \code{"e"} (error bars; default) or \code{"l"} (lines).
}
   \item{sort}{
only used in the \code{plot} method; if \code{sort="means"}, the predicted means are ploted in ascending order (default: \code{sort=NULL}); similarly, with \code{sort="fixef"} and \code{sort="ranef"} the predicted means are sorted along the fixed effects or the random effects, respectively
}
   \item{object}{
object of the class \code{fitsaemodel}; a fitted model used in the \code{residuals} method.
}
   \item{\dots}{
not used
}
}
\details{
Given the robustly estimated SAE model, one considers robustly predicting the random- and fixed effect (and the final area-specific means). The tuning constant \code{k} regulates the degree of robustness when predicting the random effects. 
}
\value{
Instance of the S3 class \code{meanssaemodel}
}
\references{
Schoch, T. (2011) \emph{The robust basic unit-level small area model. A simple and fast algorithm for large datasets}, in: Proceedings of the New Technologies and Techniques Conference (NTTS), EUROSTAT, Brussels.
}
\author{
Tobias Schoch
} 
\examples{
#generate the synthetic data/model
mymodel <- makedata()
#compute Huber M-estimation type estimates of the model "mymodel"
#robustness tuning constant k = 2
myfittedmodel <- fitsaemodel("huberm", mymodel, k=2)
myfittedmodel
#get a summary of the model
summary(myfittedmodel)
#robustly predict the random effects and the area-level means.  
#Here, we choose the robustness tuning constant k equal to 1.8
mypredictions <- robpredict(myfittedmodel, k=1.8)
mypredictions
}

