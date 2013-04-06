\encoding{latin1}
\name{autoplsVAL}
\alias{autoplsVAL}
\alias{R2.autopls}
\alias{RMSEP.autopls}
\alias{jack.test.autopls}
\alias{metaval}
\alias{repeatedCV}
\title{
  Validate a fitted autopls model
}
\description{
  Functions to extract R2 and RMSEP from \code{autopls} objects and for 
  significance testing based on jackknife variance estimates for regression 
  coefficients.
}
\usage{
  \method{R2}{autopls}(object, estimate, nc = 'inherit', ic = FALSE, \dots)
  \method{RMSEP}{autopls}(object, estimate, nc = 'inherit', ic = FALSE, \dots)
  jack.test.autopls (object, nc = 'inherit')
  metaval (object, method, estimate, ic)
  repeatedCV (object, k = 100, segments = 4)
}
\arguments{
  \item{object}{object of class \code{autopls}}
  \item{method}{character. Should be or 'R2' or 'RMSEP'}
  \item{estimate}{character vector. Which estimators to use. In \code{metaval}
    this can be \dQuote{train} or \dQuote{CV}. Additional options in \code{R2} 
    and \code{RMSEP} are \dQuote{all} and \dQuote{test}).} 
  \item{nc}{\code{'inherit'} returns values corresponding to the number of 
    latent vectors in the current model, \code{'all'} returns values for all
    numbers of latent vectors. A specific number returns values corresponding to 
    the respective number of latent vectors.}
  \item{ic}{logical. Specifies whether estimates for a model with zero 
    components should be returned}
  \item{k}{number of cross-validations used in \code{repeatedCV}}  
  \item{segments}{number of cross-validation segments used in \code{repeatedCV}}
  \item{\dots}{Arguments to be passed to methods}
}
\details{
  Some of these functions are just convenience wrappers for \code{mvrVal} 
  functions and for the \code{jack.test} function in package \pkg{pls}. 
  More details are given here: \code{\link[pls:mvrVal]{mvrVal}}, 
  \code{\link[pls:jack.test]{jack.test}}. Other functions are specific 
  \code{autopls} functions. \code{metaval} is used for a summary of validation
  results during backselection. \code{repeatedCV} is a meta cross-validation 
  (repeated ten-fold cross-validation).
}
\value{
  see \code{\link[pls:mvrVal]{mvrVal}} and \code{\link[pls:jack.test]{jack.test}}. 
  The main difference is a reduced selection of functions (see above) and the
  possibility to inherit a number of latent vectors from the \code{autopls} 
  object. 
  
  The \code{metaval} function provides a matrix overview of model results for
  all iterations and numbers of latent vectors in an \code{autopls} object.
  \code{repeatedCV} provides results and basic statistics for repeated 
  cross-validation runs.
}
\author{
  Sebastian Schmidtlein, linking to code from package \pkg{pls} by Ron Wehrens 
  and Bjørn-Helge Mevik.
}
\note{
  If you want to make full use of the \code{mvrVal} functions in the \pkg{pls} 
  package assign class \code{mvr} to the model object.
}
\seealso{
   \code{\link[pls:mvrVal]{mvrVal}}, \code{\link[pls:jack.test]{jack.test}}, 
   \code{\link{autopls}}, \code{\link[cvTools:repCV]{repCV}},
   \code{\link[chemometrics:mvr_dcv]{mvr_dcv}}
}
\examples{
  ## load predictor and response data to the current environment
  data(murnau.X)
  data(murnau.Y)
  
  ## call autopls with the standard options
  model<-autopls (murnau.Y ~ murnau.X)
  
  ## Validation
  R2 (model)
  R2 (model, nc = 'all')  
  RMSEP (model)  
  metaval (model, 'R2', 'CV', ic = FALSE)
  
  ## Jackknife test
  jack.test.autopls (model)
  
  ## Meta cross-validation
  repeatedCV (model)
}
\keyword{regression}
\keyword{multivariate}

