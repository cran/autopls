\name{summary.autopls}
\alias{summary.autopls}
\alias{summary.slim}
\alias{print.autopls}
\alias{print.slim}
\title{
  Summary and print functions for autopls objects
}
\description{
  Summary and print methods for \code{autopls} and \code{slim} objects.
}
\usage{
  \method{summary}{autopls}(object, \dots)
  \method{summary}{slim}(object, \dots)
  \method{print}{autopls}(x, \dots)
  \method{print}{slim}(x, \dots)
}
\arguments{
  \item{object}{object of class \code{autopls} or \code{slim} respectively}   
  \item{x}{object of class \code{autopls} or \code{slim}}   
  \item{\dots}{Arguments to be passed to methods}   
}
\value{
  \item{predictors}{number of predictors used in the final model}
  \item{lv}{number of latent vectors used in the final model}
  \item{rmse.cal}{root mean squared errors in calibration}
  \item{rmse.val}{root mean squared errors in validation}
  \item{r2.cal}{R2 in calibration}
  \item{r2.val}{R2 in validation}  
  
  Print returns a screen output and an unvisible object with the 
  same content
}
\author{
  Sebastian Schmidtlein
}
\seealso{
  \code{\link{autopls}}, \code{\link{slim}}
}
\examples{
  ## load predictor and response data to the current environment
  data (murnau.X)
  data (murnau.Y)
  
  ## call autopls with the standard options
  model <- autopls(murnau.Y ~ murnau.X)
  
  ## print and plot results
  print (model)
}
\keyword{regression}
\keyword{multivariate}

