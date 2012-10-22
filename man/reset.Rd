\name{reset}
\alias{reset}
\title{
  Resets the number of latent vectors and the iteration used in autopls objects
}
\description{
  Resets the number of latent vectors and the iteration used in autopls objects
  to the values originally selected by the autopls procedure 
}
\usage{
  reset(object, verbose = TRUE)
}
\arguments{
  \item{object}{object of class \code{autopls}   
}
  \item{verbose}{logical. If a summary of the resulting object should be
        printed on the screen
}
}
\value{
  Returns an object of class \code{autopls} 
}
\author{
  Sebastian Schmidtlein
}
\seealso{
  \code{\link{autopls}}, \code{\link{set.iter}}, \code{\link{set.lv}}
}
\examples{
  ## load predictor and response data to the current environment
  data(murnau.X)
  data(murnau.Y)
  
  ## call autopls with the standard options
  model<-autopls(murnau.Y,murnau.X)
  
  ## select another iteration
  newmodel <- set.iter(model,3)
  
  ## set another number of latent vectors
  evennewermodel <- set.lv(newmodel,2)

  ## return to the original values
  firstmodel <- reset(evennewermodel)

}
\keyword{regression}
\keyword{multivariate}
