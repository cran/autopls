\name{set.iter}
\alias{set.iter}
\title{
  Sets the run of an autopls backwards selection to be used
}
\description{
  Changes the run of an autopls backwards selection to be used
}
\usage{
  set.iter (object, iteration, verbose = TRUE)
}
\arguments{
  \item{object}{object of class \code{autopls}   
}
  \item{iteration}{new value for the iteration used
}
  \item{verbose}{logical. If a summary of the resulting object should be
        printed on the screen
}
}
\details{
  The number of latent vectors is set to the orinal number for this run. 
}
\value{
  Returns an object of class \code{autopls} 
}
\author{
  Sebastian Schmidtlein
}
\seealso{
  \code{\link{autopls}}, \code{\link{set.lv}}
}
\examples{
  ## load predictor and response data to the current environment
  data(murnau.X)
  data(murnau.Y)
  
  ## call autopls with the standard options
  model<-autopls(murnau.Y,murnau.X)
  
  ## set another number of latent vectors
  newmodel <- set.iter(model,3)
}
\keyword{regression}
\keyword{multivariate}
