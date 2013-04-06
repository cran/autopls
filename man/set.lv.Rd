\name{set.lv}
\alias{set.lv}
\title{
  Sets the number of latent vectors in autopls objects
}
\description{
  Sets the number of latent vectors in a \code{autopls} object
}
\usage{
  set.lv(object, lv, verbose = TRUE)
}
\arguments{
  \item{object}{object of class \code{autopls}   
}
  \item{lv}{new value for the number of latent vectors used
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
  \code{\link{autopls}}, \code{\link{set.iter}}
}
\examples{
  ## load predictor and response data to the current environment
  data (murnau.X)
  data (murnau.Y)
  
  ## call autopls with the standard options
  model<-autopls (murnau.Y ~ murnau.X)
  
  ## set another number of latent vectors
  newmodel <- set.lv (model,2)
}
\keyword{regression}
\keyword{multivariate}
