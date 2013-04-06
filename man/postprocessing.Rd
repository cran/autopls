\name{postprocessing}
\alias{postprocessing}
\alias{liability}
\alias{confine}
\title{
  Test for model extrapolations or interpolations and removal of bold 
  predictions in autopls
}
\description{
  Departures of values predicted by \code{autopls} from the original data 
  space and removal of exceedingly extrapolated predictions.
}
\usage{
  liability (object, prediction)
  confine (object, prediction, tolerance)
}
\arguments{
  \item{object}{
  object of class \code{autopls}   
  }
  \item{prediction}{
  predicted values as single vector or single layer raster image 
  (\code{RasterLayer} from package \pkg{raster}). 
  }
  \item{tolerance}{
  maximum departure of preserved prediction values
  }
}
\details{
  Takes its time with large images and many objects used in calibration.
}
\value{
  Vector or raster image depending on the type of \code{prediction}
  Uncertainties (\code{liability} function) are given in original units. 
  After \code{confine}, values exceeding \code{tolerance} are replaced by 
  \code{NA}.
}
\author{
  Sebastian Schmidtlein
}
\seealso{
  \code{\link{autopls}}, \code{\link{predict.autopls}}
}
\examples{
  ## load predictor and response data to the current environment
  data (murnau.X)
  data (murnau.Y)
  
  ## call autopls with the standard options
  model <- autopls (murnau.Y ~ murnau.X)

  ## new data
  new <- murnau.X + 500 

  ## prediction
  pred <- predict (model, new)

  ## check uncertainty
  liability (model, pred)

  ## remove predictions with uncertainty value > 5
  confine (model, pred, 5)
}
\keyword{regression}
\keyword{multivariate}

