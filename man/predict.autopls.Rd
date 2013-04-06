\name{predict.autopls}
\alias{predict.autopls}
\title{
  Prediction using a fitted autopls model 
}
\description{
  Applies a model from a \code{autopls} object to a vector, matrix or to a
  \code{stack} or \code{brick} from package \pkg{raster}.
}
\usage{
  \method{predict}{autopls}(object, dat, ...)
}
\arguments{
  \item{object}{
  object of class \code{autopls}   
  }
  \item{dat}{
  vector, matrix, dataframe or imagery (the latter as \code{stack} or 
  \code{brick} from package \pkg{raster}. 
  }
  \item{\dots}{logical. Arguments to be passed to method
  }
}
\details{
  Elements, columns or layers must have the same number and order as the input 
  predictors for \code{autopls}. The predictors resulting from autopls are 
  selected silently. In case of large image files the function
  is based on tile processing.  
}
\value{
  A new vector matrix or image depending on the type of \code{newdata}
}
\author{
  Sebastian Schmidtlein
}
\seealso{
  \code{\link{autopls}}
}
\examples{
  ## load predictor and response data to the current environment
  data (murnau.X)
  data (murnau.Y)
  
  ## call autopls with the standard options
  model<-autopls (murnau.Y ~ murnau.X)

  ## new data
  new <- murnau.X + 500 

  ## prediction
  predict (model, new)
}
\keyword{regression}
\keyword{multivariate}

