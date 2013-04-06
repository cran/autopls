\name{predict.slim}
\alias{predict.slim}
\title{
  Prediction using a condensed autopls model 
}
\description{
  Applies a model object of class \code{slim} originating from 
  \code{autopls} to a vector, matrix or to a \code{stack} or \code{brick} 
  from package \pkg{raster}.
}
\usage{
  \method{predict}{slim}(object, dat, ...)
}
\arguments{
  \item{object}{
  object of class \code{slim}   
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
  predictors for \code{autopls}. In case of large image files the function
  is based on tile processing.  
}
\value{
  A new vector matrix or image depending on the type of \code{newdata}
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
  model <- autopls (murnau.Y ~ murnau.X)

  ## condensed model object
  new.model <- slim (model)
  
  ## new data
  new <- murnau.X + 500 

  ## prediction
  predict (new.model, new)
}
\keyword{regression}
\keyword{multivariate}

