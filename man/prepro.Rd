\name{prepro}
\alias{prepro}
\title{
  Preprocessing in autopls
}
\description{
  Used for preprocessing predictor data in functions \code{autopls} and 
  \code{predict.autopls}
}
\usage{
  prepro(X, prep = NA)
}
\arguments{
  \item{X}{
  predictors as vector, matrix, raster brick or raster stack
}
  \item{prep}{
  type of preprocessing (currently only brightness normalization coded as
  \dQuote{bn})
}
}
\details{
 The function is called within \code{autopls} and \code{predict.autopls}.
 The only implemented option is currently \dQuote{bn}, which is a brightness 
 normalization according to Feilhauer et al. (2010). Raster brick and raster 
 stack are objects of package \pkg{raster}.
}
\value{
  Returns the transformed matrix or raster object.
}
\references{
  Feilhauer. H., Asner, G.P., Martin, R.E., Schmidtlein, S. (2010): 
  Brightness-normalized Partial Least Squares regression for hyperspectral 
  data. \emph{Journal of Quantitative Spectroscopy and Radiative Transfer} 
  \bold{111}: 1947--1957.
}
\author{
  Hannes Feilhauer
}
\seealso{
  \code{\link{autopls}}, \code{\link{predict.autopls}} 
}
