\name{plot.autopls}
\alias{plot.autopls}
\title{
  Plotting function for autopls objects
}
\description{
  Produces plots illustrating the outcomes of \code{\link{autopls}}: 
  predicted vs. observed values, errors vs. numbers of latent vectors, 
  regression coefficients, influences of observations regarding X and Y,
  latent vectors and R2 in backward selection 
}
\usage{
  \method{plot}{autopls}(x,type="all",wl=NULL,rcxlab = "Predictors",
  plab=FALSE, bw = FALSE, ...)
}
\arguments{
  \item{x}{
  object of class \code{autopls}     
}
  \item{type}{specifying the type of plot. (\code{"all"}: all plot; 
     \code{"ovp"}: observed vs. predicted values; 
     \code{"rmse"}: errors vs. numbers of latent vectors; 
     \code{"rc"}: regression coefficients; 
     \code{"x.inf"}: influence plot (X-variance); 
     \code{"y.inf"}: influence plot (Y-variance); 
     \code{"meta"}: latent vectors and R2 in backward selection)}   
  \item{wl}{denotes an optional vector of numerical values 
     describing the position of predictors along the x axis in the \code{rc} 
     plot. The values should refer to all bands (before backward selection) or
     to the bands that are actually used.} 
  \item{rcxlab}{Label for x axis in rc plot.} 
  \item{plab}{logical. Whether observations are labeled.} 
  \item{bw}{logical. Whether plots are given in grey-scales (partly
     realized).} 
  \item{\dots}{Arguments to be passed to methods}
}
\details{
  Red dots in the influence plots indicate potentially dangerous outliers
}
\value{
  Apart from the plots the function returns the underlying values
}
\author{
  Sebastian Schmidtlein, Carsten Oldenburg
  
  The placement of observation labels if \code{plab = TRUE} is done using 
  code borrowed from the \code{pointLabel} in \pkg{maptools}. The author 
  of this function is Tom Short (EPRI).
}
\seealso{
  \code{\link{autopls}}
}
\examples{
  ## load predictor and response data to the current environment
  data(murnau.X)
  data(murnau.Y)
  data(murnau.W)
  
  ## call autopls with the standard options
  model<-autopls(murnau.Y,murnau.X)
  
  ## plot results
  \dontrun{plot (model)}

  ## use wavelengths in rc plot
  \dontrun{plot (model, type = "rc", wl = murnau.W, rcxlab = "Wavelength (nm)")}
  
  ## predicted vs. observed
  \dontrun{x <- plot (model, type = "ovp")}
  \dontrun{x}
  
}
\keyword{regression}
\keyword{multivariate}
\keyword{hplot}
