\encoding{latin1}
\name{autopls}
\alias{autopls}
\title{
autopls
}
\description{
Partial least squares regression with backward selection of predictors
}
\usage{
autopls(resp, pred, prep = NA, val = "LOO", spectral = TRUE, 
scaling = TRUE, stingy = TRUE, verbose = TRUE, backselect = TRUE,
jump = NA, thorough = FALSE)
}
\arguments{
  \item{resp}{
     Vector with single response variable
}
  \item{pred}{
     Matrix with predictors (number and order of rows fitting \code{resp}) 
}
  \item{prep}{
     character. optional preprocessing (only one choice implemented: 
     \code{"bn"} (see details)
}
  \item{val}{
     character. Validation used (\code{dQuote{CV}} or \code{dQuote{LOO}}, 
     see details)
}
  \item{spectral}{
     logical. If \code{TRUE}, predictors are assumed to be reflectance
     (implications for methods used in backward selection)
}
  \item{scaling}{
     logical. if \code{TRUE}, predictors are scaled by dividing each variable 
     by its standard deviation. This is repeated in all validation steps
}
  \item{stingy}{
     logical. If \code{TRUE}, the number of latent vectors is kept low during 
     backward selection
}
  \item{verbose}{
     logical. If \code{TRUE}, details about the backward selection processes 
     are reported
}
  \item{backselect}{
     logical. If \code{TRUE}, backward selection of predictors takes place
}
  \item{jump}{
     numeric. If a number is given, backward selection starts with a 
     forced reduction of predictors to the given number. This reduction is 
     based on significance in jackknifing. The argument can be useful in the 
     case of large predictor matrices.
}
  \item{thorough}{
     logical. if \code{TRUE}, more methods of variable selection are used.
     This slows down the process and may sometimes increase the chance of 
     finding good solutions. May also end up with decreased model quality.
}
}
\details{
  The \code{autopls} function is a wrapper for \link[pls:mvr]{pls} 
  in package \pkg{pls} written by Bjørn-Helge Mevik, Ron Wehrens 
  and Kristian Hovde Liland. As for now, the wrapper can be cited as 
  Schmidtlein et al. (2012). \code{autopls} works only for single target 
  variables.
  
  Preprocessing options: The only implemented option is currently 
  \code{"bn"}, which is a brightness normalization according to 
  Feilhauer et al. (2010).  
  
  If validation = \dQuote{CV}, 10-fold cross-validation is performed. If 
  validation = \dQuote{LOO}, leave-one-out cross-validation is performed.

  Large data matrices: Examine the arguments \code{jump} (forced reduction of
  predictors in the first iteration) and \code{thorough} (reduced number of
  methods considered for variable selection). With very large data sets 
  \code{autopls} may run out of memory, especially if \code{val = "LOO"}. 
  \code{val = "CV"} should work in most cases. Large model objects can be 
  shrinked using the function \code{\link{slim}}.    

  During backward selection, several methods for predictor selection 
  are tested. They are coded as follows:  
  
  \emph{A) Filtering based on thresholds}

  (A0 and A1) Based on significance, A0 with user-defined threshold (see 
  argument \code{jump}); (A2) based on VIP; (A3) based on combined 
  significance and VIP; (A4) removal of 10 \% predictors with the lowest 
  significance; (A5) removal of 25 \% predictors with the lowest significance.
  
  \emph{B) Filtering followed by reduction of autocorrelation}
  
  These are only used with \code{spectral = TRUE}: (B1) filtering based on 
  significance, thinning starting with local maxima in weighted regression 
  coefficients; (B2) filtering based on significance, thinning starting with 
  local maxima in significance; (B3) filtering based on significance, thinning 
  starting with local maxima in vip; (B4) filtering based on VIP, thinning 
  starting with local maxima in weighted regression coefficients; 
  (B5) filtering based on VIP, thinning starting with local maxima in 
  significance; (B6) filtering based on VIP, thinning starting with local 
  maxima in VIP.

  \emph{C) Just reduction of autocorrelation}
  
  (C1): reduction starting with local maxima in regression 
  coefficients.
  
  With \code{thorough = FALSE} A1, A2, A3, A4 and B2 are used. With 
  \code{spectral = FALSE} only A-methods are used.
    
}  
\value{
  An object of class \code{autopls} is returned. This equals a 
  \link[pls:mvr]{pls} object and some added objects:

  \item{predictors}{logical. Vector of predictors that have been or have not
               been used in the current model}
  \item{metapls}{outcomes of the backward selection process}
  \item{iterations}{models selected during the backward selection process}
  
  The \code{$metapls} item consists of the following:

  \item{current.iter}{iteration of the backward selection procedure the
               current model is based upon}
  \item{autopls.iter}{iteration of the backward selection procedure originally
               selected by autopls}
  \item{current.lv}{number of latent vectors the current model is based upon}
  \item{autopls.lv}{number of latent vectors originally selected by autopls}
  \item{lv.history}{sequence of number of latent vectors values selected during 
               iterations in backward selection}
  \item{rmse.history}{sequence of root mean squared errors obtained during 
               iterations in backward selection. Errors are reported for
               calibration and validation. The validation errors are also 
               reported for the number of latent vectors corresponding to 
               \code{ceiling (nrow (pred) / 10)}.}
  \item{r2.history}{sequence of number of r2 values obtained during
               iterations in backward selection}
  \item{X}{original predictors}
  \item{Y}{original target variable}
  \item{preprocessing}{method used for preprocessing}
  \item{scaling}{\code{TRUE} if scaling was requested}
  \item{call}{the function call}
}
\references{
  Feilhauer. H., Asner, G.P., Martin, R.E., Schmidtlein, S. (2010): 
  Brightness-normalized Partial Least Squares regression for hyperspectral 
  data. \emph{Journal of Quantitative Spectroscopy and Radiative Transfer} 
  \bold{111}: 1947--1957.

  Schmidtlein, S., Feilhauer, H., Bruelheide, H. (2012): Mapping plant
  strategy types using remote sensing. \emph{Journal of Vegetation Science}
  \bold{23}: 395--604. Open Access.
}
\author{
  Sebastian Schmidtlein with contributions from Carsten Oldenburg and 
  Hannes Feilhauer. A code snipped for computing \code{VIP} is borrowed from 
  Bjørn-Helge Mevik.
}

\seealso{
   \code{\link[pls:mvr]{pls}}, \code{\link{set.iter}}, \code{\link{set.lv}}, 
   \code{\link{predict.autopls}}, \code{\link{plot.autopls}}
}
\examples{
  ## load predictor and response data to the current environment
  data(murnau.X)
  data(murnau.Y)
  
  ## call autopls with the standard options
  model<-autopls(murnau.Y, murnau.X)
  
  ## Plot
  \dontrun{plot (model)}
}
\keyword{regression}
\keyword{multivariate}

