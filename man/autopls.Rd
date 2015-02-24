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
autopls (formula, data, testset = NULL, tselect = "none", prep = "none", 
  val = "LOO", scaling = TRUE, stingy = TRUE, verbose = TRUE, 
  backselect = "auto", jt.thresh = 0.1, vip.thresh = 0.2, jump = NA, 
  lower = NA, method = "oscorespls")
}
\arguments{
  \item{formula}{
     model formula
}
  \item{data}{
     optional data frame with the data to fit the model 
}
  \item{testset}{
     optional vector defining a test set (row indices)
}
  \item{tselect}{
     string specifying the role of the test set in model selection 
     (\code{"none"}, \code{"passive"} or \code{"active"}, see details)
}
  \item{prep}{
     character. optional preprocessing (only one choice implemented: 
     \code{"bn"} (see details)
}
  \item{val}{
     character. Validation used (\code{"CV"} or \code{"LOO"}, 
     see details)
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
     one or more character strings defining the methods used in backwards 
     selection (see details). \code{"no"} means no backselection. Defaults 
     to \code{"auto"}
}
  \item{jt.thresh}{
     threshold used in predictor selections that are based on jackknife 
     testing (methods based on A1, see details) 
}
  \item{vip.thresh}{
     threshold used in predictor selections that are based on VIP
     (methods based on A2, see details). VIP is scaled to a maximum
     of 1.
}
  \item{jump}{
     numeric. If a number is given, backward selection starts with a 
     forced reduction of predictors to the given number (see A0 in details). 
     This reduction is based on significance in jackknifing. The argument can 
     be useful in the case of large predictor matrices.
}
  \item{lower}{
     numeric. Backward selection proceeds as long as R2 in validation reaches 
     the given value (experimental, backward selection continues further if
     models improve in other respects such as decreasing numbers of latent 
     vectors).     
}
  \item{method}{
     character string indicating what plsr method to use. \code{autopls}
     works with the orthogonal scores algorithm (\code{"oscorespls"})
	 and with the kernel algorithm (\code{"kernelpls"}). 
}
}
\details{
  The \code{autopls} function is a wrapper for \link[pls:mvr]{pls} 
  in package \pkg{pls} written by Bjørn-Helge Mevik, Ron Wehrens 
  and Kristian Hovde Liland. As for now, the wrapper can be cited as 
  Schmidtlein et al. (2012). \code{autopls} works only for single target 
  variables.
  
  If validation = \dQuote{CV}, 10-fold cross-validation is performed. If 
  validation = \dQuote{LOO}, leave-one-out cross-validation is performed.
  Test set validation takes always place if a test set has been defined. 
  \code{tselect} specifies how the test set is used in model selection. 
  \code{"none"}: just use it for external validation; \code{"passive"}: 
  use error in external validation for model selection but do not use it 
  for the determination of the number of latent vectors; \code{"active"} 
  use the error in external validation for model selection and for the 
  determination of the number of latent vectors. With \code{stingy = TRUE} 
  the errors that are used in the selection are measured at a number of latent 
  vectors that depends on the number of observations (1/10 at maximum). 
  Otherwise, the number of latent vectors is chosen where errors approach 
  a first minimum. In order to avoid minor local minima the error values are 
  first smoothed.
  
  Large data matrices: Examine the arguments \code{jump} (forced reduction of
  predictors in the first iteration). Large model objects can be 
  shrinked using the function \code{\link{slim}} but some functionality (like 
  plotting or change of the number of latent vectors) is lost. Shrinked models 
  can still be used for predictions.    

  Preprocessing options: The only implemented option is currently 
  \code{"bn"}, which is a brightness normalization according to 
  Feilhauer et al. (2010).  
  
  Several methods for predictor selection are available. In default mode
  (\code{backselect = "auto"}) the selection follows an optimization procedure 
  using methods A1 and A3. However, apart from A0 any user-defined 
  combination can be selected using the \code{backselect} argument. Note that 
  VIP-based methods (A2, A3, B3 to B6) are meant to be used with the oscorespls
  method and methods B1 to B6 and C1 do only make sense with sequences of 
  spectral bands or similar sequences of autocorrelated predictors. The methods 
  are coded as follows:  
  
  \emph{A) Filtering based on thresholds}

  (A0 and A1) Based on significance, A0 with user-defined threshold (see 
  argument \code{jump}); (A2) based on VIP; (A3) based on combined 
  significance and VIP; (A4) removal of 10 \% predictors with the lowest 
  significance; (A5) removal of 25 \% predictors with the lowest significance.
  
  \emph{B) Filtering followed by reduction of autocorrelation}
  
  (B1) Filtering based on significance, thinning starting with local maxima 
  in weighted regression coefficients; (B2) filtering based on significance, 
  thinning starting with local maxima in significance; (B3) filtering based 
  on significance, thinning starting with local maxima in VIP; (B4) filtering 
  based on VIP, thinning starting with local maxima in weighted regression 
  coefficients; (B5) filtering based on VIP, thinning starting with local 
  maxima in significance; (B6) filtering based on VIP, thinning starting with 
  local maxima in VIP.

  \emph{C) Just reduction of autocorrelation}
  
  (C1): reduction starting with local maxima in regression 
  coefficients.
    
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
  \item{X.testset}{test set: predictors}
  \item{Y.testset}{test set: target variable}
  \item{preprocessing}{method used for preprocessing}
  \item{scaling}{\code{TRUE} if scaling was requested}
  \item{val}{\code{LOO} or \code{CV}}
  \item{call}{the function call}
}
\references{
  Feilhauer. H., Asner, G.P., Martin, R.E., Schmidtlein, S. (2010): 
  Brightness-normalized Partial Least Squares regression for hyperspectral 
  data. \emph{Journal of Quantitative Spectroscopy and Radiative Transfer} 
  \bold{111}: 1947--1957.

  Schmidtlein, S., Feilhauer, H., Bruelheide, H. (2012): Mapping plant
  strategy types using remote sensing. \emph{Journal of Vegetation Science}
  \bold{23}: 395--405. Open Access.
}
\author{
  Sebastian Schmidtlein with contributions from Carsten Oldenburg 
  and Hannes Feilhauer. The code for computing \code{VIP} is borrowed from 
  Bjørn-Helge Mevik.
}

\seealso{
   \code{\link[pls:mvr]{pls}}, \code{\link{set.iter}}, \code{\link{set.lv}}, 
   \code{\link{predict.autopls}}, \code{\link{plot.autopls}}
}
\examples{
  ## load predictor and response data to the current environment
  data (murnau.X)
  data (murnau.Y)
  
  ## call autopls with the standard options
  model <- autopls (murnau.Y ~ murnau.X)
  
  ## S3 plot method
  \dontrun{plot (model)}
  \dontrun{plot (model, type = "rc")}
  
  ## Loading and score plots
  \dontrun{plot (model$loadings, main = "Loadings")}
  \dontrun{plot (model$loadings [,c(1,3)], main = "Loadings")}
  \dontrun{plot (model$scores, main = "Scores")}
  
}
\keyword{regression}
\keyword{multivariate}

