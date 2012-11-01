\encoding{latin1}
\name{extract.autopls}
\alias{extract.autopls}
\alias{coef.autopls}
\alias{coef.slim}
\alias{scores.autopls}
\alias{loadings.autopls}
\alias{fitted.autopls}
\alias{get.lv}
\alias{get.iter}
\alias{residuals.autopls}
\alias{predicted}
\alias{slim}
\title{
  Extract information from a fitted autopls model
}
\description{
  Functions to extract information from \code{autopls} objects: 
  crossvalidation, fitted values, regression coefficients,
  residuals, scores, loadings, latent vectors used, underlying run.
}
\usage{
  predicted (object)
  get.lv (object)
  get.iter (object)
  slim (object)
  \method{scores}{autopls}(object, ...)
  \method{loadings}{autopls}(object, ...)
  \method{fitted}{autopls}(object, ...)
  \method{coef}{autopls}(object, intercept = FALSE, ...)
  \method{coef}{slim}(object, intercept = FALSE, ...)
  \method{residuals}{autopls}(object, ...)
}
\arguments{
  \item{object}{object of class \code{autopls}
  }
  \item{intercept}{logical. Should intercept be given?
  }
  \item{\dots}{logical. Arguments to be passed to methods
  }
}
\details{
  Provides convenience wrappers for \code{extract} functions in package 
  \pkg{pls}. More details are given here: \link[pls:coef.mvr]{coef.mvr}.
  Other functions extract information specific for \code{autopls}
  objects: \code{get.lv}, \code{get.iter} or condense the model information
  to a memory saving object of class \code{slim} that can be used for
  predictions with \code{predict.slim}. This makes sense if large
  pedictor data sets result in huge \code{autopls} model objects that
  are difficult to handle.
}
\value{
  see \link[pls:coef.mvr]{coef.mvr}. \code{get.iter} returns the run in the 
  \code{autopls} backwards selection procedure that has been used for the 
        current model. 
  \code{get.lv} returns the number of latent vectors used for the present model. 
  \code{predicted} returns the predictions in model validation while 
  \code{fitted} returns the predictions in model calibration.
  \code{slim} returns an object of class \code{slim}.
}
\author{
  Sebastian Schmidtlein, links to code from package \pkg{pls} by 
  Ron Wehrens and Bjørn-Helge Mevik.
}
\note{
  If you want to make full use of the \code{extract} functions in the \pkg{pls} 
  package assign class \code{mvr} to the model object.
  
  Reducing a model to an object of class \code{slim} means loosing plotting
  options.  
}
\seealso{
  \code{\link{autopls}}, \code{\link{metaval}}, \code{\link{set.iter}}, 
  \code{\link{set.lv}}, \code{\link{predict.slim}}
}

\examples{
  ## load predictor and response data to the current environment
  data (murnau.X)
  data (murnau.Y)
  
  ## call autopls with the standard options
  model <- autopls (murnau.Y ~ murnau.X)
  
  ## get fitted values
  fitted(model)
}
