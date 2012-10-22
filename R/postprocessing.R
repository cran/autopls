uncertainty <- function (object, prediction)
{

  ## Determine method
  if (is.vector (prediction)) method <- 'vec'
  if (class (prediction) == 'RasterLayer') method <- 'rst'
  
  ## Observed data points
  original <- object$model$Y

  dif <- abs (prediction - original [1]) ## Initialization

  if (method == 'vec') ## Predictions resulting from a matrix of predictors
  {
    for (i in 2:length(original))
    {
      difnew <- abs (prediction - original [i])
      sm <- dif > difnew
      dif [sm] <- difnew [sm]  ## Returns a vector    
    }
  }

  if (method == 'rst') ## Predictions resulting from a multi-layer image
  {
    for (i in 2:length(original))
    {
      difnew <- abs (prediction - original [i])
      dif <- min (dif, difnew) ## Returns a RasterLayer
    }
  }
 
  return (dif)
}

cleanup <- function (object, prediction, tolerance)
{
  
  dif <- uncertainty (object = object, prediction = prediction)

  idx <- dif <= tolerance
  
  if (is.vector (dif)) prediction [!idx] <- NA
  if (class (dif) == 'RasterLayer') prediction <- idx * prediction
    
  return (prediction)
} 