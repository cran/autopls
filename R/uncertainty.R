uncertainty <- function (object, prediction)
{

  ## Determine method
  if (class (prediction) == 'RasterLayer') method <- 'rst'
  if (is.vector (prediction)) method <- 'vec'
  
  ## Observed data points
  original <- object$model$Y

  dif <- abs (prediction - original [1])

  if (method == 'vec')
  {
    for (i in 2:length(original))
    {
      difnew <- abs (prediction - original [i])
      sm <- dif > difnew
      dif [sm] <- difnew [sm]      
    }
  }

  if (method == 'rst')
  {
    for (i in 2:length(original))
    {
      difnew <- abs (prediction - original [i])
      dif <- min (dif, difnew)
    }
  }
 
  return (dif)
}