prepro <- function (X, prep = 'bn')
{
  
  method <- match.arg (prep, 'bn') ## more to come ...
  
  if (is.vector (X)) X <- X / sqrt (sum (X ^ 2))
  

  if (is.matrix (X)) X <- X / sqrt (rowSums (X ^ 2))         
  

  if (class (X) == 'RasterBrick' || class (X) == 'RasterStack')
    X <- X / sqrt (stackApply (X ^ 2, rep (1, nlayers (X)), sum))
  
  invisible (X)
}
