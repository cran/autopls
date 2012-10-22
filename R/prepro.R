prepro <- function (X, prep = NA)
{
  
  if (is.vector (X))
  {
    if (!is.na (prep) & prep == 'bn') X <- X / sqrt (sum (X^2))
  }
  
  ## ------------------- Matrix preprocessing ------------------------------- ##

  if (is.matrix (X))
  {
    if (!is.na (prep) & prep == 'bn') X <- X / sqrt (apply (X^2, 1, sum))
  }

  ## ------------------- Image preprocessing -------------------------------- ##

  if (class (X) == 'RasterBrick' || class (X) == 'RasterStack')
  {
    if (!is.na (prep) & prep == 'bn') X <- X / sqrt (stackApply (X^2,
      rep (1, nlayers (X)), sum))
  }
  invisible (X)
}
