## Summary and print functions for autopls objects

summary.autopls <- function (object, ...)
{

  ## Get parameters
  iter <- object$metapls$current.iter
  lv <- get.lv (object)
  N <- length (object$metapls$Y)
  pred <- sum (object$predictors)
  r2.all <- unlist (R2 (object, c('train', 'CV'), nc = lv, ic = FALSE))
  r2.cal <- r2.all$val1
  r2.val <- r2.all$val2
  rmse.all <- unlist (RMSEP (object, c('train', 'CV'), nc = lv, ic = FALSE))
  rmse.cal <- rmse.all$val1
  rmse.val <- rmse.all$val2
  scaling <- object$metapls$scaling
  preprocessing <- object$metapls$preprocessing

  ## Output object
  output <- list (
  predictors = pred,
  observations = N,
  lv = lv,
  run = iter,
  rmse.cal = rmse.cal,
  rmse.val = rmse.val,
  r2.cal = r2.cal,
  r2.val = r2.val,
  scaling = scaling,
  preprocessing = preprocessing)

  return (output)
}

print.autopls <- function (x, ...)
{

  ## Get parameters
  iter <- x$metapls$current.iter
  lv <- get.lv (x)
  N <- length (x$metapls$Y)
  pred <- sum (x$predictors)
  r2.all <- unlist (R2 (x, c('train', 'CV'), nc = lv, ic = FALSE))
  r2.cal <- r2.all$val1
  r2.val <- r2.all$val2
  rmse.all <- unlist (RMSEP (x, c('train', 'CV'), nc = lv, ic = FALSE))
  rmse.cal <- rmse.all$val1
  rmse.val <- rmse.all$val2

  ## Screen output
  cat ('\n')
  cat (paste ('Predictors:', pred, '  '))
  cat (paste ('Observations:', N, '  '))
  cat (paste ('Latent vectors:', lv, '  '))
  cat (paste ('Run:', iter, '\n'))
  cat (paste ('RMSE in calibration:', format (rmse.cal, digits = 6), '   '))
  cat (paste ('RMSE in validation:', format (rmse.val, digits = 6), '\n'))
  cat (paste ('R2 in calibration:', format (r2.cal, digits = 6), '   '))
  cat (paste ('R2 in validation:', format (r2.val, digits = 6), '\n\n'))

  ## Output
  output <- list (
  predictors = pred,
  lv = lv,
  rmse.cal = rmse.cal,
  rmse.val = rmse.val,
  r2.cal = r2.cal,
  r2.val = r2.val)
  invisible (output)
}

summary.slim <- function (object, ...)
{

  ## Get parameters
  iter <- object$metapls$current.iter
  lv <- get.lv (object)
  N <- object$slimobj$N
  pred <- sum (object$predictors)
  r2.cal <- object$slimobj$r2$val [1]
  r2.val <- object$slimobj$r2$val [2]
  rmse.cal <- object$slimobj$rmse$val [1]
  rmse.val <- object$slimobj$rmse$val [2]

  ## Output
  output <- list (
  predictors = pred,
  lv = lv,
  rmse.cal = rmse.cal,
  rmse.val = rmse.val,
  r2.cal = r2.cal,
  r2.val = r2.val)
  
  return (output)
}

print.slim <- function (x, ...)
{

  ## Get parameters
  iter <- x$metapls$current.iter
  lv <- get.lv (x)
  N <- x$slimobj$N
  pred <- sum (x$predictors)
  r2.cal <- x$slimobj$r2$val [1]
  r2.val <- x$slimobj$r2$val [2]
  rmse.cal <- x$slimobj$rmse$val [1]
  rmse.val <- x$slimobj$rmse$val [2]

  ## Screen output
  cat ('\n')
  cat (paste ('Predictors:', pred, '  '))
  cat (paste ('Observations:', N, '  '))
  cat (paste ('Latent vectors:', lv, '  '))
  cat (paste ('Run:', iter, '\n'))
  cat (paste ('RMSE in calibration:', format (rmse.cal, digits = 6), '   '))
  cat (paste ('RMSE in validation:', format (rmse.val, digits = 6), '\n'))
  cat (paste ('R2 in calibration:', format (r2.cal, digits = 6), '   '))
  cat (paste ('R2 in validation:', format (r2.val, digits = 6), '\n\n'))

  ## Output
  output <- list (
  predictors = pred,
  lv = lv,
  rmse.cal = rmse.cal,
  rmse.val = rmse.val,
  r2.cal = r2.cal,
  r2.val = r2.val)
  invisible (output)
}
