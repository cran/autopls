autopls <- function (resp, pred, prep = NA, val = 'LOO',
  spectral = TRUE, scaling = TRUE, stingy = TRUE, verbose = TRUE,
  backselect = TRUE, jump = NA, thorough = FALSE)
{

  ## Some data preparation
  X <- as.matrix (pred)
  if (is.null (colnames (X))) colnames (X) <- 1:ncol (X)
  Y <- resp
  if (!is.null (dim (Y)[2])) 
  {
    if (dim (Y)[2] == 1) Y <- as.vector (t (Y))
    else stop ('Found multiple target variables. PLS2 not implemented')
  }  
  
  ## Check data suitability
  if (dim (X)[1] != length (Y)) stop ('Samples do not match')
  if (min (apply (X, 2, var)) == 0) stop ('Predictor(s) with no variance')
  if (min (apply (X, 1, var)) == 0) warning ('Observation(s) with no variance')

  ## Determine acceptable number of latent vectors
  maxnlv <- ceiling (nrow (X) * 0.1)
  ## With large matrices, results are better if autopls searches solutions 
  ## with less than n*0.1 latent vectors:
  if (maxnlv > 5) maxnlv <- log (nrow (X), base = 2)

  ## --- FUNCTIONS ---------------------------------------------------------- ##

  ## Extract RMSE
  getrmse <- function (model)
  {
    rmse <- RMSEP (model, c('train', 'CV'), intercept = FALSE)$val [1:2,,]
    if (ncol (rmse) > 40) rmse <- rmse [,1:40]
    return (rmse)
  }

  ## Extract R2
  getr2 <- function (model)
  {
    r2 <- R2 (model, c('train', 'CV'), intercept = FALSE)$val [1:2,,]
    if (ncol (r2) > 40) r2 <- r2 [,1:40]
    return (r2)
  }

  ## Determine number of latent vectors where RMSEval minimizes (nlv)
  getbest <- function (rmse)
  {
    v <- rmse [2,]
    ## Check (1) if RMSE decreases monotonically, if (2) rmse are too short for
    ## filtering
    slope <- v [2:length (v)] - v [1:(length (v) - 1)]
    if (max (slope) > 0 | length (v) < 4)
    {
      ## Variable determining 'conservatisvism' of proposed nlv
      cnsv <- 0.005
      criterion <- cnsv * (max (v) - min (v))
      dif <- v [2:length (v)] + criterion - v [1:length (v) - 1]
      nlv <- match (FALSE, dif < 0)
      if (is.na (nlv)) nlv <- 1  ## should not happen
    }
    else
    {
      ## Filter for smoothing
      flt <- c(1, 2, 1)
      ## Low-pass filter (avoiding conflicts with package signal)
      lp <- stats::filter (v, flt)
      
      ## First minimum
      dlp <- diff (lp)
      minflt <- match (FALSE, dlp < 0)
      if (is.na (minflt)) nlv <- 1
      ## First minimum in filter region
      else
      {
        dwdw <- diff (v [(minflt - 1):(minflt + 1)])
        nlv <- c(minflt-1, minflt, minflt+1) [match (FALSE, dwdw < 0)]
        if (is.na (nlv)) nlv <- 1
      }
    }
    return (nlv)
  }

  ## Validation
  checkout <- function (xdat, selcode)
  {
    chkX <- xdat [,selcode == 1]
    chk <- plsr (Y ~ chkX, jackknife = TRUE,
      scale = scaling, validation = 'LOO', method = 'oscorespls')
    chk.rmse <- getrmse (chk)
    chk.best <- getbest (chk.rmse)
    ## returns (1.) error at maxnlv, (2.) nlv where error minimizes and
    ## (3.) error at that point
    if (ncol (chk.rmse) >= maxnlv)
      chk.val <- c(chk.rmse [2, maxnlv], chk.best, chk.rmse [2, chk.best])
    ## If the number of available LV's is smaller than the critical value
    ## use the lowest absolute error instead of error at maxnlv
    else 
      chk.val <- c(chk.rmse [2, chk.best], chk.best, chk.rmse [2, chk.best])        
    return (chk.val)
  }

  ## Reduce autocorrelation
  ac <- function (vec, selcode, pow = 1)
  {
    ## Function for finding local maxima
    locmax <- function (vec)
    {
      kern <- 5
      stp <- (kern - 1) / 2
      rca <- as.vector (vec)
      lthrca <- length (rca)

      ## Treat positive and negative values separately
      rcapos <- rca
      rcaneg <- rca
      rcapos [rca <= 0] <- NA
      rcaneg [rca >= 0] <- NA
      rcaneg <- abs (rcaneg)

      ## Function
      stepthrough <- function (vec, stp, lthrca)
      {
        loc2 <- vector ()
        for (i in (stp+1):(lthrca-stp))
        {
          getit <- vec [i] == max (vec [(i-stp):(i+stp)], na.rm = TRUE)
          loc2 <- c (loc2, getit)
        }
        ## Deal with the ends
        loc1 <- vector ()
        for (j in 1:stp)
        {
          getit <- vec [j] == max (vec [1:(j+stp)], na.rm = TRUE)
          loc1 <- c (loc1, getit)
        }
        loc3 <- vector ()
        for (k in (lthrca-stp+1):lthrca)
        {
          getit <- vec [k] == max (vec [(k-stp):lthrca], na.rm = TRUE)
          loc3 <- c (loc3, getit)
        }
        loc <- c(loc1, loc2, loc3)
        loc <- which (loc == TRUE)
        return (loc)
      }
      locpos <- suppressWarnings (stepthrough (rcapos, stp, lthrca))
      locneg <- suppressWarnings (stepthrough (rcaneg, stp, lthrca))

      ## Local maxima:
      loc <- sort (c(locpos, locneg))
      return (loc)
    }

    ## Compute correlation matrix of predictors
    cmat <- cor (subX)

    ## Mean nearest neighbor
    diag (cmat) <- -1
    nn <- apply (cmat, 1, max)
    mnn <- median (nn)  ## Change this if another threshold is needed

    vec.thinout <- vec
    vec.thinout [selcode == 0] <- NA
    thin <- vector ()
    prog <- TRUE
    while (prog == TRUE)
    {
      ## Get local maxima from regression coefficients
      loc <- locmax (vec.thinout)
      ## Add this selection to overall selection
      thin <- c (thin, loc)
      ## Correlation to local maxima larger than criterion?
      snn <- cmat [loc,] >= (mnn / pow)
      ## Remove such predictors
      if (length (loc) > 1) vec.thinout [apply (snn, 2, sum) > 0] <- NA
      else vec.thinout [snn == TRUE] <- NA
      ## Remove previous local maxima
      vec.thinout [thin] <- NA
      ## Check if there are at least two predictors left
      if (length (which (!is.na (vec.thinout))) < 3) prog <- FALSE
    }

    newvec <- rep (0, length (vec))
    newvec [thin] <- 1
    return (newvec)
  }

  ## Dynamic significance filter
  dynp <- function (vlth, jkn, stp)
  {
    th <- max (jkn)
    stopcrit <- vlth
    while (stopcrit > stp)
    {
      vec <- rep (1, vlth)
      vec [jkn >= th] <- 0
      th <- th - 0.0001
      stopcrit <- sum (vec)
    }
    return (vec)
  }
      
tryaround <- function (subX, reg.coef, vip, jt)
  {
    seln <- ncol (subX)

    ## Thresholds (change if appropriate)
    jt.thresh <- 0.1
    vip.thresh <- 0.2
    comb.thresh <- - 0.1

    ## In jumpmode, try to reduce drastically in the first iteration
    jumpmode <- FALSE
    if (!is.na (jump) & counter == 1 & seln > jump) jumpmode <- TRUE
    
    if (jumpmode)
    {
      selmat <- matrix (1, nrow = 1, ncol = seln)
      selmat <- dynp (vlth = seln, jkn = jt, stp = jump)        
      result <- list (selmat, 'A0')
    }

    else ## if not in jumpmode
    {
      if (!thorough)  ## if not thorough (default)
      {
        if (spectral)
        {        
          compn <- c('A1','A2','A3','A4','B2')
          selmat <- matrix (1, nrow = 5, ncol = seln)
          rownames (selmat) <- compn
          selmat [1, jt >= jt.thresh] <- 0 ## A1
          selmat [2, vip < vip.thresh] <- 0  ## A2
          selmat [3, vip - jt < comb.thresh] <- 0  ## A3
          selmat [4,] <- dynp (vlth = seln, jkn = jt, stp = seln * 0.9) ## A4            
          selmat [5,] <- ac (jt, selmat [1,]) ## B2
        }
        else ## if not thorough and not spectral
        {        
          compn <- c('A1','A2','A3','A4')
          selmat <- matrix (1, nrow = 4, ncol = seln)
          rownames (selmat) <- compn
          selmat [1, jt >= jt.thresh] <- 0 ## A1
          selmat [2, vip < vip.thresh] <- 0  ## A2
          selmat [3, vip - jt < comb.thresh] <- 0  ## A3
          selmat [4,] <- dynp (vlth = seln, jkn = jt, stp = seln * 0.9) ## A4            
        }
        
      }
      else ## if thorough (comprehensive mode)
      { 
        if (spectral)
        {
          compn <- c('A1','A2','A3','A4','A5','B1','B2','B3','B4','B5','B6','C1')
          selmat <- matrix (1, nrow = 12, ncol = seln)
          rownames (selmat) <- compn
          selmat [1, jt >= jt.thresh] <- 0  ## A1                          
          selmat [2, vip < vip.thresh] <- 0  ## A2
          selmat [3, vip - jt < comb.thresh] <- 0  ## A3
          selmat [4,] <- dynp (vlth = seln, jkn = jt, stp = seln * 0.9) ## A4            
          selmat [5,] <- dynp (vlth = seln, jkn = jt, stp = seln * 0.75) ## A5
          selmat [6,] <- ac (reg.coef, selmat [1,]) ## B1
          selmat [7,] <- ac (jt, selmat [1,]) ## B2
          selmat [8,] <- ac (vip, selmat [1,]) ## B3
          selmat [9,] <- ac (reg.coef, selmat [2,]) ## B4
          selmat [10,] <- ac (jt, selmat [2,]) ## B5
          selmat [11,] <- ac (vip, selmat [2,]) ## B6
          selmat [12,] <- ac (reg.coef, rep (1, seln), pow = 2) ## C1
        }                        
        else
        {
          compn <- c('A1','A2','A3','A4','A5')
          selmat <- matrix (1, nrow = 5, ncol = seln)        
          rownames (selmat) <- compn        
          selmat [1, jt >= jt.thresh] <- 0  ## A1
          selmat [2, vip < vip.thresh] <- 0  ## A2
          selmat [3, vip - jt < comb.thresh] <- 0  ## A3
          selmat [4,] <- dynp (vlth = seln, jkn = jt, stp = seln * 0.9) ## A4            
          selmat [5,] <- dynp (vlth = seln, jkn = jt, stp = seln * 0.75) ## A5
        }        
      }    
    
      ## Remove solutions without or only one selection or no reduction
      sums <- apply (selmat, 1, sum)
      use <- (sums != ncol (selmat)) == (sums > 1)            
      selmat <- selmat [use,]      
      compn <- compn [use]

      ## Unique solutions
      if (is.null (dim (selmat)))
      { 
        unq <- 1
        checkmat <- matrix (NA, nrow = 1, ncol = 3)
        checkmat <- checkout (subX, selmat)       
        ## errors at maxnlv
        errors <- checkmat [1]    
        ## nlv where errors in validation minimize
        comps <- checkmat [2]    
        ## errors at minimum
        errors.min <- checkmat [3]
        result <- list (selmat, names (use) [use])        
      }
      else 
      {
        code <- apply (selmat, 1, FUN = function(x) paste (x, collapse = ''))
        grps <- as.numeric (as.factor (code))
        selmat <- selmat [order (grps),]
        grps <- sort (grps)
        unq <- match (1:max(grps), grps)
        
        ## Matrix for the results of checkout ()
        checkmat <- matrix (NA, nrow = nrow (selmat), ncol = 3)
        rownames (checkmat) <- rownames (selmat)
        
        ## Use function checkout only if selcode is met for the first time
        for (i in 1:nrow(checkmat))
        {
          if (i %in% unq) checkmat [i,] <- checkout (subX, selmat [i,])
          else checkmat [i,] <- checkmat [i-1,]
        }
        
        ## Rearrange checkmat and selmat
        checkmat <- checkmat [order (rownames (checkmat)),]
        selmat <- selmat [order (rownames (selmat)),]
        
        ## errors at maxnlv
        errors <- checkmat [,1]    
        ## nlv where errors in validation minimize
        comps <- checkmat [,2]    
        ## errors at minimum
        errors.min <- checkmat [,3]
        
        
        if (stingy)
        {
          ## Scale errors
          scalederrors <- errors / sd (Y)
          ## Select predictor sets leading to a scaled error that does not exceed
          ## the minimum scaled error + x
          lowerrors <- scalederrors <= min (scalederrors) + 0.025
          ## Remove models with higher error
          comps [!lowerrors] <- 9999
          ## Get the models with the lowest nlv where error minimizes
          these <- comps == min (comps)
          scalederrors [!these] <- 9999
          ## Among these models, get the first one with the lowest error
          thisone <- which.min (scalederrors)
        }
        else ## several solutions and not stingy
        {
          ## Get the models with the lowest error
          these <- errors.min == min (errors.min)          
          ## Among these models, get the one with the lowest nlv
          comps [!these] <- 9999
          thisone <- which.min (comps)
        }        
        
        result <- list (selmat [thisone,], compn [thisone])
      }
    
    }

    return (result)
  }

  ## VIP.R: Implementation of VIP (variable importance in projection)(*) for the
  ## `pls' package.
  ## Copyright © 2006,2007 Bjørn-Helge Mevik
  ## This program is free software; you can redistribute it and/or modify
  ## it under the terms of the GNU General Public License version 2 as
  ## published by the Free Software Foundation.
  ##
  ## This program is distributed in the hope that it will be useful,
  ## but WITHOUT ANY WARRANTY; without even the implied warranty of
  ## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ## GNU General Public License for more details.

  ## A copy of the GPL text is available here:
  ## http://www.gnu.org/licenses/gpl-2.0.txt

  ## Contact info:
  ## Bjørn-Helge Mevik
  ## bhx6@mevik.net
  ## Rødtvetvien 20
  ## N-0955 Oslo
  ## Norway

  ## (*) As described in Chong, Il-Gyo & Jun, Chi-Hyuck, 2005, Performance of
  ## some variable selection methods when multicollinearity is present,
  ## Chemometrics and Intelligent Laboratory Systems 78, 103--112.

  ## VIP returns all VIP values for all variables and all number of components,
  ## as a ncomp x nvars matrix.

  VIP <- function (object)
  {
      SS <- c(object$Yloadings) ^ 2 * colSums (object$scores^2)
      Wnorm2 <- colSums (object$loading.weights ^ 2)
      SSW <- sweep (object$loading.weights ^ 2, 2, SS / Wnorm2, "*")
      sqrt (nrow (SSW) * apply (SSW, 1, cumsum) / cumsum (SS))
  }  

  ## -- END FUNCTIONS ------------------------------------------------------- ##

  ## At this point, a loop begins that stores model results after predictor
  ## selection and turns back here.

  counter <- 1
  loop <- TRUE

  while (loop == TRUE)
  {

    ## --- Variables from previous run and update of X ---------------------- ##

    if (counter > 1)
    {
      ## Previous selection of predictors
      prev.selection <-
        as.logical (unlist (tmp [[6 + counter]][length (model) + 5]))
      ## Previous RMSE
      prev.rmse <- tmp [[6 + counter]] [[length (model) + 1]]
      ## Previous R2
      prev.r2 <- tmp [[6 + counter]] [[length (model) + 2]]
      ## Previous nlv
      prev.nlv <- tmp [[6 + counter]] [[length (model) + 3]]
      ## Previous croperror
      prev.croperror <- tmp [[6 + counter]] [[length (model) + 4]]
      ## Construct new X
      subX <- X [,prev.selection] ##selection
    }
    else subX <- X ## else rename X nevertheless

    ## original predictor names
    subX.names <- colnames (subX)

    ## --- Data preparation ------------------------------------------------- ##

    ## X preprocessing
    if (!is.na (prep)) subX <- prepro (subX, prep)

    ## --- Run PLSR --------------------------------------------------------- ##

    model <- plsr (Y ~ subX, jackknife = TRUE, scale = scaling,
      validation = val, method = 'oscorespls')

    ## Model parameters
    rmse <- getrmse (model)
    nlv <- getbest (rmse)
    r2 <- getr2 (model)
    r2.cal <- r2 [1,nlv]
    r2.val <- r2 [2,nlv]
    reg.coef <- as.vector (coef (model, ncomp = nlv))
    jt <- as.vector (suppressWarnings (jack.test (model, ncomp = nlv)$pvalues))
    vip <- VIP (model) [nlv,]
    vip <- (vip - min (vip)) / max (vip - min (vip))

    ## ---- Reporting ------ ##
    if (verbose)
    {
      cat (counter,
        paste ('  Pred:', ncol (subX)),
        paste ('  LV:', nlv),
        paste ('  R2val:', round (r2.val, 6)),
        paste ('  RMSEval:', round (rmse [2, nlv], 6)))
      if (counter > 1) cat (paste ('  Criterion:',  res [[2]]))
      cat ('\n')
      flush.console ()
    }
     
    ## ---------- Selection of predictors for subsequent runs --------------- ##

    selcode <- NA

    if (backselect == TRUE)
    {
      res <- tryaround (subX, reg.coef, vip, jt)
      selcode <- unlist (res [[1]])    
    }
    
    if (counter > 1)
    {
      new.selection <- prev.selection      
      new.selection [prev.selection == TRUE] <- selcode
      selcode <- new.selection
    }

    ## --- Result object ---------------------------------------------------- ##

    ## --- Construct tables of RMSE, R2 and a vector of nlv

    ## Compute RMSE and R2
    ## ... in case of large nlv
    if (nlv > maxnlv) 
    {
      rmse.crop <- rmse [2, maxnlv]
      r2.crop <- r2 [2, maxnlv]
    }
    else 
    {
      rmse.crop <- rmse [2, nlv]
      r2.crop <- r2.val
    }  
    metarmse <- c(rmse [1, nlv], rmse [2, nlv], rmse.crop)
    names (metarmse) <- c('RMSEcal', 'RMSEval', 'RMSEcrop')
    metar2 <- c(r2.cal, r2.val, r2.crop)
    names (metar2) <- c('R2cal', 'R2val', 'R2crop')
    metanlv <- nlv

    if (counter > 1)
    {
      ## append RMSE to the existing table of RMSE's
      iter.1 <- prev.rmse
      metarmse <- cbind (iter.1, metarmse)
      colnames (metarmse) [counter] <- paste ('iter.', counter, sep = '')

      ## same for R2
      iter.1 <- prev.r2
      metar2 <- cbind (iter.1, metar2)
      colnames (metar2) [counter] <- paste ('iter.', counter, sep = '')

      ## same for nlv
      iter.1 <- prev.nlv
      metanlv <- c (iter.1, nlv)
      names (metanlv) [counter] <- paste ('iter.', counter, sep = '')
    }

    ## --- First time here and in backselection mode? Prepare result object

    if (counter == 1)
    {
      ## First add information about the best model:
      tmp <- list ( ## will be added as 'metapls' object
        'preprocessing' = prep,       ## [[1]]
        'scaling' = scaling,          ## [[2]]
        'RMSE' = NA,                  ## [[3]]
        'R2' = NA,                    ## [[4]]
        'nlv' = NA,                   ## [[5]]
        'best.model' = NA,            ## [[6]]
        'predictors' = NA)            ## [[7]]
        ## This is followed by details about single iterations (see below)

      predictors <- rep (TRUE, ncol (subX)) ## All predictors have been used
      names (metanlv) [1] <- 'iter.1'

    }
    else predictors <- prev.selection

    ## List element containing the iteration results

    index <- 7 + counter
    lth <- length (model)

    tmp [[index]] <- model
    names (tmp)[index] <- paste ('iter.', counter, sep = '')

    ## Element [[1]]: RMSE
    tmp [[index]][[lth + 1]] <- metarmse
    names (tmp [[index]])[lth + 1] <- 'RMSE'

    ## Element [[2]]: R2
    tmp [[index]][[lth + 2]] <- metar2
    names (tmp [[index]])[lth + 2] <- 'R2'

    ## Element [[3]]: Number of latent vectors (nlv)
    tmp [[index]][[lth + 3]] <- metanlv
    names (tmp [[index]])[lth + 3] <- 'nlv'

    ## Element [[4]]: Predictor selection for current run
    tmp [[index]][[lth + 4]] <- predictors
    names (tmp [[index]])[lth + 4] <- 'predictors'
    
    ## Element [[5]]: Predictor selection for next run
    tmp [[index]][[lth + 5]] <- selcode
    names (tmp [[index]])[lth + 5] <- 'selection'

    ## Continue?
    if (counter > 1)
    {
      lth <- ncol (metar2)
      
      ## Continue if the last model is not worse than "limit"
      limit <- max (metar2 [2,]) - 0.05 
      if (limit > 0 & metar2 [2,lth] > limit) loop <- TRUE
      
      else
      {
        if (metanlv [lth] < metanlv [lth - 1]) loop <- TRUE
        else loop <- FALSE
      }
    }

    ## exit if backselect = FALSE
    if (backselect == FALSE) loop <- FALSE
    
    ## exit if less than 10 predictors have been selected
    else if (sum (selcode) < 10) loop <- FALSE

    ## Final action in case of no return
    if (loop == FALSE)
    {
      ## Parameters relating to the best model
      ## Case of no selection but backselect == TRUE or numeric
      if (counter == 1)
      {
        best <- 1
        tmp [[7]] <- predictors                           ## selection
        tmp [[6]] <- best                                 ## best model
        tmp [[5]] <- tmp [[8]][length (model) + 3]        ## nlv
        tmp [[4]] <- metar2                               ## R2
        tmp [[3]] <- metarmse                             ## RMSE
      }
      
      ## Case of prev. selection and backselect == TRUE or numeric
      else
      {
        ## ---- Final model selection out of all iterations:

        meta.error <- unlist (metarmse [2,])
        meta.crop <- unlist (metarmse [3,])

        ## this is the first model with the global, minimum RMSEval:
        errorsel <- which (meta.error == min (meta.error))[1]

        if (stingy)
        {
          ## nlv in this model:
          errorsel.nlv <- tmp [[7 + errorsel]][[length (model) + 3]][errorsel]
          ## If the model with the lowest RMSEval has a nlv <= maxnlv take that:
          if (errorsel.nlv <= maxnlv) best <- errorsel
          ## otherwise take the model with the lowest RMSE at maxnlv:
          else best <- which (meta.crop == min (meta.crop))[1]
        }
        else best <- errorsel
        names (best) <- NULL
      }
  
      ## Model parameters
      opt.lv <- metanlv [best] ## sequence of lv
      names (opt.lv) <- NULL
      preds <- unlist (tmp [[7 + best]][length (model) + 4]) ## last predictors
      names (preds) <- NULL      
      
      tmp [[7]] <- preds ## last predictors
      tmp [[6]] <- best  ## selected model
      tmp [[5]] <- opt.lv ## nlv

      if (is.null (dim (metarmse) [2]))
      {
        opt.rmse <- metarmse [2]
        opt.r2 <- metar2 [2]
        tmp [[4]] <- opt.r2 ## R2
        tmp [[3]] <- opt.rmse ## RMSE
        catrmse <- format (opt.rmse, digits = 6)
        catrmse.crop <- format (metarmse [3], digits = 6)
        catr2 <- format (opt.r2, digits = 6)
        catr2.crop <- format (metar2 [3], digits = 6)
      }
      else
      {
        opt.rmse <- metarmse [2, best]
        opt.r2 <- metar2 [2, best]
        tmp [[4]] <- metar2 [,best]                       ## R2
        tmp [[3]] <- metarmse  [,best]                    ## RMSE
        catrmse <- format (opt.rmse, digits = 6)
        catrmse.crop <- format (metarmse [3, best], digits = 6)
        catr2 <- format (opt.r2, digits = 6)
        catr2.crop <- format (metar2 [3, best], digits = 6)
      }

      ## The result object consists of three parts: 
      ## (1) The first part represents the "best" model object and has the 
      ##     structure of an mvr object in package pls.
      ## (2) The second is a list object called metapls and contains desriptive
      ##     statistics based on the selected number of latent vectors,
      ##     predictors used, summarizing statistics from the other 
      ##     iterations and the entire, unreduced set of predictors
      ## (3) The third consists of other iterations and can be used to produce 
      ##     a new autopls object based on another iteration.
      
      ## Get "best" model and put it in part 1
      if (scaling) part1 <- tmp [[7 + best]][-c(21:23,25)]      
      else part1 <- tmp [[7 + best]][-c(20:22,24)]      
      
      ## Construct part 2 ($metapls)  
      part2 <- list (
        current.iter = best,
        autopls.iter = best,
        current.lv = opt.lv,       
        autopls.lv = opt.lv,
        lv.history = tmp [[length (tmp)]] $nlv,
        rmse.history = tmp [[length (tmp)]] $RMSE,
        r2.history = tmp [[length (tmp)]] $R2,
        X = X,
        Y = Y,
        preprocessing = prep,
        scaling = scaling,
        call = sys.call ()
        )
      
      ## Get data for part 3 ($iterations)
      part3 <- tmp [8:length(tmp)]

      ## In order to save memory, remove original data from all iterations
      ## This can be reconstructed using $X and $Y ($metapls) and $predictors
      for (i in 1:length(part3))
      { 
        if (scaling) part3 [[i]] <- part3 [[i]] [-c(21:23,25)]
        else  part3 [[i]] <- part3 [[i]] [-c(20:22,24)]
        part3 [[i]] $model$Y <- NULL
        part3 [[i]] $model$subX <- NULL
      }
      
      ## Put parts together
      result <- c(part1)
      result [[length (result) + 1]] <- part2
      names (result) [length (result)] <- 'metapls'
      result [[length (result) + 1]] <- part3
      names (result) [length (result)] <- 'iterations'
    }
    
    ## Update counter
    counter <- counter + 1
  }
  
  class (result) <- 'autopls'

  ## Reporting
  if (verbose) print (result)  
  invisible (result)
}
