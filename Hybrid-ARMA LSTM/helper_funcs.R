dataprep_nn_arma <- function(y, ntrain = NULL, p.max = 0, q.max = 0, fixed = FALSE, criterion = "bic") {
  n <- length(y)
  #ntrain = floor(n * ntrain + 0.5)
  ntest = n - ntrain
  y_arma <- y[1:ntrain]
  ytest <- y[(ntrain + 1):(ntrain + ntest)]
  mutrain <- mean(y_arma)
  
  if (isFALSE(fixed)) {
    BICmat <- smoots::critMatrix(y_arma - mutrain, p.max = p.max, q.max = q.max, criterion = criterion)
    
    pq <- suppressMessages(smoots::optOrd(BICmat))
    p <- unname(pq[1])
    q <- unname(pq[2])
  }
  if(isTRUE(fixed)) {
    p <- p.max
    q <- q.max
  }
  
  m_arma <- arima(y_arma - mutrain, order=c(p, 0, q))
  c_arma <- m_arma$coef[1:(p + q)]  
  
  resid <- m_arma$resid
  yarma_fit <- y_arma - resid
  ytest_arma <- ytest * 0
  
  if (q == 0) {
    for(i in 1:ntest) {
      yi <- y[(ntrain + i - 1):(ntrain + i - p)] - mutrain
      ytest_arma[i] <- sum(c_arma * c(yi)) + mutrain
    }
  }
  
  if (q != 0) {
    for(i in 1:ntest) {
      yi <- y[(ntrain + i - 1):(ntrain + i - p)] - mutrain
      resi <- resid[(ntrain + i - 1):(ntrain + i - q)]
      ytest_arma[i] <- sum(c_arma * c(yi, resi)) + mutrain
      resid <- c(resid, y[ntrain + i] - ytest_arma[i])
    }
  }
  
  ### creating input variables for NN
  y_nn <- y[(p + 1):n]
  datnorm <- minmax(y_nn)
  ystd = datnorm$xstd
  ymin = datnorm$xmin
  dminmax = datnorm$dminmax
  
  armaord <- p + q
  if (armaord == 0) armaord <- 1
  X <- matrix(0, n - p, armaord)
  
  if (p != 0) {
    for (j in 1:p) {
      X[, j] <- minmax(y[(p - j + 1):(n - j)])$xstd
    }
  }
  
  if (q != 0) {
    res_std <- minmax(resid)$xstd
    for(k in 1:q){
      X[, p + k] <- res_std[(p + 1 - k):(n - k)]
    }
  }
  
  Xtrain <- X[1:(ntrain - p),]
  Ytrain <- ystd[1:(ntrain - p)]
  Xtest <- X[(ntrain - p + 1):(n - p), ]
  Ytest <- ystd[(ntrain - p + 1):(n - p)]
  
  # workaround as `mlp()`only accepts matrices
  if (p == 1) {
    Xtrain = as.matrix(Xtrain)
    Xtest = as.matrix(Xtest)
  }
  
  return(list(Xtrain = Xtrain,
              Xtest = Xtest,
              Ytrain = Ytrain,
              Ytest = Ytest,
              Y0test = ytest,
              Ytrain_arma = yarma_fit,
              Ytest_arma = ytest_arma,
              Ymin = ymin,
              Ydiff = dminmax,
              m_arma = m_arma)
  )
}

dataprep_nn_garch <- function(y, ntrain = NULL, p.max = 0, q.max = 0, fixed = FALSE, criterion = "bic") {
  n <- length(y)

  #ntrain = floor(n * ntrain + 0.5)
  ntest = n - ntrain
  ytest <- y[(ntrain + 1):n]^2
  
  if (isFALSE(fixed)) {
    BIC_GARCH = matrix(0, p.max, q.max + 1)
    for (i in 1:p.max) {
      for (j in 0:q.max) {
        spec <- rugarch::ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(i, j)),
                                    mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
                                    distribution.model = "std")
        GARCHij = rugarch::ugarchfit(spec = spec, data = y, out.sample = ntest)
        BIC_GARCH[i, j + 1] = rugarch::infocriteria(GARCHij)["Bayes", ]
      }
    }
    p = which(BIC_GARCH == min(BIC_GARCH), arr.ind = TRUE)[1]
    q = which(BIC_GARCH == min(BIC_GARCH), arr.ind = TRUE)[2] - 1
  }
  if(isTRUE(fixed)) {
    p <- p.max
    q <- q.max
  }
  spec <- rugarch::ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(p, q)),
                              mean.model = list(armaOrder = c(0, 0)),
                              distribution.model = "std")
  
  m_garch = rugarch::ugarchfit(spec = spec, data = y, out.sample = ntest)
  GARCH_fc = rugarch::ugarchforecast(m_garch, n.ahead = 1, n.roll = ntest - 1)
  sig2.tr = as.numeric(rugarch::sigma(m_garch))^2
  sig2.te = as.numeric(rugarch::sigma(GARCH_fc))^2
  
  ### creating input variables for NN
  y_nn <- y[(p + 1):n]^2
  datnorm <- minmax(y_nn)
  ystd = datnorm$xstd
  ymin = datnorm$xmin
  dminmax = datnorm$dminmax
  
  garchord <- p + q
  if (garchord == 0) garchord <- 1
  X <- matrix(0, n - p, garchord)
  
  if (p != 0) {
    for (j in 1:p) {
      X[, j] <- minmax(y[(p - j + 1):(n - j)])$xstd
    }
  }
  
  if (q != 0) {
    resid = c(sig2.tr, sig2.te)
    res_std <- minmax(resid)$xstd
    for(k in 1:q){
      X[, p + k] <- res_std[(p + 1 - k):(n - k)]
    }
  }
  
  Xtrain <- X[1:(ntrain - p),]
  Ytrain <- ystd[1:(ntrain - p)]
  Xtest <- X[(ntrain - p + 1):(n - p), ]
  Ytest <- ystd[(ntrain - p + 1):(n - p)]
  
  # workaround as `mlp()`only accepts matrices
  if (p == 1) {
    Xtrain = as.matrix(Xtrain)
    Xtest = as.matrix(Xtest)
  }
  
  return(list(Xtrain = Xtrain,
              Xtest = Xtest,
              Ytrain = Ytrain,
              Ytest = Ytest,
              Y0test = ytest,
              Ytrain_garch= sig2.tr,
              Ytest_garch = sig2.te,
              Ymin = ymin,
              Ydiff = dminmax,
              m_garch = m_garch)
  )
}

ASE_calc = function(y, yhat) {
  ASE = mean((y -yhat)^2)
  return(ASE)
}

AAE_calc = function(y, yhat) {
  AAE = mean(abs(y - yhat))
  return(AAE)
}

minmax <- function(x) {
  xmin = min(x)
  xmax = max(x)
  dminmax = xmax - xmin
  xstd = (x - xmin) / dminmax
  return(list(xstd = xstd,
              dminmax = dminmax,
              xmin = xmin)
  )
}