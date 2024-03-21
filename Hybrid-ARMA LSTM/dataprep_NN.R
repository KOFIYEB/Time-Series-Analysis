dataprep_nn <- function(y, ntrain = NULL, p.max, q.max, fixed = FALSE, criterion = "bic") {
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
  for(i in 1:ntest) {
    yi <- y[(ntrain + i - 1):(ntrain + i - p)] - mutrain
    resi <- resid[(ntrain + i - 1):(ntrain + i - q)]
    ytest_arma[i] <- sum(c_arma * c(yi, resi)) + mutrain
    resid <- c(resid, y[ntrain + i] - ytest_arma[i])
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
  
  return(list(Xtrain = Xtrain,
              Xtest = Xtest,
              Ytrain = Ytrain,
              Ytest = Ytest,
              Y0test = ytest,
              Ytrain_arma = yarma_fit,
              Ytest_arma = ytest_arma,
              Ymin = ymin,
              Ydiff = dminmax)
         )
}