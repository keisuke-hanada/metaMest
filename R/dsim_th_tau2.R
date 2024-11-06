#'Simultaneous distribution of overall treatment effect and between-study variance
#'
#'@param y value of overall effect.
#'@param x value of between-study variance.
#'@param theta value of true overall effect that define 0 in default.
#'@param tau2 between-study variance.
#'@param K number of studies.
#'@param vi vector of length K with the corresponding sample variances.
#'@param method character string to specify the estimating method of between-study variance. By default (method="DL"), DerSimonian-Laird estimator is obtained. When method="SJ", Sidik-Jonkman estimator is obtained.
#'@param n_ftau2 value of number of sampling for estimates of between-study variance. The default is 100000.
#'@param obj character string to specify the objective. "Estimate", which set to default, is overall treatment effect. "Predict" is predictive treatment effect.
#'
#'@return value of simultaneous distribution of overall treatment effect and between-study variance.
#'
#'@export
dsim_th_tau2 <- function(y, x, theta=0, tau2, K, vi, method="DL", n_ftau2=100000,
                         obj = "Estimate"){
  v2 <- 1 / sapply(x, function(xval) sum(1/(vi+xval)))
  if(obj == "Estimate"){
    dth <- sapply(y, function(yval) dnorm(yval, mean=theta, sd=sqrt(v2)))
  }else if(obj == "Predict"){
    dth <- sapply(y, function(yval) dnorm(yval, mean=theta, sd=sqrt(v2+x)))
  }
  dtau2val <- dtau2(x=x, tau2=tau2, K=K, vi=vi, method=method, n_ftau2=n_ftau2)
  dyx <- array(dim=c(length(y), length(x), length(tau2)))
  for(i in 1:length(tau2)){
    dyx[,,i] <- dth * dtau2val[,i]
  }
  dyx[is.na(dyx)] <- 0
  if(length(tau2)==1) dyx <- dyx[,,1]
  return(dyx)
}


#'Density function of between-study variance
#'
#'@param x value of between-study variance.
#'@param tau2 between-study variance.
#'@param K number of studies.
#'@param vi vector of length K with the corresponding sample variances.
#'@param method character string to specify the estimating method of between-study variance. By default (method="DL"), DerSimonian-Laird estimator is obtained. When method="SJ", Sidik-Jonkman estimator is obtained.
#'@param n_ftau2 value of number of sampling for estimates of between-study variance. The default is 100000.
#'@param obj character string to specify the objective. "Estimate", which set to default, is overall treatment effect. "Predict" is predictive treatment effect.
#'
#'@return value of simultaneous distribution of overall treatment effect and between-study variance.
#'
#'@export
dtau2 <- function(x, tau2, K, vi, method="DL", n_ftau2=100000){
  if(method == "DL"){
    rtau2 <- rand_tau2dl_tau2(n=n_ftau2, K=K, tau2=tau2, sig2k=vi)
    dtau2 <- lapply(1:length(tau2), function(i) density(rtau2[,i]))
    dtau2fun <- lapply(dtau2, function(fun) approxfun(x=fun$x, y=fun$y))
    dx <- sapply(dtau2fun, function(fun) fun(x))
  }else if(method == "SJ"){
    dx <- sapply(tau2, function(val) (K-1)/val * dchisq((K-1)*x/val, df = K-1))
  }
  dx[is.na(dx)] <- 0
  dx <- array(dx, dim=c(length(x), length(tau2)))
  return(dx)
}
