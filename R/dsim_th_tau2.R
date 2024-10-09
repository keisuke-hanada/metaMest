#'Simultaneous distribution of overall treatment effect and between-study variance
#'
#'@param theta value of true overall effect that define 0 in default.
#'@param tau2 between-study variance.
#'@param K number of studies.
#'@param vi vector of length K with the corresponding sample variances.
#'@param method character string to specify the estimating method of between-study variance. By default (method="DL"), DerSimonian-Laird estimator is obtained. When method="SJ", Sidik-Jonkman estimator is obtained.
#'@param n_ftau2 value of number of sampling for estimates of between-study variance. The default is 1000.
#'
#'@return value of simultaneous distribution of overall treatment effect and between-study variance.
#'
#'@export
dsim_th_tau2 <- function(y, x, theta=0, tau2, K, vi, method="DL", n_ftau2=100000){
  v2 <- 1 / sapply(x, function(xval) sum(1/(vi+xval)))
  dth <- dnorm(y, mean=theta, sd=sqrt(v2))
  if(method == "DL"){
    rtau2 <- rand_tau2dl_tau2(n=n_ftau2, K=K, tau2=tau2, sig2k=vi)
    dtau2 <- density(rtau2)
    dtau2fun <- approxfun(x=dtau2$x, y=dtau2$y)
    dyx <- dth * dtau2fun(x)
  }else if(method == "SJ"){
    dtau2 <- (K-1)/tau2 * dchisq((K-1)*x/tau2, df = K-1)
    dyx <- dth * dtau2
  }
  dyx[is.na(dyx)] <- 0
  return(dyx)
}
