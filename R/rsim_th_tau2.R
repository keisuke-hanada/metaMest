#'Simultaneous distribution of overall treatment effect and between-study variance
#'
#'@param n number of random variables.
#'@param theta value of true overall effect that define 0 in default.
#'@param tau2 between-study variance.
#'@param K number of studies.
#'@param vi vector of length K with the corresponding sample variances.
#'@param method character string to specify the estimating method of between-study variance. By default (method="DL"), DerSimonian-Laird estimator is obtained. When method="SJ", Sidik-Jonkman estimator is obtained.
#'@param obj character string to specify the objective. "Estimate", which set to default, is overall treatment effect. "Predict" is predictive treatment effect.
#'
#'@return value of simultaneous distribution of overall treatment effect and between-study variance.
#'
#'@export
rsim_th_tau2 <- function(n, theta=0, tau2, K, vi, method="DL", obj="Estimate"){
  if(method == "DL"){
    x <- rand_tau2dl_tau2(n=n, K=K, tau2=tau2, sig2k=vi)
  }else if(method == "SJ"){
    x <- sapply(tau2, function(tau2) tau2/(K-1) * rchisq(n=n, df=K-1))
  }
  v2 <- array(1 / sapply(x, function(xval) sum(1/(vi+xval))), dim=c(n,length(tau2)))
  if(obj == "Estimate"){
    y <- array(rnorm(n=n*length(tau2), mean=theta, sd=sqrt(v2)), dim=c(n,length(tau2)))
  }else if(obj == "Predict"){
    y <- array(rnorm(n=n*length(tau2), mean=theta, sd=sqrt(v2+x)), dim=c(n,length(tau2)))
  }
  names.x <- paste("tau2=", tau2, sep="")
  dimnames(x) <- dimnames(y) <- list(NULL, names.x)

  return(list(theta=y, tau2=x))
}
