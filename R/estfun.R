#'Estimation function of between-study variance
#'
#'@param tau2m M-estimator of between-study variance.
#'@param tau2h estimator of between-study variance by setting method. The default is "DL" method.
#'@param K number of studies.
#'@param vi vector of length K with the corresponding sample variances.
#'@param method character string to specify the estimating method of between-study variance. By default (method="DL"), DerSimonian-Laird estimator is obtained. When method="SJ", Sidik-Jonkman estimator is obtained.
#'@param n_ftau2 value of number of sampling for estimates of between-study variance. The default is 1000.
#'@param func If 1 (default), the convex function is absolute function. If 2, the convex function is square function.
#'@param seed seed value of random numbers.
#'
#'@return value of estimation function for between-study variance.
#'
#'@export
estfun <- function(tau2m, tau2h, K, vi, method="DL", n_ftau2=1000, func=1, seed=1){
set.seed(seed)
if(method == "DL"){
  x <- rand_tau2dl_tau2(n=n_ftau2, K=K, tau2=tau2m, sig2k=vi)
}else if(method == "SJ"){
  x <- tau2m / (K-1) * stats::rchisq(n_ftau2, K-1)
}
x <- sapply(x, function(y) max(y,1e-100))
Rx <- sum(rho(log(x)-log(tau2h), func))
return(Rx)
}
