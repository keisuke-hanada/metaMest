#'Simultaneous distribution of overall treatment effect and between-study variance
#'
#'@param yi vector of length K with the observed effect sizes or outcomes.
#'@param vi vector of length K with the corresponding sample variances.
#'@param K number of studies.
#'@param method character string to specify the estimating method of between-study variance. By default (method="DL"), DerSimonian-Laird estimator is obtained. When method="SJ", Sidik-Jonkman estimator is obtained.
#'@param rho the character to specify the convex function. By default (rho="ML"), maximum likelihood is used. When rho="L1", absolute function is used. When rho="L2", square function is used.
#'@param n_ftau2 value of number of sampling for estimates of between-study variance. The default is 1000.
#'@param tau2.min value to be added so that the logarithm does not go to infinity. By default, 1e-10 is used.
#'@param tau2.max maximum search range in minimizing the estimation equation for M-estimation. The default is 100, but a smaller value may be appropriate depending on the between-study variance.
#'@param seed random number seed used within the function. By default, 1234 is used.
#'
#'@return M-estimator of between-study variance by a specific method.
#'
#'@export
est_tau2m <- function(yi, vi, K, method="DL", rho="ML", n_ftau2, tau2.min=1e-10, tau2.max=100, seed=1234){
  Ropt <- function(tau2, K, yi, vi, method="DL", rho="ML", n_ftau2=1000, tau2.min=1e-10, seed=1234){
    set.seed(seed)
    if(method=="DL"){
      wk <- 1/vi
      thetabar0 <- sum(yi*wk) / sum(wk)
      Q <- sum(wk*(yi-thetabar0)^2)
      tau2x <- max((Q-(K-1))/(sum(wk)-sum(wk^2)/sum(wk)), 0)
    }else if(method=="SJ"){
      vk <- 1 + vi/tau2
      theta_v <- sum(yi/vk) / sum(1/vk)
      tau2x <- sum((yi-theta_v)^2/vk) / (K-1)
    }

    if(rho=="ML"){
      val <- - dtau2(x=tau2x, tau2=tau2, K=K, vi=vi, method=method, n_ftau2=n_ftau2)
    }else if(rho=="L1"){
      x <- rsim_th_tau2(n=n_ftau2, theta=0, tau2=tau2, K=K, vi=vi, method=method)$tau2
      val <- sum(abs(log(x+tau2.min) - log(tau2x+tau2.min)))
    }else if(rho=="L2"){
      x <- rsim_th_tau2(n=n_ftau2, theta=0, tau2=tau2, K=K, vi=vi, method=method)$tau2
      val <- sum((log(x+tau2.min) - log(tau2x+tau2.min))^2)
    }
  }

  x <- optimize(Ropt, K=K, yi=yi, vi=vi, method=method, rho=rho, n_ftau2=n_ftau2,
                tau2.min=tau2.min, seed=seed, interval=c(0,tau2.max), maximum=FALSE)
  return(x$minimum)
}
