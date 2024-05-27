#'Probability distribution function of overall effect by M-estimator based approach
#'
#'@param qval value of overall effect
#'@param yi vector of length K with the observed effect sizes or outcomes.
#'@param vi vector of length K with the corresponding sample variances.
#'@param theta value of true overall effect that define 0 in default.
#'@param tau2 value of estimated between-study variance.
#'@param tau2b vector of candidate values of M-estimator based parameter. If the default is insufficient, the search area can be increased by increasing the upper limit or narrowing the interval.
#'@param n.tau2 value of number of sampling for estimates of between-study variance.
#'@param method character string to specify the estimating method of between-study variance. By default (method="DL"), DerSimonian-Laird estimator is obtained. When method="SJ", Sidik-Jonkman estimator is obtained. When method="HKSJ-DL", DerSimonian-Laird method with HKSJ approach is obtained. When method="HKSJ-SJ", Sidik-Jonkman method with HKSJ approach is obtained.
#'
#'@return Probability of qval and M-estimator based estimates of location parameter.
#'
#'@export
ptheta_with_M_estimator <- function(qval, yi, vi, theta=0, tau2, tau2b=c(0.001,1:300/10), n.tau2=10000, method="DL"){
  rho <- function(u, func=1){
    if(func==1){
      return(abs(u))
    }else if(func==2){
      return(u^2)
    }else{
      return(u)
    }
  }
  estfun <- function(tau2b, x, func=1){
    n.tau2b <- length(tau2b)
    R <- array(10000, dim=n.tau2b)
    for(i.tau2m in 1:n.tau2b){
      R_tau2b <- array(10000, dim=n.tau2b)
      for(i.tau2b in 1:n.tau2b){
        R_tau2b[i.tau2b] <- mean(rho(x[,i.tau2m] - tau2b[i.tau2b], func))
      }
      R[i.tau2m] <- tau2b[which.min(R_tau2b)]
    }
    return(R)
  }
  K <- length(vi)
  if(method=="DL"){
    x <- rand_tau2dl_tau2(n=n.tau2, K=K, tau2=tau2b, sig2k=vi)
    R <- estfun(tau2b, x, func=1)
    tau2m <- tau2b[which.min(abs(R-tau2))]
    x.m <- x[,which.min(abs(R-tau2))]
    v2 <- numeric(n.tau2)
    for(k in 1:K){
      v2 <- v2 + 1/(vi[k]+x.m)
    }
    v2 <- 1/v2
    qy <- sapply(qval, function(yval){mean(stats::pnorm(yval, mean=theta, sd=sqrt(v2)))})
  }else if(method=="SJ"){
    x <- sapply(tau2b, function(x){ x/(K-1)*stats::rchisq(n.tau2, df=K-1) })
    R <- estfun(tau2b, x, func=1)
    vk_m <- sapply(1:K, function(k){1 + vi[k]/tau2b})
    theta_vm <- apply(yi/t(vk_m), 2, sum) / apply(1/vk_m, 1, sum)
    tau2sj_m <- apply( (t(array(yi, dim=c(K,length(tau2b))))-theta_vm)^2/vk_m, 1, sum) / (K-1)
    tau2m_sj_num <- which.min(abs(R-tau2sj_m))
    tau2m <- tau2b[tau2m_sj_num]
    x.m <- x[,tau2m_sj_num]
    v2 <- numeric(n.tau2)
    for(k in 1:K){
      v2 <- v2 + 1/(vi[k]+x.m)
    }
    v2 <- 1/v2
    qy <- sapply(qval, function(yval){mean(stats::pnorm(yval, mean=0, sd=sqrt(v2)))})
  }else if(method=="HKSJ-DL"){
    x <- rand_tau2dl_tau2(n=n.tau2, K=K, tau2=tau2b, sig2k=vi)
    R <- estfun(tau2b, x, func=1)
    tau2m <- tau2b[which.min(abs(R-tau2))]
    x.m <- x[,which.min(abs(R-tau2))]
    wk_m <- sapply(1:K, function(k){1/(vi[k]+x.m)})
    theta_hksj_m <- apply(t(wk_m)*yi, 2, sum) / apply(wk_m, 1, sum)
    v2_hksj_m <- apply(wk_m*(t(array(yi,dim=c(K, n.tau2)))-theta_hksj_m)^2, 1, sum) / ((K-1)*apply(wk_m, 1, sum))
    qy <- sapply(qval, function(yval){mean(stats::pt(q=(yval-theta_hksj_m)/sqrt(v2_hksj_m), df=K-1))})
  }else if(method=="HKSJ-SJ"){
    x <- sapply(tau2b, function(x){ x/(K-1)*stats::rchisq(n.tau2, df=K-1) })
    R <- estfun(tau2b, x, func=1)
    vk_m <- sapply(1:K, function(k){1 + vi[k]/tau2b})
    theta_vm <- apply(yi/t(vk_m), 2, sum) / apply(1/vk_m, 1, sum)
    tau2sj_m <- apply( (t(array(yi, dim=c(K,length(tau2b))))-theta_vm)^2/vk_m, 1, sum) / (K-1)
    tau2m_sj_num <- which.min(abs(R-tau2sj_m))
    tau2m <- tau2b[tau2m_sj_num]
    x.m <- x[,tau2m_sj_num]
    wk_m <- sapply(1:K, function(k){1/(vi[k]+x.m)})
    theta_hksj_m <- apply(t(wk_m)*yi, 2, sum) / apply(wk_m, 1, sum)
    v2_hksj_m <- apply(wk_m*(t(array(yi,dim=c(K, n.tau2)))-theta_hksj_m)^2, 1, sum) / ((K-1)*apply(wk_m, 1, sum))
    qy <- sapply(qval, function(yval){mean(stats::pt(q=(yval-theta_hksj_m)/sqrt(v2_hksj_m), df=K-1))})
  }else{
    stop("Method should be DL, SJ, HKSJ-DL or HKSJ-SJ")
  }

  return(list(qy=qy, tau2m=tau2m))
}
