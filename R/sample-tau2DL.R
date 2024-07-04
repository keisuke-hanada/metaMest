#'Sampling from exact distribution of DerSimonian-Laird between-study variance estimator
#'
#'@param n number of observations.
#'@param K number of studies.
#'@param tau2 value of between-study variance parameter. vector is available.
#'@param sig2k vector of within-study variance in each study.
#'
#'@return sampling values from exact distribution of DerSimonian-Laird between-study variance estimator
#'
#'@export
rand_tau2dl_tau2 <- function(n, K, tau2, sig2k){
  s1 <- sum(1/sig2k)
  s2 <- sum(1/sig2k^2)
  chi_sim <- array(0, dim=c(K,n))
  for(i in 1:K){
    chi_sim[i,] <- stats::rchisq(n=n, df=1)
  }
  W <- diag(x=1/sig2k,nrow=K,ncol=K)

  x_ftau20 <- array(0, dim=c(length(tau2), n))
  for(i.b in 1:length(tau2)){
    V <- U <- array(0, dim=c(K,K))
    lambda <- array(0, dim=c(K))
    for(i in 1:K){
      for(j in 1:K){
        V[i,j] <- -1/s1 + (s2/s1^2 - (1/sig2k[i]+1/sig2k[j])/s1)*tau2[i.b]
      }
      V[i,i] <- (sig2k[i] - 1/s1) + (1 + s2/s1^2 - 2/sig2k[i]/s1)*tau2[i.b]
    }
    U <- V %*% W
    lambda <- eigen(U)$value
    Q_chis <- c(lambda %*% chi_sim)

    x_ftau20[i.b,] <- (Q_chis - (K-1)) / (s1-s2/s1)
  }
  x_ftau2 <- x_ftau20
  x_ftau2[x_ftau2<0] <- 0
  return(t(x_ftau2))
}
