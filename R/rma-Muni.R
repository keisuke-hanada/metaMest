#'Random-effects meta-analysis with M-estimator
#'
#'@param yi vector of length K with the observed effect sizes or outcomes.
#'@param vi vector of length K with the corresponding sample variances.
#'@param method character string to specify the estimating method of between-study variance. By default (method="DL"), DerSimonian-Laird estimator is obtained. When method="SJ", Sidik-Jonkman estimator is obtained.
#'@param HKSJ If TRUE, HKSJ approach is applied. HKSJ=FALSE is default.
#'@param alpha significance level of confidence intervals, and simultaneous plot.
#'@param n.rand sampling number of random values for distribution of between-study variance.
#'@param tau2.max upper bound for M-estimate of tau2.
#'@param theta.min lower bound for M-estimate of theta.
#'@param n.grid number of grid for contour plot of treatment effect and between-study variance.
#'@param obj character string to specify the objective. "Estimate", which set to default, is overall treatment effect. "Predict" is predictive treatment effect.
#'@param contour_val vector of the thresholds for contour plot. By default, the line of (0.01, 0.05, 0.2, 0.4, 0.6, 0.8) is plotted.
#'@param plot If TRUE, the contour plot is created.
#'
#'@return overall results for M-estimated random-effects meta-analysis with logarithmic absolute convex function.
#'
#'
#'@export
rma.Muni <- function(yi, vi, method="DL", HKSJ=FALSE, alpha=0.05, n.rand=10000, tau2.max=100,
                     epsilon=1e-10, theta.min=-100, theta.max=100, n.grid=300, obj="Estimate",
                     contour_val=c(0.01, 0.05, 0.2, 0.4, 0.6, 0.8), plot=TRUE){
  K <- length(yi)
  s1 <- sum(1/vi); s2 <- sum(1/vi^2)

  ### each method estimate
  if(method=="DL"){
    wk <- 1/vi
    thetabar0 <- sum(yi*wk)/sum(wk)
    Q <- sum(wk*(yi-thetabar0)^2)
    tau2u <- (Q-(K-1))/(sum(wk)-sum(wk^2)/sum(wk))
    tau2h <- max(0,tau2u)
  }else if(method=="SJ"){
    thetabar <- mean(yi)
    vk <- 1+ K*vi/sum((yi - thetabar)^2)
    theta_v <- sum(yi/vk) / sum(1/vk)
    tau2h <- sum((yi-theta_v)^2/vk) / (K-1)
  }
  wkh <- 1/(vi + tau2h)
  thetah <- sum(wkh*yi) / sum(wkh)
  v2h <- 1/sum(wkh)


  ### M-estimator
  r_chisq <- array(rchisq(K*n.rand, df=1), c(n.rand,K))
  tau2m <- optimize(rfunc, tau2h=tau2h, vi=vi, r_chisq=r_chisq, method=method, interval=c(0,tau2.max))$minimum


  ### CI of between-study variance
  eigen_tau2 <- dtau2_eigen(tau2=tau2h, vi=vi, method=method)
  if(method=="DL"){
    tau2_rand <- (r_chisq %*% eigen_tau2 - (K-1))/(s1-s2/s1)
  }else if (method=="SJ"){
    tau2_rand <- r_chisq %*% eigen_tau2
  }
  tau2_rand[tau2_rand<=epsilon] <- epsilon
  cil.tau2 <- quantile(tau2_rand, probs=alpha/2)
  ciu.tau2 <- quantile(tau2_rand, probs=1-alpha/2)


  ### CI of overall treatment effect
  if(!HKSJ){

    v2m <- numeric(n.rand)
    for(k in 1:K){
      v2m <- v2m + 1/(vi[k]+tau2_rand)
    }
    v2m <- 1/v2m

    ltheta <- function(theta, thetah, v2m, level){
      mean(pnorm(thetah, theta, sqrt(v2m))) - level
    }

    cil <- uniroot(ltheta, thetah=thetah, v2m=v2m, level=1-alpha/2, interval=c(theta.min,theta.max))$root
    ciu <- uniroot(ltheta, thetah=thetah, v2m=v2m, level=alpha/2, interval=c(theta.min,theta.max))$root

  }else{

    t0val <- t1val <- array(0, n.rand)
    for(k in 1:K){
      t0val <- t0val + (yi[k]/(vi[k]+c(tau2_rand)))
      t1val <- t1val + (1/(vi[k]+c(tau2_rand)))
    }
    theta_hksj <- t0val / t1val

    ## v2hksj
    v2hksj <- v2hksj0 <- array(0, n.rand)
    for(k in 1:K){
      v2hksj0 <- v2hksj0 + (yi[k]-theta_hksj)^2/(vi[k]+c(tau2_rand))
    }
    v2hksj <- v2hksj0 / ((K-1)*t1val)

    ## CI of HKSJ-m
    pt_hksj <- function(q, theta, v2, K, level){
      tval <- (q-theta)/sqrt(v2)
      return(mean(pt(tval, K-1)) - level)
    }
    cil <- uniroot(pt_hksj, theta=theta_hksj, v2=v2hksj, K=K, level=alpha/2, interval=c(theta.min,theta.max))$root
    ciu <- uniroot(pt_hksj, theta=theta_hksj, v2=v2hksj, K=K, level=1-alpha/2, interval=c(theta.min,theta.max))$root

  }




  ### simultaneous distribution plot
  if(plot){
    p.tau2 <- quantile(tau2_rand, probs=c(0.0001, 0.9999))
    if(obj=="Estimate"){
      pmin.theta <- uniroot(ltheta, thetah=thetah, v2m=v2m, level=0.9999, interval=c(theta.min,theta.max))$root
      pmax.theta <- uniroot(ltheta, thetah=thetah, v2m=v2m, level=0.0001, interval=c(theta.min,theta.max))$root
    }else if(obj=="Predict"){
      p.tau2[2] <- 2*p.tau2[2]
      pmin.theta <- uniroot(ltheta, thetah=thetah, v2m=v2m+max(tau2_rand), level=0.9999, interval=c(theta.min,theta.max))$root
      pmax.theta <- uniroot(ltheta, thetah=thetah, v2m=v2m+max(tau2_rand), level=0.0001, interval=c(theta.min,theta.max))$root
    }

    yval <- seq(pmin.theta, pmax.theta, length.out=n.grid)
    xval <- seq(p.tau2[1], p.tau2[2], length.out=n.grid)
    zval <- dsim_th_tau2(y=yval, x=xval, theta=thetah, tau2=tau2h, K=K, vi=vi, method=method, obj=obj)
    df <- expand.grid(x = yval, y = xval)
    df$z <- as.vector(t(zval))

    # threshold
    v2 <- 1 / sapply(tau2_rand, function(xval) sum(1/(vi+xval)))
    if(obj == "Estimate"){
      yrand <- rnorm(n=n.rand, mean=thetah, sd=sqrt(v2))
    }else if(obj == "Predict"){
      yrand <- rnorm(n=n.rand, mean=thetah, sd=sqrt(v2+tau2_rand))
    }

    dxys <- dsim_th_tau2(y=yrand, x=c(tau2_rand), theta=thetah, tau2=tau2h, K=K, vi=vi, method=method, obj=obj)
    z_threshold <- quantile(dxys, probs = alpha)
    z_contour <- quantile(dxys, probs = contour_val)
    z_threshold <- quantile(dxys, probs = alpha)

    # labeling
    contour_plot <- ggplot2::ggplot(df, aes(x = x, y = y, z = z)) +
      ggplot2::geom_contour(breaks = z_contour)
    contour_data <- ggplot2::ggplot_build(contour_plot)$data[[1]]
    label_data <- aggregate(cbind(x,y) ~ level, data = contour_data, FUN = max)
    label_data$x <- thetah
    label_data$z <- label_data$level
    label_data$label <- contour_val

    g <- ggplot2::ggplot(df, aes(x = x, y = y, z = z)) +
      ggplot2::geom_contour(breaks = z_contour) +
      ggplot2::geom_contour(breaks = z_threshold, color = "green", linewidth = 1.2) +
      ggplot2::geom_text(data=label_data, aes(x=x, y=y, label=label), color="black") +
      ggplot2::labs(
        title = paste("Joint density function of treatment effect and between-study variance: ", method, " method", sep=""),
        x = "treatment effect",
        y = "between-study variance"
      )
  }else{
    g <- NULL
  }


  return(
    list(theta=thetah, tau2=tau2h, tau2m=tau2m, ci.lb=cil, ci.ub=ciu,
         cil.tau2=cil.tau2, ciu.tau2=ciu.tau2, plot=g, method=method, HKSJ=HKSJ, alpha=alpha)
  )

}






#'Random-effects meta-analysis with M-estimator
#'
#'@param yi vector of length K with the observed effect sizes or outcomes.
#'@param vi vector of length K with the corresponding sample variances.
#'@param theta value of estimator for the overall treatment effect.
#'@param tau2 value of estimator for the between-study variance.
#'@param alpha significance level of confidence intervals, and simultaneous plot.
#'@param method character string to specify the estimating method of between-study variance. By default (method="DL"), DerSimonian-Laird estimator is obtained. When method="SJ", Sidik-Jonkman estimator is obtained.
#'@param interval search range of confidence interval for overall treatment effect.
#'@param h search grid of confidence interval for overall treatment effect.
#'
#'@return confidence interval for M-estimated random-effects meta-analysis with absolute convex function.
#'
#'
#'@export
rma.Muni.m1 <- function(yi, vi, theta, tau2, alpha=0.05, method="DL", interval=c(-10,10), h=1000){
  qval <- interval[1] + 0:h*(interval[2]-interval[1]) / h
  pvals <- ptheta_with_M_estimator(qval=qval, yi=yi, vi=vi, theta=theta, tau2=tau2, method=method)
  pval_fun <- approxfun(x=qval, y=pvals$qy, rule=2)
  cifun <- function(qval, alpha){
    pval_fun(qval) - alpha
  }
  if(method == "DL" | method == "SJ"){

    v2val <- uniroot(cifun, alpha=1-alpha/2, interval=c(min(qval),max(qval)))$root
    ci <- theta + c(-1,1)*v2val

  }else if(method == "HKSJ-DL" | method == "HKSJ-SJ"){

    cil <- uniroot(cifun, alpha=alpha/2, interval=c(min(qval),max(qval)))$root
    ciu <- uniroot(cifun, alpha=1-alpha/2, interval=c(min(qval),max(qval)))$root
    ci <- c(cil, ciu)

  }

  return(ci)
}





# internal functions

dtau2_eigen <- function(tau2, vi, method="DL"){
  K <- length(vi)
  if(method=="DL"){
    s1 <- sum(1/vi)
    s2 <- sum(1/vi^2)
    W <- diag(x=1/vi,nrow=K,ncol=K)
    V <- U <- array(0, dim=c(K,K))
    lambda <- array(0, dim=c(K))
    for(i in 1:K){
      for(j in 1:K){
        V[i,j] <- -1/s1 + (s2/s1^2 - (1/vi[i]+1/vi[j])/s1)*tau2
      }
      V[i,i] <- (vi[i] - 1/s1) + (1 + s2/s1^2 - 2/vi[i]/s1)*tau2
    }
    U <- V %*% W
    lambda <- eigen(U)$value
    return(lambda)
  }else if(method=="SJ"){
    lambda <- numeric(K) + tau2/(K-1)
    lambda[K] <- 0
    return(lambda)
  }
}

rfunc <- function(tau2, tau2h, vi, r_chisq, method="DL", epsilon=1e-10){
  K <- length(vi)
  s1 <- sum(1/vi); s2 <- sum(1/vi^2)
  eigen_tau2 <- dtau2_eigen(tau2=tau2, vi=vi, method=method)
  if(method=="DL"){
    tau2_rand <- (colSums(t(r_chisq) * eigen_tau2) - (K-1))/(s1-s2/s1)
  }else if(method=="SJ"){
    tau2_rand <- colSums(eigen_tau2 * t(r_chisq))
  }
  tau2_rand <- sapply(tau2_rand, function(x) max(epsilon,x))
  tau2h_u <- max(epsilon, tau2h)
  return(mean(abs(log(tau2_rand) - log(tau2h_u))))
}
