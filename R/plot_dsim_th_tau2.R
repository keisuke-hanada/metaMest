#'Plot function for simultaneous distribution of overall treatment effect and between-study variance
#'
#'@param contour_val vector of the thresholds for contour plot. By default, the line of (0.01, 0.05, 0.2, 0.4, 0.6, 0.8) is plotted.
#'@param theta value of true overall effect that define 0 in default.
#'@param tau2 between-study variance.
#'@param K number of studies.
#'@param vi vector of length K with the corresponding sample variances.
#'@param alpha significance level of simultaneous plot.
#'@param theta.limit lower and upper limits for plotting theta.
#'@param tau2.limit lower and upper limits for plotting tau2.
#'@param length.out x-axis and y-axis step widths
#'@param n number of sampling for calculating the thresholds. The default is 100000.
#'@param method character string to specify the estimating method of between-study variance. By default (method="DL"), DerSimonian-Laird estimator is obtained. When method="SJ", Sidik-Jonkman estimator is obtained.
#'
#'@return contour plot of simultaneous distribution of overall treatment effect and between-study variance.
#'
#'
#'@export
plot_dsim_th_tau2 <- function(contour_val=c(0.01, 0.05, 0.2, 0.4, 0.6, 0.8),
                              theta=0, tau2, K, vi, alpha=0.05,
                              theta.limit=c(-10,10), tau2.limit=c(0,10), length.out=1000,
                              n=100000, method="DL"){
  yval <- seq(theta.limit[1], theta.limit[2], length.out = length.out)
  xval <- seq(tau2.limit[1], tau2.limit[2], length.out = length.out)
  zval <- outer(X=yval, Y=xval, dsim_th_tau2, theta=theta, tau2=tau2, K=K, vi=vi, method=method)
  df <- expand.grid(x = yval, y = xval)
  df$z <- as.vector(zval)

  # threshold
  xyrand <- rsim_th_tau2(n=n, theta=theta, tau2=tau2, K=K, vi=vi, method=method)
  dxys <- dsim_th_tau2(y=xyrand$theta, x=xyrand$tau2, theta=theta, tau2=tau2, K=K, vi=vi, method=method)
  z_threshold <- quantile(dxys, probs = alpha)
  z_contour <- quantile(dxys, probs = contour_val)
  z_threshold <- quantile(dxys, probs = alpha)

  # labeling
  contour_plot <- ggplot2::ggplot(df, aes(x = x, y = y, z = z)) +
    ggplot2::geom_contour(breaks = z_contour)
  contour_data <- ggplot2::ggplot_build(contour_plot)$data[[1]]
  label_data <- aggregate(cbind(x,y) ~ level, data = contour_data, FUN = max)
  label_data$x <- theta
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
  res <- list(data=df, label=label_data, contour=z_contour, qval=z_threshold, alpha=alpha, plot=g)
  return(res)
}
