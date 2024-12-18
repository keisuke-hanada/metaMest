% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pdf-theta-Mest.R
\name{ptheta_with_M_estimator}
\alias{ptheta_with_M_estimator}
\title{Probability distribution function of overall effect by M-estimator based approach}
\usage{
ptheta_with_M_estimator(
  qval,
  yi,
  vi,
  theta = 0,
  tau2,
  tau2b = c(0.001, 1:300/10),
  n.tau2 = 10000,
  method = "DL"
)
}
\arguments{
\item{qval}{value of overall effect}

\item{yi}{vector of length K with the observed effect sizes or outcomes.}

\item{vi}{vector of length K with the corresponding sample variances.}

\item{theta}{value of true overall effect that define 0 in default.}

\item{tau2}{value of estimated between-study variance.}

\item{tau2b}{vector of candidate values of M-estimator based parameter. If the default is insufficient, the search area can be increased by increasing the upper limit or narrowing the interval.}

\item{n.tau2}{value of number of sampling for estimates of between-study variance.}

\item{method}{character string to specify the estimating method of between-study variance. By default (method="DL"), DerSimonian-Laird estimator is obtained. When method="SJ", Sidik-Jonkman estimator is obtained. When method="HKSJ-DL", DerSimonian-Laird method with HKSJ approach is obtained. When method="HKSJ-SJ", Sidik-Jonkman method with HKSJ approach is obtained.}
}
\value{
Probability of qval and M-estimator based estimates of location parameter.
}
\description{
Probability distribution function of overall effect by M-estimator based approach
}
