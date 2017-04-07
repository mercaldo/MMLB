\name{datrand}
\alias{datrand}
\docType{data}
\title{Simulated data set}
\description{ A simulated data set. Data were created using fixed marginal mean parameters (beta0, beta1, beta2, beta3) = (-1.85, -0.15, 1.00, 0.15) and association parameters (gamma, sigma) = (1.5, 0.0). These data were created assuming an autocorrelation dependence structure.}
\usage{data(datrand)}
\format{
  A data frame with 24999 observations on the following 4 variables.
  \describe{
    \item{\code{id}}{a patient identifier}
    \item{\code{Y}}{a binary outcome}
    \item{\code{time}}{a time-varying covariate}
    \item{\code{binary}}{a time-invariant covariate}
  }
}


