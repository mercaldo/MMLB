\name{mm}
\alias{mm}
\title{
Fit Marginalized Transition and/or Latent Variable Models
}
\description{
Fit a marginalzed transition and/or latent variable models (mTLV) as described by Schildcrout and Heagerty 2007.  
}
\usage{
mm(mean.formula, lv.formula = NULL, t.formula = NULL, id, data,
   inits = NULL, samp.probs = c(1, 1, 1), samp.probi=NULL, offset = 0, q = 10, 
   cond.like = FALSE, step.max = 1, step.tol = 1e-06, hess.eps = 1e-07, 
   adapt.quad = FALSE, verbose = FALSE,iter.lim=100)
}
\arguments{
  \item{mean.formula}{ Mean model formula in which a binary variable is regressed on covariates.}
  \item{lv.formula}{ Latent variable model formula (right hand side only)
}
  \item{t.formula}{Transition model formula (right hand side only)
}
  \item{id}{ a vector of cluster identifiers (it should be the same length nrow(data)).}
  \item{data}{ a required data frame }
  \item{inits}{an optional vector for initial values for marginal mean parameters and all dependence parameters. If NULL, initial values will be automatically generated.}
  \item{samp.probs}{a vector of 3 values that denote the sampling probability of non-responders, any-responders, and all-responders.}
  \item{samp.probi}{a vector of sampling probabilities - if using weighted estimating equations.}
  \item{offset}{an optional offset term.}
  \item{q}{a scalar to denote the number of quadrature points used for GH numerical integration. Only values of 3, 5, 10, 20 and 50 are applicable.}
  \item{cond.like}{indicator to denote if the conditional likelihood should be maximized.}
  \item{step.max}{a scalar.}
  \item{step.tol}{a scalar.}
  \item{hess.eps}{a scalar.}
  \item{adapt.quad}{an indicator if adaptive quadrature is to be used.  NOT CURRENTLY IMPLEMENTED.}
  \item{verbose}{an indicator if model output should be printed to the screen during maximization (or minimization of negative log-likelihood).}
  \item{iter.lim}{a scalar to denote the maximum iteration limit. Default value is 100.}
}

\value{
This function returns marginal mean (beta) and dependence parameters (alpha) along with the associated model and empirical covariance matricies.
}

\author{
Nathaniel Mercaldo
}

\seealso{
\code{\link{MMLongit}}
}

\examples{
\dontrun{
data(datrand)
fit <- mm(Y~time*binary, t.formula=~1, data=datrand, id=id, step.max=4, verbose=FALSE)}
}

