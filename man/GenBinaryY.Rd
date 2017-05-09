\name{GenBinaryY}
\alias{GenBinaryY}
\title{
Generate binary response data.
}
\description{
Generate binary response data for a given marginal mean and association model. Association models include those with a random intercept and/or a first-order time lag. 
}
\usage{
GenBinaryY(mean.formula, lv.formula = NULL, t.formula = NULL, 
    beta = NULL, sigma = NULL, gamma = NULL, id, data, q = 10, 
    Yname = "Y")
}
\arguments{
  \item{mean.formula}{ Right hand side of mean model formula.}
  \item{lv.formula}{ Latent variable model formula (right hand side only) }
  \item{t.formula}{Transition model formula (right hand side only) }
  \item{beta}{ a vector of values for mean.formula.}
  \item{sigma}{ a vector of values for the latent variable portion of the association model (else NULL).}
  \item{gamma}{ a vector of values for the transition porition of the association model (else NULL).}
  \item{id}{ a vector of cluster identifiers (it should be the same length nrow(data)).}
  \item{data}{ a required data frame }
  \item{q}{a scalar to denote the number of quadrature points used for GH numerical integration.}
  \item{Yname}{a character string of the name of new binary variable; default=Y}
}

\value{
This function returns the data object augmented with a binary response vector named Yname. 
}

\author{ Jonathan Schildcrout and Nathaniel Mercaldo }

\examples{
\dontrun{
set.seed(1)
N      = 1500
nclust = sample( seq(10,10), N, replace=TRUE)
id     = rep(seq(N), nclust)
Xe     = rep(rbinom(N,size=1,prob=.5), nclust) # binary exposure
time   = unlist( sapply( as.list(nclust), function(ZZ) seq(ZZ)-1 ) )
data   = data.frame(id, time, Xe)
data   = data[order(data$id, data$time),]

newdata = GenBinaryY(mean.formula=~time*Xe, lv.formula=~1, t.formula=~1, 
          beta=c(-2.5, .25, .25, .1), sigma=1, gamma=1, id=id, data=data, q=20, 
          Yname = "binY") 

mod_mtlv = mm(binY~time*Xe,lv.formula=~1, t.formula=~1, data=newdata,id=id, q=20)
summary(mod_mtlv)
}
}