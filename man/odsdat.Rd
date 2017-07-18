\name{odsdat}
\alias{odsdat}
\docType{data}
\title{Simulated ODS data set}
\description{Simulated data set 2500 subjects according a mTLV model with marginal mean parameters (intercept, time, Xe (binary exposure), time)=(-1.75, 0.25, 0.25, 0.1), and dependence model parameters (sigma, gamma)=(1,1).  The marginal exposure prevalence of Xe was assumed to be 0.35, and all subjects had 10 observations. A two-phase ODS design was implemented using the phase 1 ODS design D1[25,50,25], and a phase 2 ODS design of D2[0,300,0]. The code to generate odsrand is provided in the Example section.}
\usage{data(odsdat)}
\format{
  A data frame with 3830 observations on the following 11 variables.
  \describe{
    \item{\code{id}}{a subject identifier}
    \item{\code{Y}}{a binary outcome}
    \item{\code{time}}{a time-varying covariate}
    \item{\code{Xe}}{a binary, time-invariant covariate (e.g., exposure of interest)}
    \item{\code{nobs}}{number of observations per subject; 10 in this example}
    \item{\code{sumY}}{subject-specific sum of Y; used to define sampling strata}
    \item{\code{ss}}{sampling strata; values 1-3}
    \item{\code{sp1}}{sampling probability for stratum 1}
    \item{\code{sp2}}{sampling probability for stratum 2}
    \item{\code{sp3}}{sampling probability for stratum 3}    
    \item{\code{sample}}{phase sampling indicator; 1=sampled in phase 1, 2=sampled in phase 2, NA=otherwise}   
  }
}

\examples{
\dontrun{
# Identify ids by sampling strata (ss) and compute basic stratum-specific summaries 
idsXss <- function(DAT) {
  DAT <- split(DAT, DAT$id)
  
  # Generate ids by each sampling stratum
  N      <- length(DAT)
  DAT_ss <- split(seq(N), factor(sapply(DAT, function(ZZ) ZZ[1,'ss']), 1:3))
  
  # Calculate frequencies in each sampling statrum and frequencies of sum(Y) and mean(Y)
  ss_sum  <- sapply(DAT_ss, length)
  sy_sum  <- table(sapply(DAT, function(ZZ) sum(ZZ$Y)))
  sy_mean <- table(sapply(DAT, function(ZZ) mean(ZZ$Y)))
  
  out <- DAT_ss
  attr(out, 'freq') <- list('ss'=ss_sum, 'sy'=sy_sum, 'sm'=sy_mean)
  out
}

# For a given design (e.g., D[n0, n1, n2]), use probability sampling to identify
# subjects to sample. get_sample() returns a list that contains sampling probabilities,
# and lists of sampled ids by ss and remaining ids by ss. When performing 2-phase ODS
# designs, we use the list of remaining ids (sampling without replacement), and thus
# our phase 2 sampling probabilities are conditioned on phase 1 data.

get_sample <- function(DESIGN, IDSXSS, RSEXACT, REPLACE=FALSE) {
  # DESIGN  = vector of expected sample counts
  # IDSXSS  = list of ids by sampling stratum (length=3)
  # RSEXACT = indicator to force exact random sampling
  # REPLACE = indicator if sampling with replacement should be performed

  tmp_design  <- as.list(DESIGN)
  blah_samp   <- vector('list',length(tmp_design))
  
  if( RSEXACT ) {
    blah <- do.call(rbind, mapply( function(AA,BB) list(cbind(AA,BB)), 
                                   AA=as.list( as.numeric(names(IDSXSS)) ),BB=IDSXSS))
    blah_tmp <-  data.frame( blah[ sort( sample(x=seq(nrow(blah)),round( sum(DESIGN),0) , replace=FALSE)) ,] )
    if(ncol(blah_tmp)==1) blah_tmp <- data.frame(t(blah_tmp))
    blah_samp     <- blah_tmp
    blah_samp[,1] <- factor(blah_samp[,1], levels=1:3)
    blah_samp     <- split(blah_samp$BB, blah_samp$AA)
  }
  
  tmp_samp    <- mapply( function(AA,BB,CC) {
    nn <- length(AA)
    if(nn < BB) BB <- nn
    sp <- ifelse(nn==0,0, BB/nn)
    sampled <- which( rbinom( nn, size=1, prob=sp )== 1 )
    if(length(sampled)>0 | !is.null(CC)) {
      if(is.null(CC)) {
        sampled <- AA[sampled]
      } else {
        sampled <- CC
      }
    }
    attr(sampled,'sp') <- sp
    list(sampled) }, AA=IDSXSS ,BB=tmp_design, CC=blah_samp)
  
  tmp_sp <- sapply( tmp_samp, function(ZZ) attr(ZZ,'sp') )
  
  tmp_remain <- mapply( function(AA,BB) {
    out <- AA
    if(length(BB)>0)  out <- AA[!(AA%in%BB)]
    out
  }, AA=IDSXSS, BB=tmp_samp)
  
  list('sp'=tmp_sp,'idsXss_samp'=tmp_samp,'idsXss_remain'=tmp_remain)
}

# Generate data
set.seed(1)
N      <- 2500
nclust <- rep(10,N)
id     <- rep(seq(N), nclust)
Xe     <- rep(rbinom(N,size=1,prob=.35), nclust) # binary exposure
time   <- unlist( lapply( as.list(nclust), function(ZZ) seq(ZZ)-1 ) )
data   <- data.frame(id, time, Xe)
data   <- data[order(data$id, data$time),]

# Generate response data
newdata <- GenBinaryY(mean.formula=~time*Xe, lv.formula=~1, t.formula=~1, 
                     beta=c(-1.75, .25, .25, .1), sigma=1, gamma=1, id=id, data=data, Yname = "Y") 

# Generate sampling strata and sample indicator
newdata <- do.call(rbind, lapply( split(newdata, newdata$id), function(ZZ) { 
  ZZ$nobs   <- nrow(ZZ)
  ZZ$sumY   <- sum(ZZ$Y)
  ZZ$ss     <- ifelse(ZZ$sumY==ZZ$nobs, 3, ifelse(ZZ$sumY==0, 1, 2))
  ZZ$sample <- NA
  ZZ}))

# Create list of ids by ss. 
tmp_ss <- idsXss(DAT=newdata) 
(n_ss  <- sapply(tmp_ss, length))

# Create sampling probabilities, and identify ODS sample. Suppose it is of 
# interset to perform a single phase (or a phase one) ODS design of the form
# D[25,50,25].

# Obtain phase 1 ODS sample
s1d      <- apply( cbind( c(25,50,25), attr(tmp_ss,'freq')$ss ), 1, min )
s1info   <- get_sample(DESIGN=s1d, IDSXSS=tmp_ss, RSEXACT=FALSE)
s1info$sp

s1id     <- sort( unique(unlist( s1info$idsXss_samp )) )
s1idX    <- which(newdata$id%in%s1id)
newdata$sample[s1idX] <- 1

# Obtain phase 2 ODS sample
tmp_ss_2      <- s1info$idsXss_remain
tmp_ss_2n     <- sapply(tmp_ss_2,length)
prob_sampled  <- matrix( s1info$sp,nrow=1 )
colnames(prob_sampled) <- c('sp1','sp2','sp3')

s2d      <- apply( cbind( c(0,300,0), sapply(tmp_ss_2,length)), 1, min )
s2info   <- get_sample(DESIGN=s2d, IDSXSS=tmp_ss_2, RSEXACT=FALSE) 
s2info$sp
s2idX    <- which(newdata$id %in% unlist(s2info$idsXss_samp))
newdata$sample[s2idX] <- 2 

# Estimate phase 2 sampling probabilities
prob_sampled        <- data.frame( rbind(s1info$sp,s2info$sp) )
prob_not_sampled    <- apply( 1-prob_sampled, 2, cumprod)                  # cumulative prob of not being sampled
prob_not_sampled    <- rbind(1, prob_not_sampled[-nrow(prob_not_sampled),])   
prob_sampled        <- prob_sampled * prob_not_sampled 
prob_sampled$sample <- c(1,2)
colnames(prob_sampled)[1:3] <- c('sp1','sp2','sp3')
prob_sampled

# Update newdata with sampling probabilities
newdata2 <- merge(newdata, prob_sampled, by='sample',all.x=TRUE)
odsdat  <- newdata2[ order(newdata2$id, newdata2$time), ]
odsdat  <- odsrand[!is.na(odsdat $sample), ]
}
}


