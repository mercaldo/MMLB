mm <- function(mean.formula, lv.formula = NULL, t.formula = NULL, id, data, inits = NULL,
         samp.probs = c(1, 1, 1), samp.probi = NULL, offset = NULL, q = 30, 
         cond.like = FALSE, step.max = 1, step.tol = 1e-06, hess.eps = 1e-07, 
         adapt.quad = FALSE, verbose = FALSE, iter.lim=100, return_args=FALSE) {
  
  if(is.null(lv.formula) & is.null(t.formula)) {stop('Specify association model (both lv.formula and t.formula arguments cannot be NULL.')}
  if(!is.data.frame(data)) {
    data <- as.data.frame(data)
    warning('data converted to data.frame.')
  }
  
  terms = unique( c(all.vars(mean.formula), all.vars(lv.formula), all.vars(t.formula), 
                    as.character(substitute(id))) )
  data  = data[,terms]
  if(any(is.na(data))) data = na.omit(data)
  id    = data$id    = data[ , as.character(substitute(id)) ]
  
  if(q<=1) q <- 2
  if(is.null(lv.formula)) q = 1
  
  mean.f = model.frame(mean.formula, data)
  mean.t = attr(mean.f, "terms")
  y  = model.response(mean.f,'numeric') 
  uy = unique(y)
  x  = model.matrix(mean.formula,mean.f)
  
  x.t = x.lv = matrix(0, ncol=1, nrow=length(y))
  if(!is.null(t.formula))   x.t  = model.matrix(t.formula,model.frame(t.formula, data)) 
  if(!is.null(lv.formula))  x.lv = model.matrix(lv.formula, model.frame(lv.formula, data)) 
  
  if (!is.matrix(samp.probs)) {
    samp.probs = matrix(samp.probs, nrow = length(y), ncol = 3, byrow = TRUE)
  }
  
  useROBCOV <- TRUE
  if(is.null(samp.probi)) {
    samp.probi = matrix(1,nrow=length(y),ncol=1) 
    useROBCOV <- FALSE
  }
  
  if (!is.null(inits)) {
    inits <- unlist(inits)
  }
  
  if(is.null(inits)) {
    inits = c(glm(mean.formula,family='binomial',data=data)$coef, rep(1, ncol(x.t) + ncol(x.lv)))
    if(any(is.na(inits))) {
      omit_dup_col = which(is.na(inits))
      x            = x[,-c(omit_dup_col)]
      inits        = inits[-c(omit_dup_col)]
    }
  }
  
  if(is.null(offset)) {
    offset <- rep(0, length(y))
  }
  
  mm.fit = MMLongit(params=inits, id=id, X=x, Y=y, Xgam=x.t, Xsig=x.lv, Q=q, condlike=cond.like,
                    sampprobs=samp.probs, sampprobi=samp.probi, offset=offset, 
                    stepmax=step.max, steptol=step.tol, hess.eps=hess.eps, 
                    AdaptiveQuad=adapt.quad, verbose=verbose,iterlim=iter.lim)
  
  nms = list()
  nms$beta = colnames(x)
  nms$alpha = c( if(!is.null(t.formula)){paste('gamma',colnames(x.t),sep=':')}, 
                 if(!is.null(lv.formula)){paste('log(sigma)',colnames(x.lv),sep=':')})
  
  out = NULL
  out$call    = match.call() 
  out$logLik = mm.fit$logL
  out$beta  = mm.fit$beta
  out$alpha = mm.fit$alpha
  out$mod.cov = mm.fit$modelcov
  out$rob.cov = mm.fit$empiricalcov
  names(out$beta) = nms$beta
  names(out$alpha) = nms$alpha
  colnames(out$mod.cov) = rownames(out$mod.cov) = colnames(out$rob.cov) = rownames(out$rob.cov) = unlist(nms)
  out$control = with(mm.fit, c(condlike, AdaptiveQuad, code, niter, length(table(id)), max(table(id)),
                               useROBCOV))
  names(out$control) <- c('cond.like','adaptive.quad','convergence_code','n_iter','n_subj','max_n_visit','useRobCov')

  aic = function(l=mm.fit$logL,k=nrow(mm.fit$modelcov)) 2*k-2*l
  bic = function(l=mm.fit$logL,k=nrow(mm.fit$modelcov),n=length(table(id))) -2*l + k *log(n) 
  deviance = function(l=mm.fit$logL) -2*l
  out$info_stats = c(aic(),bic(),mm.fit$logL,deviance())
  names(out$info_stats) = c('AIC','BIC','logLik','Deviance')
  out$m.formula = mean.formula
  out$t.formula = t.formula
  out$lv.formula = lv.formula
  out$LogLikeSubj = mm.fit$LogLikeSubj
  out$ObsInfoSubj = mm.fit$ObsInfoSubj
  out$ACSubj      = mm.fit$ACSubj
  out$LLSC_args   = mm.fit$LLSC_args
  
  if(return_args) attr(out,'args') <- list('mean_formula'=mean.formula, 't_formula'=t.formula, 'lv_formula'=lv.formula,
                                           'id'=id, 'data'=data, 'samp.probs'samp.probs=,'samp.probi'=samp.probi,'offset'=offset)
  class(out) = 'MMLongit'
  out
}
