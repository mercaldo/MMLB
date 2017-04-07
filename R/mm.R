mm <- function(mean.formula, lv.formula = NULL, t.formula = NULL, id, data, inits = NULL,
         samp.probs = c(1, 1, 1), samp.probsi = NULL, offset = NULL, q = 10, 
         cond.like = FALSE, step.max = 1, step.tol = 1e-06, hess.eps = 1e-07, 
         adapt.quad = FALSE, verbose = FALSE, iter.lim=100) {

  if(is.null(id)) {stop('Provide id vector or variable name!')}
#  if(as.character(as.name(id))%in%colnames(data)) id = data[,as.character(as.name(id))]
  
#  id.f = model.frame(id, data)
#  id   = model.extract(id.f, id)
  
  mean.f = model.frame(mean.formula, data)
  mean.t = attr(mean.f, "terms")
  y  = model.response(mean.f,'numeric') # check if its binary?
  uy = unique(y)
  if(length(uy)!=2 | !setequal(uy,c(0,1))) stop('Problem: Check outcome! Binary? Response variation?')
  x  = model.matrix(mean.formula,mean.f)
  
  x.t = x.lv = matrix(0, ncol=1, nrow=length(y))
  if(!is.null(t.formula))   x.t  = model.matrix(t.formula,model.frame(t.formula, data)) 
  if(!is.null(lv.formula))  x.lv = model.matrix(lv.formula, model.frame(lv.formula, data)) 
  
  samp.probs = matrix(samp.probs,nrow=length(y),ncol=3,byrow=TRUE) #(never, any, always)
  
#  if(cond.like | length(unique(samp.probs))==1) {
  if(is.null(samp.probsi)) {
    samp.probsi = matrix(1,nrow=length(y),ncol=1) 
  } #else {
    #samp.probsi = matrix(unlist(lapply(split(y,id),function(Z) {rep(samp.probs[1]*(sum(Z)==0) + samp.probs[3]*(sum(Z)==length(Z)) +
    #                                                                 samp.probs[2]*(sum(Z)!=0 & sum(Z)!=length(Z) ) ,length(Z))})),ncol=1)
  #}
  
  if(is.null(inits)) {
    inits = c(glm(mean.formula,family='binomial',data=data)$coef, rep(1, ncol(x.t) + ncol(x.lv)))
  }
  
  if(is.null(offset)) {
    offset <- rep(0, length(y))
  }
  
  mm.fit = MMLongit(params=inits, id=id, X=x, Y=y, Xgam=x.t, Xsig=x.lv, Q=q, condlike=cond.like,
                    sampprobs=samp.probs, sampprobi=samp.probsi, offset=offset, 
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
  out$emp.cov = mm.fit$empiricalcov
  names(out$beta) = nms$beta
  names(out$alpha) = nms$alpha
  colnames(out$mod.cov) = rownames(out$mod.cov) = colnames(out$emp.cov) = rownames(out$emp.cov) = unlist(nms)
  out$control = with(mm.fit, c(condlike, AdaptiveQuad, code, niter, length(table(id)), max(table(id))))
  
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
  class(out) = 'MMLongit'
  out
}
