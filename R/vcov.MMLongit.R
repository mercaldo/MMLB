vcov.MMLongit <- function (object, robust=FALSE) {
  n_beta     = length(object$beta)
  n_alpha    = length(object$alpha)
  beta_ix    = seq(n_beta)
  alpha_ix   = seq(n_beta+1, n_beta+n_alpha)
  object_vc  = object$mod.cov
  if(robust) object_vc  = object$emp.cov
  list('beta'=object_vc[beta_ix,beta_ix], 'alpha'=object_vc[alpha_ix,alpha_ix])
}