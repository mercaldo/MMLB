get.GH <- function(q, scale_abscissa = sqrt(2), scale_weight=1/sqrt(pi)) {
  rule = gaussHermiteData(q)
  if(scale_abscissa!=1) rule$x = rule$x*scale_abscissa
  if(scale_weight!=1)   rule$w = rule$w*scale_weight
  list(z=rule$x,w=rule$w)
}
