print.summary.MMLongit <-
function(x,...) {
  cat('\nClass:\n',x$class,'\n',sep='')
  cat("\nCall:\n", paste(x$call, sep = "\n", collapse = "\n"),"\n", sep = "")
  cat('\nInformation Criterion:\n')
  print.default(format(x$info, ...), print.gap = 2L, quote = FALSE)
  cat("\nMarginal Mean Parameters:\n")
  printCoefmat(x$mean.table,signif.stars = FALSE)
  cat('\n')
  cat("Association Parameters:\n")
  printCoefmat(x$assoc.table,signif.stars = FALSE)
  cat('\n')
  cat('Number of clusters:            ',x$control[5],'\n')
  cat('Maximum cluster size:          ',x$control[6],'\n') 
  cat('Convergence status (nlm code): ',x$control[3],'\n')
  cat('Number of iterations:          ',x$control[4])
  cat('\n')
}
