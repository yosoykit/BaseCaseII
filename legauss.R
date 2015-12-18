legauss <- function(xl,xu,ngs) {
  tmp <- .Fortran("legaus",xl=as.double(xl),xu=as.double(xu),ngs=as.integer(ngs),xg=numeric(ngs),wg=numeric(ngs))
  return(list(ngs=ngs,mesh=tmp$xg,weights=tmp$wg))
}

