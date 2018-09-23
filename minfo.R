
#' Estimate entropy by KDE-fitted distribution 
#' @export
entropise <- function(x,eps=1E-3){
  ds <- density(x)
  dx <- diff(ds$x[1:2])
  lp <- ds$y*dx
  ent <- sum(lp*-log(lp+eps))
  # eps =
  # ent = -sum(ds$y*log(ds$y+eps))*dx
  # ent <- do.call(integrate,c(f=approxfun(ds$x,ds$y*log(ds$y+eps)),as.list(range(ds$x))))
}


#' Estimate total entropy of 2 r.v. by 2d-KDE
#' @export
entropy2d<-function(x,y,nbin=c(25,25),eps=1E-3){
  ds <- MASS::kde2d(x,y,n = nbin)
  # eps = 0
  dx <- diff(ds$x[1:2])
  dy <- diff(ds$y[1:2])
  lp <- ds$z*dx*dy
  -sum(lp*log(lp+eps))
}

wrong.entropy2d<-function(x,y,nbin=c(25,25),eps=1E-3){
  ds <- MASS::kde2d(x,y,n = nbin)
  # eps = 0
  dx <- diff(ds$x[1:2])
  dy <- diff(ds$y[1:2])
  ent <- -sum(ds$z*(log(ds$z+eps)))*dx*dy
}


#' Make pairwise MI matrix from a dataset
#' @export
make.mi <- function(expr.dat){
  d.e2d <- pdist(expr.dat,partial(entropy2d,nbin=c(30,30)))
  # d.e2d <- proxy::dist(expr.dat
  #                      ,method=partial(entropy2d,nbin=c(30,30))
  #                      ,diag = T,by_rows = F) %>%as.matrix
  expr.h <- apply(expr.dat,2,entropise)
  d.h <- outer(expr.h,expr.h,'+')
  d.mi = d.h - d.e2d 
}

#' Make pairwise distance matrix
#' @export
pdist<-function(X,distF){
  D <- proxy::dist(X,method =distF,by_rows = F,diag = T)%>%as.matrix
  diag(D)<-apply(X,2,function(x)distF(x,x))
  D
}