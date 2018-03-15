{
  d.e2d <- pdist(expr.dat,partial(entropy2d,nbin=c(50,50)))
  d.hx <- cbind(diag(d.e2d))
  # d.hx <- cbind(apply(expr.dat,2,entropise))
  d.H <- bsxfun(t(d.hx),d.hx,'+')
  d.Hmax <- bsxfun(t(d.hx),d.hx,pmax)
  d.Hmin <- bsxfun(t(d.hx),d.hx,pmin)
  d.CH <-bsxfun(d.e2d,d.hx,'-') 
  mean(d.CH>0)%>%print
  d.mi <- bsxfun(t(d.hx),d.CH,'-')
  # d.mi.norm <- d.mi/d.Hmax
  d.mi.norm <- d.mi/d.Hmin
  image(d.CH)
  image(d.mi.norm)
}

# mat.sqrt <- function(m,inv = MASS::ginv){
#   ip <- inv(m)
#   # sqrtm
#   pa <- inv(ip%*%ip)
# }
{
  par(mfrow=c(1,3))
  image(log(expm(array2array(adj.true,as.numeric))))
  image(log(expm(exp(-distances(g.true)))))
  dmat <- array2array(adj.true,as.numeric)
  dmat[!adj.true] <- Inf
  C = 10
  image(dmat.t <- log(expm(exp(-C*dmat))))
}

idx <- which(adj.true,arr.ind = T)
g.dist <- distances(g.true)
# plot(expr.dat[,idx[11,]])
# load.assignment.data()

