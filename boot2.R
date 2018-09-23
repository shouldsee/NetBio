library(mixtools)
library(boot)
#### Bootstrap is sensitive to the mixture nature of a distribution

# Fc <- function(x,i){colSums(x[i,])}
Fc <- function(x,i){apply(x[i,],2,median)}
N = 1200
R = 2000
{
  mu1= c(0,0)
  mu2= c(6,6)
  mu0= (mu1+mu2)/2
  Cmu <- (outer(mu1-mu0,mu1-mu0) + outer(mu2-mu0,mu2-mu0))/2
  C0 = matrix(c(4,2,2,4),ncol=2)
  C = Cmu + C0
  r1 <- rmvnorm(N/2,mu1,C0)
  # plot(r[,1],r[,2])
  r2 <- rmvnorm(N/2,mu2,C0)
  rmix <- rbind(r1,r2)
  
  
  plot(r1,xlim=c(-6,12),ylim=c(-6,12))
  points(r2,col=2)
  hist(rmix[,1])
  grid()
  # C <- cov(rmix)
}

# cor(rbind(r1,r2))
{
  r3 <- rmvnorm(N,c(2,2),C)
  plot(r3)
  hist(r3[,1])
  cor(r3)
}


{
  dat <- rmix
  res <- boot(dat,Fc,R=R)
  dev <- t(res$t)- as.vector(res$t0)
  SD <- apply(dev,1,sd)
  # SD <- matrix(SD,dim(res$t0)) #* sqrt(res$R)
  Z  =res$t0 / SD
  P = pnorm(abs(Z),lower.tail = F)
  # diag(Z) <-0
  print(cov(dat))
  print(cor(dat))
  print(SD)
  print(Z)
}


{
  dat <- r3
  res <- boot(dat,Fc,R=R)
  dev <- t(res$t)- as.vector(res$t0)
  SD <- apply(dev,1,sd)
  # SD <- matrix(SD,dim(res$t0))
  # SD <- matrix(SD,dim(res$t0))# * sqrt(res$R-1)
  Z  =res$t0 / SD
  P = pnorm(abs(Z),lower.tail = F)
  # diag(Z) <-0
  print(cov(dat))
  print(cor(dat))
  print(SD)
  print(Z)
}
# ?boot.se()
stop()
{
  R = 100
  # dat <- expr.dat
  dat <- with(env.sub,expr.dat)
  # res <- boot(dat,Fc,R=R)
  dat <- as.matrix(dat)
  res <- boot(dat,Gc,R=R,parallel = 'snow',cl=clu)
  dev <- t(res$t)- as.vector(res$t0)
  C <- res$t0
  SD <- apply(dev,1,sd)
  SD <- matrix(SD,dim(res$t0)) #* sqrt(res$R)
  Z  =res$t0 / SD
  P = pnorm(abs(Z),lower.tail = F)
  # diag(Z) <-0
  print(cov(dat))
  print(cor(dat))
  print(SD)
  print(Z)
  image(SD)
}
pipeline(abs(C))
pipeline(abs(Z))
# image(-)

with(env.sub,pipeline(-abs(C)) )
with(env.sub,pipeline(-abs(Z)) )
expand_dims
Gc <- function(x,i){
  FUN = '*'
  m<-x[i,]; M <- bsxfun(t(m),expand_dims(m,1),FUN);
  mm<-apply(m,2,median);
  # mvar <- apply(m^2,2,mean)
  (apply(M,c(1,3),median) - outer(mm,mm,FUN))
  # (apply(M,c(1,3),median) - outer(mm,mm,FUN))/sqrt(outer(mvar,mvar))
}
image(Gc(dat,1:10))

library(GGally)
# ggpairs(as.datadat[,1:5])

idx=25:30
idx=21:25
idx=16:20
idx=11:15
idx=5:10
idx = 1:5
x <- apply(dat, 2, entropise)
od <- order(x,decreasing = T)
# ggpairs(as.data.frame(dat[,od[idx]]))
# with(env.sub,{


# plot(expr)
# apply(dat[,od[idx]],2,hist)
plot(apply(dat,2,sd),SD)
plot(apply(dat,2,entropise),SD)
grid()
# source('minfo.R')
stop()
# boot.comp(expr.dat,max.comp = 10,N=50)
{
  x <- expr.dat[,1:2]
  res = mvnormalmixEM(x,k=2)
  # res$loglik
  res$all.loglik
}
plot(x)
# res$loglik
# library(ggplot2)
# plotmatrix(x)
pairs(expr.dat[,1:5])
with(env.sub,
     pairs(expr.dat[,15:20])
)
# hist(x)
plot(res$posterior)

rdata(RTdata)
set.seed(100)
# x <- as.matrix(RTdata[, 1:3])
# y <- makemultdata(x, cuts = quantile(x, (1:9)/10))$y
# x
