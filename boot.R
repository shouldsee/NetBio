library(mixtools)
library(boot)
#### Bootstrap is sensitive to the mixture nature of a distribution

N = 1200
{
  mu1= c(0,0)
  mu2= c(6,6)
  mu0= (mu1+mu2)/2
  Cmu <- (outer(mu1-mu0,mu1-mu0) + outer(mu2-mu0,mu2-mu0))/2
  C0 = matrix(c(4,2,2,4),ncol=2)
  C = Cmu + C0
  r1 <- rmvnorm(N/2,mu1,C0)
  # plot(r[,1],r[,2])
  plot(r1,xlim=c(-6,12),ylim=c(-6,12))
  r2 <- rmvnorm(N/2,mu2,C0)
  points(r2,col=2)
  grid()
  rmix <- rbind(r1,r2)
  # C <- cov(rmix)
}

# cor(rbind(r1,r2))
{
  r3 <- rmvnorm(N,c(2,2),C)
  plot(r3)
  cor(r3)
}


{
  dat <- rmix
  res <- boot(dat,function(x,i)cor(x[i,]),R=8000)
  dev <- t(res$t)- as.vector(res$t0)
  SD <- apply(dev,1,sd)
  SD <- matrix(SD,dim(res$t0)) #* sqrt(res$R)
  Z  =res$t0 / SD
  P = pnorm(abs(Z),lower.tail = F)
  diag(Z) <-0
  print(cov(dat))
  print(cor(dat))
  print(SD)
  print(Z)
}


{
  dat <- r3
  res <- boot(dat,function(x,i)cor(x[i,]),R=8000)
  dev <- t(res$t)- as.vector(res$t0)
  SD <- apply(dev,1,sd)
  SD <- matrix(SD,dim(res$t0))
  SD <- matrix(SD,dim(res$t0))# * sqrt(res$R-1)
  Z  =res$t0 / SD
  P = pnorm(abs(Z),lower.tail = F)
  diag(Z) <-0
  print(cov(dat))
  print(cor(dat))
  print(SD)
  print(Z)
}
# ?boot.se()
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
