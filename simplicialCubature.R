ddiri <- function(x,eta=rep(1,length(x)),log=T,check=F){
  if (check){
    if(sum(x)>1){return(0)}
    
  }
  if (log){
    cst = sum(lgamma(eta))-lgamma(sum(eta))
    lp = -cst + sum((eta-1)*log(x))
    if(is.na(lp)){
      lp=-Inf
    }
    lp
  }else{
    cst = prod(gamma(eta))/gamma(sum(eta))
    p = 1/cst*prod(x^(eta-1))
    # exp(lp)
  }
}
log(10^59)

{
  # eta = c(2,2,2,2)
  eta = c(1,1,1,1)
  library(SimplicialCubature)
  main<-function(dat.dec){
    print('')
  lst = c()
  ds = 2^seq(-2,2,length.out = 10)
  for (d in ds){
    
    S <- UnitSimplexV(4)
    # measure <-  adaptIntegrateSimplex( function(x)1, S)$integral
    # sf = SimplexSurfaceArea(S)
    # d = 20
    # wt <-function(x){dnorm(x[1]*x[3]-x[2]*x[4],sd=d)}
    # wt <- function(x) dnorm((x[1]*x[3]-x[2]*x[4])/(x[1]*x[3]+x[2]*x[4]),sd=d)
    # wt <- function(x) ddiri(x,log=F,eta,check=F) * dnorm((x[1]*x[3]-x[2]*x[4]),sd=d)
    # wt <- function(x) ddiri(x,log=F,eta,check=F) * dnorm(MI(x),sd=d)
    wt <- function(x) ddiri(x,log=F,eta,check=F) * dnorm((x[1]*x[3]-x[2]*x[4])/(x[1]*x[3]+x[2]*x[4]),sd=d)
    # wt <- function(x) ddiri(x,log=F,eta,check=F) * dnorm((x[1]*x[3]-x[2]*x[4]),sd=d)
    cst <-  (res<-adaptIntegrateSimplex( function(x) wt(x),S,maxEvals = 10000))$integral
    # print(res$estAbsError)
    # cst2<- adaptIntegrateSimplex( function(x) prod(x[dat.dec+1]),S)$integral
    
    # print(bquote(Normalising_constant==.(cst)))
    
    # print(cst2)
    # S = SimplicialCubature::CanonicalSimplex(5)
    # sf = SimplexVolume(S)
    # eta = c(2,2,2,2)
    # eta = rep(10,nrow(S))
    
    # eta12=sum(eta[1:2])
    # eta34=sum(eta[3:4])
    # eta13=sum(eta[c(1,3)])
    # eta24=sum(eta[c(2,4)])
    f1 <- function( x ) { 
      # r <- x[1]+x[2]+x[3]
      # r <- dbeta(x[1],1,2)
      # print(sum(x))
      x12 <- x[1]+x[2]
      x13 <- x[1]+x[3]
      # r = dbeta(x12,eta12,eta34)
      # r = r*dbeta(x13,eta13,eta24)
      
      # r = dbeta(x12,eta12,eta12)
      # r = r*dbeta(x13,eta12,eta12)
      
      # print(x)
      # r = dbeta(x[1]+x[2],)
      # r = ddiri(x,log=F,eta,check=F) 
      r = 1
      # r= r * dnorm((x[1]*x[3]-x[2]*x[4]),sd=d)
      p <- exp(mean(log(x)[dat.dec+1]))
      r= r * wt(x)*p
      # r = r * 
      
      # /gamma(nrow(S))
      
      # r/sf
      
      # r/sqrt(nrow(S))
      # dbeta
    }
    
    res2 = adaptIntegrateSimplex( f1, S,maxEvals = 20000)
    r = res2$integral/sqrt(nrow(S))/cst
    print(bquote(Likelihood==.(r)~val1==.(res$integral)~err1==.(res$estAbsError)~err2==.(res2$estAbsError)))
    lst = c(lst,r)
    # plot(lst,type='b')
  }
  lst
  }
  
# plot(lst,type='b')
}

lst <- main(dat.dec.list[1,])
plot(ds,lst,type='b',log='x')
lst

library(Rutil)
dat.logL <- apply(dat.dec.list[1:3,],1,main)
{
  par(mfrow=c(1,3))
  for (i in 1:ncol(dat.logL)){
    plot(ds,dat.logL[,i],type='b',log='x')
  }
  # matplot(log(dat.logL),type='b')
}
# .Machine$double.eps
{
  dat0=dat2
  dat.dec.list = cbind()
  dat <- dat0[,1:2]
  dat.dec = bin2dec(as.matrix(dat)=='1')
  dat.dec.list <- rbind(dat.dec.list,dat.dec)
  dat <- dat0[,2:3]
  dat.dec = bin2dec(as.matrix(dat)=='1')
  dat.dec.list <- rbind(dat.dec.list,dat.dec)
  dat <- dat0[,c(1,3)]
  dat.dec = bin2dec(as.matrix(dat)=='1')
  dat.dec.list <- rbind(dat.dec.list,dat.dec)
  dat.dec.list
}

# dat0
hist(dat.dec.list[1,])
bn <-bnlearn::hc(dat1,score = 'bde')
plot(bn)
# dat.dec.list

dnorm(x[1]*x[3]-x[2]*x[4],sd=d)
{
  library(MCMCpack,lib.loc = '/data/subliminal/RLib/temp/')
  N = 3
  # eps = 1E-3
  eps = 1E-1
  alpha  = rep(10,N)
  eta = 6
  lo = rep(eps,N);hi = rep(1-eps,N)
  f = function(x){
    # print(x)
    # dbeta(x[1],2,2)*dbeta(x[2],2,2)
    # ddirichlet(c(x,min(0,1-sum(x)) ),alpha)
    # r = ddiri(c(x,min(0,1-sum(x))) ,alpha,log=F)
    x = order(x)/max(x)
    x = c(x[1],diff(x))
    # r = ddirichlet(x,alpha)
    # x12 <- x[1]+x[2]
    # x13 <- x[1]+x[3]
    # r = dbeta(x12,eta,eta)
    # r = r*dbeta(x13,eta,eta)

    
    r= ddiri(x,alpha,log=F,check=F)
    
    # r = ddiri(c(x,min(0,1-sum(x))) ,alpha,log=F)
    # x = e_to_p(x)
    # r = ddiri(c(x),alpha[1:N])
    # r = ddiri(c(x,1-sum(x)),alpha)
    # hyper2::
  }
  # print(f(c(3:1)*5))
  cubature::hcubature(f,lo,hi)
}

# dhyperdirichlet_e( c(.5,.8),dirichlet(1:3))
?hyperdirichlet::p_to_e
library(hyperdirichlet)
warnings()
detach('package:MCMCpack',unload=T)
library(Rutil)
install.packages.lazy('MCMCpack')
# cubature::
# dbeta()
S = SimplicialCubature::CanonicalSimplex(3)

SimplexVolume(S)

install.packages('aylmer')
install.packages('abind')
install.packages('https://cran.r-project.org/src/contrib/Archive/hyperdirichlet/hyperdirichlet_1.5-1.tar.gz',repos=NULL,dependencies = T)
res$integral
sqrt(2)
sqrt(3/4)
# sqrt()

# 0.8660254^2
# (3^(1/3))/3
# lp.diri
MCMCpack::
  install.packages.lazy('MCMCpack')

install.packages("hyper2",dependencies = T)
a = ddiri(c(0,.5,0),c(3,2,2),log=F)
a
library(MCMCpack,lib.loc = '/data/subliminal/RLib/temp/')
# MCMCpack::diri  
# ddirichlet(c(1,),c(2,2))
# ddiri
# ddiri
integrate(function(x)dbeta(x,2,2),0,1)
xs = (0:50)/50
plot(dbeta(xs,2,2))
plot(sapply(xs,function(x)ddiri(c(x,1-x),log=F,c(2,2) )))
# res
# rbeta()
install.packages("MCMCpack",'/data/subliminal/RLib/temp/')
res$integral
