source('header.R')
library(boot)
load.assignment.data()
load('e-coli-subnet.rdf5')
{
  # clu = parallel::makeCluster(10)
  env.sub <- subnet(env.data,ngene = 30)
  # with(env.sub,
  with(.GlobalEnv,
       {
    par(mfrow=c(4,3),
        mar=c(2.5,1,2.5,1))
    expr.df <- data.frame(expr.dat)
    {
      
      res <- routine.bnlearn.bootstrap(expr.dat,cluster=clu)
      pipeline(res);
      title('Bayes',line = -1)
    }
    
    {
      # res <- cor(expr.dat,method='pearson')
      res <- cor(expr.dat,method='spearman')
      # res <- make.mi(expr.dat)
      res <- abs(res)
      pipeline(res)
      title('spearman',line = -1)
      # title('PPMCC',line = -1)
    }
  
    # {
    #   res <- boot.strength(expr.dat.bin,algorithm = 'aracne',cluster = clu)
    #   pipeline(res)
    #   title('Aracne-Boot',line = -1)
    # }
    
    {
      out <- boot(as.matrix(expr.dat),function(x,i){cor(x[i,])},R=800)
      dev <- t(out$t)- as.vector(out$t0)
      SD <- apply(dev,1,sd)
      SD <- matrix(SD,dim(out$t0)) #* sqrt(out$R)
      C = out$t0
      Z  =out$t0 / SD
      P = pnorm(abs(Z),lower.tail = F)
      # dat <- -log(P)
      # dat <- Z
      # diag(Z) <-0
      # res <- boot.strength(as.data.frame(dat),algorithm = 'chow.liu',cluster = clu)
      # res <- minet::minet(dat)
      # res <- minet::minet(abs(Z))
      # res <- cor(Z,method='spearman')
      res = abs(Z)
      pipeline(-log(P))
      pipeline(abs(C)*SD)
      # pipeline(res)
      
      # SD = SD* sqrt(out$R)
      # Z  =out$t0 / SD
      # diag(Z) <-0
      # pipeline(abs(Z))
    }
    # expr.dat.bin <- bnlearn::discretize(expr.dat)
    # pipeline(minet::aracne(-log(P)))
    {
      
      res <- boot.strength(expr.df,algorithm = 'chow.liu',cluster = clu)
      # pipeline(res)
      # title('Chow-Liu-Boot',line = -1)
    }
    
  })
}
stop('end')

{
  res <- boot(as.matrix(expr.dat),function(x,i){cor(x[i,])},R=600)
  C <- res$t0
  dev <- t(res$t)- as.vector(res$t0)
  SD <- apply(dev,1,sd)
  SD <- matrix(SD,dim(res$t0)) * sqrt(res$R)
  Z  =res$t0 / SD
  P = pnorm(abs(Z),lower.tail = F)
  # dat <- -log(P)
  dat <- Z
  diag(dat) <-0
  # res <- boot.strength(as.data.frame(dat),algorithm = 'chow.liu',cluster = clu)
  # res <- minet::minet(dat)
  # res <- minet::minet(abs(Z))
  # res <- cor(Z,method='spearman')
  res = abs(Z)
  pipeline(res)
}
with(env.sub,
     {
       out <- boot(as.matrix(expr.dat),function(x,i){cor(x[i,])},R=800)
       dev <- t(out$t)- as.vector(out$t0)
       SD <- apply(dev,1,sd)
       SD <- matrix(SD,dim(out$t0)) #* sqrt(out$R)
       C = out$t0
       Z  =out$t0 / SD
       P = pnorm(abs(Z),lower.tail = F)
       # dC = cov(t(dev))
       dC = cor(t(dev))
       dC[is.nan(dC)]<-0
       arr<-array(dC,rep(sqrt(nrow(dC)),4) )
       M <- apply(arr,c(1,3),mean)
       D <-apply(arr,c(1,3),sd)
       diag(Z) <- 0
    par(mfrow=c(3,3))
    # image(abs(C))
    image(M)
    pipeline(abs(C))    
    pipeline(M)

    mdl <- lm(as.vector(abs(Z))~as.vector(abs(C)))
    plot(abs(C),abs(Z))
    abline(mdl)
    mdl$residuals
    dC
    # dev
}) ->dC
image(dC)
dC[is.nan(dC)]<-0
arr<-array(dC,rep(sqrt(nrow(dC)),4) )
M <- apply(arr,c(1,3),mean)
D <-apply(arr,c(1,3),sd)
# which(is.na(dC),arr.ind = 1)
image(M/D)
with(env.sub,{
  image(C)
  image(M)
})
image(abs(dC[1:120,1:120]))
hist(dC)
