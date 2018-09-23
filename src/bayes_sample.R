# options(error=recover)
# qc <- pipeline(res)

rbinom
# dmultinom
# rmultinom
# bn <- graph
source('dag.R')
aM=all.graphs[[5]]
aM = rbind(c(0,1),
           c(0,0)
           )
aM.dep <- aM
aM.idep<-  matrix(0,2,2)
dat <- dat1[,1:2]

states <- c()
bn <- amat2bn(aM)
sim.dat<-bnlearn::rbn(bn,10,data=dat)
# sim.dat
array2array(sim.dat)
plot(c(0,1),c(0,1))
points(sim.dat[,1],sim.dat[,2])
plot(bn)

{
bin2dec <- function(dat){
  if(!is.array(dat)){
    callback<-function(x){
      if (is.factor(x)){
        as.integer(x)-1
      }else{
        as.integer(x)
      }
    }
    dat<-array2array(sim.dat,callback = callback)
  }
  raiser<-2^(1:ncol(dat)-1)
  rowSums(bsxfun(dat,t(cbind(raiser)),'*'))
}


dec2bin <- function(dat)
  {
  m <- sapply(dat, function(x){ as.integer(intToBits(x))})
  md <- ceiling(log2(max(dat)))
  t(m[1:md,])
}
}

sim.dat.dec <- bin2int(sim.dat)
int2bin(sim.dat.dec)
dat <- int2bin(rep(0:3,each=10))
bn.dat <-as.data.frame(array2array(dat,callback = as.factor))
colnames(bn.dat) <- c('A','B')
bfit <- bn.fit(bn,bn.dat)
sim.dat<-bnlearn::rbn(bn,1000,data=bn.dat)
hist(bin2int(sim.dat),0:10-.5)

aM <- matrix(0,2,2)
aM <- aM.dep
bn <- amat2bn(aM)
bfit <- bn.fit(bn,bn.dat)
sim.dat<-bnlearn::rbn(bn,1000,data=bn.dat)
sim.dat<-bnlearn::rbn(bn,1000,data=table2df(sess$tb))
hist(bin2int(sim.dat),0:10-.5)
data.frame(sess$tb)
table2df <- function(tb){
  df <- data.frame(tb)
  ridx <- rep(1:nrow(df),df$Freq)
  df[ridx,-ncol(df)]
}

sess$dat <-data.frame(A=1,B=0)
# sess$make_table()
{
  sess$tb <-as.table( matrix(0,2,2,dimnames = list(A=c(0,1),B=c(0,1))))
  sess$tb[1,1] <- 1
  # sess$tb[2,2] <- 1
  sess$mdlgraph <- igraph::graph_from_adjacency_matrix(aM)
  sess$logL()
}

m <- matrix(1,2,2)

pmat <- prop.table(m)
{
dep_sampler <- function(pmat,n=1){
  pmat <- as.table(matrix(pmat,2,2))
  # if (!is.matrix(pmat)){
  #   pmat <- as.matrix(pmat)
  # }
  p <- pmat
  stopifnot(abs(sum(p)-1)<1E-8)
  # t(rmultinom(n,size=1,prob=p))
  t(rmultinom(1,n,prob=p))
}
idep_sampler <- function(pmat,n=1,as.p =F){
  pmat <- as.table(matrix(pmat,2,2))
  rowmarg <- margin.table(pmat,1)
  colmarg <- margin.table(pmat,2)
  p <- outer(rowmarg,colmarg,'*')
  if (as.p){
    return(p)
  }
  stopifnot(abs(sum(p)-1)<1E-8)
  # t(rmultinom(n,size=1,prob=p))
  # t(rmultinom(1,n,prob=p))
  t(rmultinom(1,n,prob=p))
}

loopy_sampler <- function(pmat,n=1){
  if (!is.table(pmat)){
    pmat <- as.table(pmat)
  }
  rowprop <- prop.table(pmat,1)
  colprop <- prop.table(pmat,2)
  p <- rowprop * colprop
  p <- p/sum(p)
  stopifnot(abs(sum(p)-1)<1E-8)
  # t(rmultinom(n,size=1,prob=p))
  t(rmultinom(1,n,prob=p))
}
}

loopit <-function(pmat){
  pmat <- as.table(matrix(pmat,2,2))
  rowprop <- prop.table(pmat,1)
  colprop <- prop.table(pmat,2)
  p <- rowprop * colprop
  p <- p/sum(p)
  as.vector(p)
}

overcor <- function(pmat,dat,as.mi = T){
  pmat <- as.table(matrix(pmat,2,2))
  rowprop <- prop.table(pmat,1)
  colprop <- prop.table(pmat,2)
  lpmat <- log(rowprop) + log(colprop) - log(pmat)
  # pmat <- rowprop * colprop/pmat
  sum(lpmat[(as.matrix(dat)=="1")+1])
}


MI <-function(p){
  pmat <- as.table(matrix(p,2,2))
  rowmarg <- margin.table(pmat,1)
  colmarg <- margin.table(pmat,2)
  sum(bsxfun(log(rowmarg),t(log(colmarg)),'+') -log(pmat))
}
# log(rowmarg)
# dim(t(log(colmarg)))
MI(p)

{
  {
    pmat <- runif(4)
    pmat <- pmat/sum(pmat)
    # pmat <- matrix(pmat,2,2)
  }
  pold <- pmat
  for ( i in 1:20)
  {
    dlogL =overcor(pmat,dat)
    print(dlogL)
    print(pmat)
    print(logL)
    if (dlogL>0){
      pold <- pmat
      pmat<-loopit(pmat)
      logL=sum(log(pmat[(as.matrix(dat)=='1')+1]))
    }else{
      p <- pmat[c(3:4,1:2)]
      dlogL =overcor(p,dat)
      print(bquote(dlogL==.(dlogL)))
      if(dlogL<0){
        break
      }else{
        pmat <-p
      }
      # break
      # pmat <- t(matrix(pmat,2,2))
      # pmat <- idep_sampler( pmat,as.p = T)
      
    }
  }
  par(mfrow=c(1,2))
  p <- pold
  # p <- pmat
  logL=sum(log(p[(as.matrix(dat)=='1')+1]))
  print(bquote(logL==.(logL)))
  # p <- pa
  # plot(as.table(matrix(dep_sampler(p,200),2,2)))
  MI.cur <- MI(p)
  print(bquote(MI==.(MI.cur)))
  MI.true <- MI(prop.table(table(dat)))
  print(bquote(TRUEMI==.(MI.true)))
  plot(as.table(matrix(p,2,2)))
  plot((table(dat)))
}
matrix(p,2,2)

# library(GeneNet)
# installed.packages()
plot(as.vector(dep_sampler(pa,100)),type='b')
pa <- pmat
pmat
states <- matrix(1:4,2,2)
# states[(as.matrix(dat)=='1')+1]
# plot(unlist(as.table(dec2bin(dep_sampler(pa,100)))))
# dev.off()
{
  par(mfrow=c(1,2))
  p <- pmat
  # p <- pa
  # plot(as.table(matrix(dep_sampler(p,200),2,2)))
  plot(as.table(matrix(p,2,2)))
  plot((table(dat)))
}
# a = hist(bin2dec(array2array(dat,as.numeric)),1:11-.5)

overcor(pmat,dat)
{
  par(mfrow = c(1,2))
  # p <- 1:4
  ps <-list(c(1,1,1,0.9),
           c(1,1,1,1.1))
  lapply(ps,
         {function(p) {
           p0<-p
    plot(0,0,xlim=c(0,5),ylim=c(0,1))
    p <- p/sum(p);for (i in 1:20){
    lines(p,col=i)
    points(p,col=i)
    p<-loopit(p)
    }
    title(bquote(
      p==.(paste(round(p0,3),collapse=','))))
    # paste(round(p,3),collapse=','))
         }
  })
  dev.copy(jpeg,filename="overcorrelation.jpg")
  # dev.print()
  dev.off()
}

# install.packages.lazy('cubature')
# library(cubature)
loss <- function(y){
  y0 = c(0.25,0.125,0.375,0.25)
  Metrics::mae(loopit(y),y0) + abs(sum(y)-1) + sum(-pmin(y,0))
}
{
  y0=runif(4)
  y0 = y0/sum(y0)
  res <- optim(y0,loss
             ,method='BFGS'
             ,control = list(maxit=1000)
             )
  print(res$par)
  print(res$value)
}
res
loopit(res$par)
cubature::hcubature(0,1)
# {
#   png() # on a Unix-alike
#   plot(rnorm(10), main = "Plot 1")
#   dev.copy(device = png,file='test.png')
#   mtext("Copy 1", 3)
#   dev.print(width = 6, height = 6, horizontal = FALSE) # prints it
#   dev.off(dev.prev())
#   dev.off()
# }

plot(1,1)
pmat
margin.table(pmat,1)
dimnames(pmat)<-list(A=c(0,1),B=c(0,1))


{
  par(mfrow = c(1,3))
  # p<-c(1,1,1,1)
  # p <- 1:4
  p <- 1:4
  # p <- c(1,4,3,2)
  # plot()
  # N=7000
  pmat <- matrix(p,2,2)%>%prop.table
  pmat%>%print
  # pmat <- matrix(c(1,2,4,3),2,2)%>%prop.table
  sim.dat <-dep_sampler(pmat,N)
  plot(1:4,sim.dat,type='l')
  sim.dat <-idep_sampler(pmat,N)
  lines(1:4,sim.dat,col=2)
  sim.dat <-loopy_sampler(pmat,N)
  lines(1:4,sim.dat,col=3)
  # bindat <- 
}
# image(as.matrix(sess$tb))
prop.table(sess$tb)
# exp(-1.386294 - -2.120264)
# is(sess)

# hist(bin2int(dat),0:10-.5)
# hist(bin2int(sim.dat),0:10-.5)
# hist(rowSums(bsxfun(array2array(sim.dat,callback = as.numeric),t(cbind(2^(1:ncol(sim.dat)-1))),'*')))
# decompose()

# binarise.aggF
# hist(sim.dat)

# bnlearn::rbn
# ?rbind_list
# ?bnlearn::rbn
m<-methods(rbn,class='default')
# print(m)
attr(m,'info')
class(bn) <- 'bn.fit'
bnlearn:::rbn.backend
bnlearn:::call_rbn_master
methods("axis")
