# datadir <- 'dream5/training data/Network 3 - E. coli/'
res.list <- list()
{
  fname = 'dream5/training data/Network 3 - E. coli/net3_expression_data.tsv'
  f <- readLines(fname)
  dat <- read.csv(fname,sep='\t')
  expr.dat.big <- as.matrix(dat)
  
  
  fname <- 'dream5/test data/DREAM5_NetworkInference_GoldStandard_Network3 - E. coli.tsv'
  dat <- read.csv(fname,sep='\t')
  true.pairs.big<-dat[dat[,3]==1,1:2]
  genes.big <- colnames(expr.dat.big)
  adj.true.big <- pair2adj(true.pairs.big ,genes = genes.big)
}

{
  fname = 'dream5/custom/InSilicoSize50-Ecoli2_proteins_multifactorial.tsv'
  f <- readLines(fname)
  dat <- read.csv(fname,sep='\t')
  expr.dat <- as.matrix(dat)
  
  fname <- 'dream5/custom/InSilicoSize50-Ecoli2_goldstandard.tsv'
  dat <- read.csv(fname,sep='\t')
  true.pairs<-dat[dat[,3]==1,1:2]
  genes <- colnames(expr.dat)
  adj.true <- pair2adj(true.pairs ,genes = genes)
  g.true <- graph_from_adjacency_matrix(adj.true)
}



dim(adj.true.big)
{
  genes<-randomWalk(adj.true.big,30,directed = F)
  genes <- sample(genes,replace = F)
  # genes <- sample(genes.big,30,replace = F) 
  # genes <- do.call(c,neighborhood(g.big,nodes=genes,order = 1))
  # true.pairs <- true.pairs.big
  adj.true <- adj.true.big[genes,genes]
  expr.dat<-expr.dat.big[,genes]
  plot(g.true <- igraph::graph_from_adjacency_matrix(adj.true)
       ,main=bquote(N==.(length(genes))))
  # length(genes)
}

normalise <- function(x){
  (x-mean(x))/sd(x)
}
renorm <- function(x){
  x/sum(x)
}



{
  pc <- prcomp(t(expr.dat))
  par(mfrow=c(1,2))
  biplot(pc)
  plot(apply(expr.dat,2,var)%>%sort)
}


links <- adj.true

# library(MRNet)
install.packages.lazy('MRNET')
library(GeneNet)
genes <- colnames(expr.dat)


# funcion
read.tsv <- function(fname,...){
  read.csv(fname,sep='\t',...)
}
# BiocInstaller::biocLite('MRnet')
{
  {
    # fname='dream5/custom/InSilicoSize50-Ecoli2_multifactorial_perturbations.tsv'
    fname='dream5/custom/InSilicoSize50-Ecoli2_proteins_multifactorial.tsv'
    # fname='dream5/custom/InSilicoSize50-Ecoli2_nonoise_multifactorial.tsv'
    dat <- read.tsv(fname)
    expr.dat <- as.matrix(dat)
  }
  
  {
    fname='dream5/custom/InSilicoSize50-Ecoli2_goldstandard.tsv'
    dat <- read.tsv(fname)
    true.pairs<-dat[dat[,3]==1,1:2]
  }
  genes <- colnames(expr.dat)
  adj.true <- pair2adj(true.pairs,genes = genes)
}




NPT=200
hist(as.vector(expr.dat))
{
  # dat <- log(expr.dat-min(expr.dat)+0.1)
  dat <- t(apply(t(expr.dat),2,renorm))
  # dat <- expr.dat
  # dat <- apply(exp)
  PKGMethod='GeneNet'
  fig.cap=PKGMethod
  PKG = strsplit(PKGMethod,'\\.')[[1]][1]
  library(PKG,character.only = T)
  # ?ggm.estimate.pcor
  # ?network.test.edges
  pcor.dyn <- ggm.estimate.pcor(dat)
  df <- network.test.edges(pcor.dyn,direct = T,plot=F)
  pairs = as.matrix(df[,c('node1','node2')])
  # score = 
  res <- pair2adj(pairs,genes=genes,is.indexed = T,symmetric = F
                  # ,fill= df$pval
                  # ,fill= -df$qval.dir
                  # ,fill = -log(df$pval)
                  # ,fill = -log(df$qval)
                  # ,fill=sqrt(df$pval)
                  # ,fill = log(df$prob.dir+0.0001)
                  ,fill = -log(df$qval.dir+0.0001)
                  # ,fill = -log(df$pval.dir+0.0001)
  )
  print('QC started')
  qc <- pipeline(res,NPT)
  qc$method =  PKGMethod
  qc.list[[PKG]] <- qc
}

require(Rutil)
# Rutil::compose

library(ggplot2)

entropise <- entropy
{
  VAR <- apply(expr.dat,2,var)
  ENT <- apply(expr.dat,2,entropise)
  
  odvar = order(VAR,decreasing = T)
  od = order(ENT,decreasing = T)
  idx <-1:5
  # idx <-9:16
  print(ENT[od[idx]])
  expr.df <- as.data.frame(expr.dat[,od[idx]])
  expr.df <- reshape2::melt(expr.df,variable.name='gene',value.name='expr') 
  ggplot(expr.df) + geom_density(aes(x=expr,color=gene,y=1+..count..)) + xlab('expression') + 
    NULL #scale_y_log10()
}


# dst <- rDPGibbs(Data = list(y=expr.dat),Mcmc = list(R=500),Prior = list())
library(bayesm)
library(DPpackage)
clu =makeCluster(9)
install.packages.lazy('DPpackage')

{
  deg <- neighborhood.size(g.big,nodes = genes,order=2)
  plot(ENT,log(deg))
  title(lm2eqn(mdl<-lm(log(deg)~ENT),add.rsq = T))
  abline(mdl,col=2)
}

plot(VAR,ENT)
points(VAR[id],ENT[id],col=2)
plot(sqrt(VAR),log(deg))
id <- which.max(sqrt(VAR)-ENT)
# plot(density(expr.dat[,id]))

install.packages.lazy('entropy')
x<-expr.dat[,od[1]]

ds <- density(x)
plot(ds)
sum(ds$y)

# ds$entropy
# diff(ds$x)

entropy <- function(x){
  ds <- density(x)
  eps = 0
  ent = -sum(ds$y*log(ds$y+eps))*diff(ds$x[1:2])
  # ent <- do.call(integrate,c(f=approxfun(ds$x,ds$y*log(ds$y+eps)),as.list(range(ds$x))))
}


{
  # dat <- log(expr.dat)
  dat <- expr.dat
  # dat <- log(expr.dat-min(expr.dat)+0.1)
  # dat <- t(apply(t(expr.dat),2,renorm))

  PKGMethod='GENIE3'
  fig.cap=PKGMethod
  PKG = strsplit(PKGMethod,'\\.')[[1]][1]
  library(PKG,character.only = T)
  library(doParallel)
  set.seed(0)
  res <- GENIE3(exprMatrix = t(as.matrix(dat)),nCores = 10,verbose = T)
  res[is.nan(res)]=0
  print('QC started')
  qc <- pipeline(res,NPT)
  qc$method =  PKGMethod
  qc.list[[PKG]] <- qc
}




# expr.dat.std <- apply(expr.dat,2,normalise)
{
  # dat <- expr.dat.std
  dat <- expr.dat
  PKGMethod='xgboost.importance'
  cap=PKGMethod
  PKG = strsplit(PKGMethod,'\\.')[[1]][1]
  library(PKG,character.only = T)
  # source('xgb.R')
  df <- Gxgb.fit(dat
                 ,model.dir = 'qc/',max_depth=length(genes)%>%floor
                 ,alpha=2.5/10
                 ,lambda=2.5/10)
  pairs = as.matrix(df[,c('Feature','output')])
  VAR <- apply(expr.dat,2,var)
  res <- pair2adj(pairs,genes=genes,is.indexed = F,symmetric = F,
                  fill =df$Gain)
                  # fill =df$Gain*VAR[df$output])
                  # fill = df$Gain*VAR[df$output]*VAR[df$Feature])
  qc <- pipeline(res)
  qc$method =  PKGMethod
  qc.list[[PKG]] <- qc
}
library(bnlearn)
plot(amat2bn(as.symmetric(adj.true)))
# xgb.plot.multi.trees(mdls[[1]])
sum(df$Gain[df$output==genes[2]])
df[!complete.cases(df),]
# expr.dat

# g.true <- graph_from_adjacency_matrix(adj.true,)
# distance_table(g.true, directed = F)
# g.dist <-distances(g.true,)
# res <- pair2adj(pairs,genes=genes,is.indexed = T,symmetric = F
#                 ,fill = -log(df$prob.dir))
# res <- pair2adj(pairs,genes=genes,is.indexed = T,symmetric = F
#                 ,fill = -log(df$prob))
# 


{
  library(igraph)
  g.big <- graph_from_adjacency_matrix(adj.true.big)
  g.big.dist <- distances(g.big)
  
  
}
plot(g.sub.dist,res)
plot(g.dist,res)


lres <- -log(res)
lres.dist <- distances(graph_from_adjacency_matrix(1-res,weighted = T))
lres.dist <- distances(graph_from_adjacency_matrix(lres,weighted = T))
heatmap((expr.dat))
marg.mi <- function(expr.dat,mimat=make.mi(expr.dat)){
  
}

d.mi <- make.mi(expr.dat)
image(d.mi)
hist()
image(log(abs(X)))
image(-g.sub.dist)
image(-g.dist)

{
  {
    par(mfrow=c(1,2))
    g.sub.dist<-g.big.dist[genes,genes]
    g.dist <- distances(g.true)
    res.tsd=(res)
    # res.tsd=-log(1-res)
    # res.tsd <- -log(res)
    # res.tsd=abs(res)
    g.res <- graph_from_adjacency_matrix(res.tsd,weighted = T)
    plot(g.res)
    res.dist <- distances(g.res)
    xs = upper.tri.get(g.sub.dist)
    # xs <- 1-exp(-xs)
    # ys = upper.tri.get(-log(1-res))
    ys = upper.tri.get(res.dist)
    
  }

  {
    idx = is.finite(ys)
    mdl = lm(ys[idx]~xs[idx])
    # ys[is.infinite(ys)]<- 10
    plot(xs,
         ys,ylim=c(0,min(10,max(ys[is.finite(ys)]))),
         xlab='distance on true network',
         ylab='distance on predicted network')
    abline(mdl,col=2)
    title(lm2eqn(mdl,add.rsq = T))
  }
  abs(cor(xs,ys,method = 'spearman'))
}

{
  
  C = cov(expr.dat)
  COR <- cor(expr.dat)
  P = chol2inv(C)
  pdg <- sqrt(diag(P))
  pdg <- outer(pdg,pdg,'*')
  pcor <- abs(P/pdg)
  # image(abs(C/out))
  image(log(abs(COR)))
  image(log(pcor))
  image(log(abs(COR)))
  image(log(abs(C)))
  image(as.matrix(as_adjacency_matrix(g.true)))
  image(log(abs(C)))
  d <- d.mi
  # d <- exp(-d.mi)
  # d <- log(abs(pcor)+0.0001)
  d <- d.mi
  plot(xs<-upper.tri.get(g.sub.dist),
       ys<-upper.tri.get(d))
  l =lm(as.vector(ys)~as.vector(xs))
  title(lm2eqn(l,add.rsq = T))
  image(d.mi)
  image(g.sub.dist)
  image(-g.dist)
  image(-g.sub.dist+g.dist)
  plot(g.sub.dist,C)
  plot(g.sub.dist,pcor)
  dis.dat <- discretize(as.data.frame(expr.dat))
}
plot(density(d.mi[,1]))
{
  g1 = 16
  g2 = which.max(d.mi[-g1,g1])
  g2 = g2 + (g2 >= g1)
  plot(expr.dat[,c(g1,g2)])
}
df.mi <- reshape2::melt(d.mi)
od.mi <- order(df.mi$value,decreasing = T)

{
  df.mi <- reshape2::melt(X)
  od.mi <- order(df.mi$value,decreasing = T)
  
  par(mfrow=c(1,3))
  for (i in 1:3){
    l=df.mi[od.mi[i],]
    g1 = l[[1]];g2=l[[2]]
    plot(expr.dat[,c(g1,g2)])
    title(bquote(MI==.(l[[3]])))
  }
}

identity(1)
iv <- chol2inv(d.mi)
sdg<-(diag(iv))
sdg <- (sdg)
niv <- iv/outer(sdg,sdg,'*')
niv.abs <- abs(niv)
image(log(niv.abs))
pipeline(res=niv.abs)
image(log(d.mi))
image(-g.sub.dist)
image(log(niv))
pipeline(res=d.mi)
pipeline(res=niv.abs)
d.mi.partial <- niv.abs

# pipeline(d.mi.partial)
d.mi
d.cov <- cov(expr.dat)
image(log(abs(d.mi)))
image(log(abs(d.cov)))
f <- partial(mat.invert,post=abs)
plot(f(d.mi),f(d.cov))
{
  pipeline((mat.invert(cor(expr.dat),post=abs)))
  f = function(mat){(mat.invert(mat,post=function(x)x^2))}
  f = function(mat){(mat.invert(mat,post=function(x)abs(x)))}
  plot(g.true)
  image(distances(g.true)%>%as.matrix)
  X <- cov(expr.dat)
  X = abs(X)
  diag(X) <- 1
  pipeline((X))
  pipeline(Fcompose(f)(X))
  res <- f(X)
  plot(upper.tri.get(adj.true),upper.tri.get(res))
  # plot(g.true)
  
  X = make.mi(expr.dat)
  pipeline((X),dmat.prep = log )
  pipeline(Fcompose(f)(X),dmat.prep = log)
  pipeline(Fcompose(f,f)(X))
  pipeline(Fcompose(f,f,f)(X))
  pipeline(Fcompose(f,f,f,f)(X))

  X <- minet::build.mim(expr.dat);diag(X)<-1
  # pipepipeline()
  qc.list$mrnet$adjmat <- minet::mrnet(make.mi(expr.dat))
  res <- minet::mrnet(minet::build.mim(expr.dat))
  res <- minet::minet(minet::build.mim(expr.dat))
  pipeline(res,dmat.prep = abs)
  pipeline(res,dmat.prep = Fcompose(function(x)x+0.1,log))
  pipeline(qc.list$GENIE3$adjmat)
  pipeline((normalise(res)+normalise(qc.list$GENIE3$adjmat)))
  pipeline(qc.list$mrnet$adjmat)
  # pipeline(Fcompose(f,f,f,f,f,f,f,f,f,f)(cov(expr.dat)))
  # pipeline(mat.invert(mat.invert(mat.invert(cor(expr.dat),post=abs),post=abs),post=abs))
  pipeline(mat.invert(cov(expr.dat)))
  pipeline(mat.invert(d.mi,post=abs))
  pipeline(qc.list$GENIE3$adjmat)
  
  pipeline(qc.list$GENIE3$adjmat <- GENIE3(t(expr.dat),nCores = 6))
}
# pipeline(minet::minet(expr.dat),dmat.prep=identity)
# pipeline(d.mi)
pipeline(qc.list$GENIE3$adjmat)
row.names(qc.list$GENIE3$adjmat)

pipeline(niv.abs)
pipeline(qc.list$GENIE3$adjmat)
{
  options(error=recover)
  pipeline(qc.list$GENIE3$adjmat)
  options(error=stop)  
}
pipeline(qc.list$GENIE3$adjmat)
pipeline
nan.fill=0
npt=200
res[is.nan(res)]<-nan.fill
thres.lst <- select.cutoff(res,npt)

hist(d.mi)
BiocInstaller::biocLite('minet')

y2d = rbind( c(1,2,3), c(6,5,4) )
entropy::mi.empirical
entropy::freqs.empirical(expr.dat)
y1 = c(4, 2, 3, 1, 10, 4)
y2 = c(2, 3, 7, 1, 4, 3)
{
  idx = c(1,15)
  entropy::mi.empirical(expr.dat[,idx])
}
entropy::mi.empirical(t(expr.dat[,1:10]))
dist(expr.dat,method= function()entropy::mi.empirical )

entropy::entropy.empirical(expr.dat[,1],unit = 'log')
{
  entropise <- function(x,eps=1E-3){
    ds <- density(x)
    # eps = 
    ent = -sum(ds$y*log(ds$y+eps))*diff(ds$x[1:2])
    # ent <- do.call(integrate,c(f=approxfun(ds$x,ds$y*log(ds$y+eps)),as.list(range(ds$x))))
  }
  ds =density(rbind(x,x+1))
  entropy2d <- function(x,y,nbin=c(25,25),eps=1E-3){
    ds <- MASS::kde2d(x,y,n = nbin)
    # eps = 0
    dx <- diff(ds$x[1:2])
    dy <- diff(ds$y[1:2])
    ent <- -sum(ds$z*(log(ds$z+eps)))*dx*dy
  }
  d.e2d <- proxy::dist(expr.dat
                       ,method=partial(entropy2d,nbin=c(30,30))
                       ,diag = T,by_rows = F) %>%as.matrix
  expr.h <- apply(expr.dat,2,entropise)
  d.h <- outer(expr.h,expr.h,'+')
  d.mi = d.h - d.e2d 
  # image(d.mi)
  image(log(d.mi))
}

entropise(expr.dat[,1])%>%print
vectorized_pdist <- function(A,B) {
  an = apply(A, 1, function(rvec) crossprod(rvec,rvec))
  bn = apply(B, 1, function(rvec) crossprod(rvec,rvec))
  
  m = nrow(A)
  n = nrow(B)
  
  tmp = matrix(rep(an, n), nrow=m) 
  tmp = tmp +  matrix(rep(bn, m), nrow=m, byrow=TRUE)
  sqrt( tmp - 2 * tcrossprod(A,B) )
}

install.packages.lazy(c('proxy','pdist'))
library(pdist)

f2d <- function(f){
  function(x,y) f(rbind(x,y))
}
d.mi <- proxy::dist(expr.dat,method = f2d(entropy::mi.empirical),by_rows = F,diag = T) %>% as.matrix
image(d)


f <- partial(array2array,callback=as.numeric)
im <- apply(dis.dat,2,Fcompose(factor,as.numeric))
heatmap(im,distfun = partial(dist,method='minkowski'))


##### Regress the distance function 