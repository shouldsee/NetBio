d.true <- distances(g.true)
d.mi <- make.mi(expr.dat)
d.mi.sp <-cor(expr.dat,method = 'spearman')%>%abs
pdist<-function(X,distF){
  D <- proxy::dist(X,method =distF,by_rows = F,diag = T)%>%as.matrix
  diag(D)<-apply(X,2,function(x)distF(x,x))
  D
}

{
  distF =  entropy::KL.empirical
  d.kl <- pdist(expr.dat,distF)
  # d.kl <- proxy::dist(expr.dat,method =distF,by_rows = F,diag = T)%>%as.matrix
  # diag(d.kl) <-
}

{
  genes<-randomWalk(adj.true.big,60,directed = F)
  genes <- union(genes,
                 randomWalk(adj.true.big,60,directed = F))
  # genes <- sample(genes.big,30,replace = F)
  # genes <- do.call(c,neighborhood(g.big,nodes=genes,order = 1))
  # true.pairs <- true.pairs.big
  adj.true <- adj.true.big[genes,genes]
  expr.dat<-expr.dat.big[,genes]
  plot(g.true <- igraph::graph_from_adjacency_matrix(adj.true)
       ,main=bquote(N==.(length(genes))))
  # length(genes)
}

# genes <- genes[]
subdata <- function(genes){
  with(.GlobalEnv,
       {adj.true <- adj.true[genes,genes]
       expr.dat<-expr.dat[,genes]
       g.true <- subgraph(g.true,genes)
       })
}

{
  distF <- transport::wasserstein1d
  d.wt <- pdist(expr.dat,distF)
}

routine.GENIE3<-function(dat,nCores=6,silent=1,...){
  require(GENIE3)
  require(doParallel)
  qc <- pipeline(GENIE3(t(dat),nCores=nCores,...),silent=silent)
  qc$method =  'GENIE3'
  .GlobalEnv$qc.list$GENIE3 <- qc
  qc$adjmat
}
routine.GENIE3(expr.dat)
# .GlobalEnv$
image(d.kl)
# summary.table(log(d.kl))
# image(-d.kl)
X <- -log(d.kl)
res <-d.kl
pipeline(res)
image(log(d.kl))
plot(density(t(tail(expr.mat,2))))
plot(density(tail(t(expr.mat),2)),ylim=c(0,1))
plot(density(tail(t(expr.mat),3)[1,]),ylim=c(0,1))
plot(expr.mat[,ncol(expr.mat)])
plot(expr.mat[,ncol(expr.mat)-1])
plot(expr.mat[,ncol(expr.mat)-2])
plot(expr.mat[,ncol(expr.mat)-3])
# plot(apply(expr.dat,2,var))
density.plot(expr.dat[,1],ylim=c(0,1),bw=.1)
density.plot(expr.dat[,2],ylim=c(0,1),bw=.1)
density.plot(expr.dat[,3],ylim=c(0,1),bw=.1)

density.plot <- function(x,bw='nrd0',...){
  plot(density(x,bw=bw,...),...)
}
# density.default()
# dens
# plot((apply(expr.mat,2,var)),type='b')


image(log(d.kl))
image(-g.dist)
make.mi


load.assignment.data()
subdata(genes[od==1])

{
  C<-cov(expr.dat)
  VAR <- diag(C)
  lVAR = log(VAR)
  max.lvar <- outer(lVAR,lVAR,pmax)
  mean.lvar <- outer(lVAR,lVAR,'+')
  inv.cov <- mat.invert(C)%>%abs
  l.inv.cov <- log(inv.cov)
  H <- apply(expr.dat,2,entropise)
  d.h <- outer(H,H,'+')
  d.maxh <- outer(H,H,pmax)
  d.mi <- make.mi(expr.dat)
  d.e2d <- pdist(expr.dat,entropy2d);diag(d.e2d) <- H
  d.kl <- pdist(expr.dat,entropy::KL.empirical)
  # d.e2d <- pdist(expr.dat,entropy2d);diag(d.e2d) <- H
  inv.d.e2d <- mat.invert(d.e2d,as.norm = T,post = abs)
  f <- function(mat) mat.invert(mat,post = abs)
  # pipeline(f(d.h/d.mi))
  pipeline(d.mi)
  pipeline(1/d.mi)
  # pipeline(f(d.e2d))
}
pipeline(d.mi)
pipeline(1/d.mi)
pipeline((f(d.e2d)))
pipeline((f(1/d.e2d)))
pipeline((dmat <- inv(d.kl)))
pipeline((dmat <- log(f(1-C))))
pipeline(log(f(d.mi)))
pipeline(inv(d.h/d.mi))
diag((d.mi)-diag(H)*2)
image(diag(1/d.mi))
dmat <-1/d.mi
C = cov(expr.dat)
inv = MASS::ginv
load.assignment.data()
pipeline(log(inv(cor(expr.dat))%>%abs))
pipeline(log(inv(C%>%abs)%>%abs))
pipeline(log(inv(1/d.e2d)%>%abs))
pipeline(log(inv(1/(d.mi-diag(H)))%>%abs))
pipeline(log(inv(1/(d.mi-diag(H)))%>%abs))
pipeline(log(dmat <- inv(d.h/d.mi)%>%abs))
pipeline(log(dmat <- inv(d.mi/d.maxh)%>%abs))

pipeline(dmat<-mat.invert(cor(expr.dat,method = 'pearson'),invertor = inv)%>%abs)
pipeline(dmat<-mat.invert(cor(expr.dat,method = 'spearman'),invertor = inv)%>%abs)
pipeline(-dmat)
pipeline(dmat<-cor(expr.dat,method = 'spearman')%>%abs)
pipeline(dmat<-log(inv(cor(expr.dat,method = 'spearman'))%>%abs))
pipeline(log(inv(cor(expr.dat,method = 'pearson'))%>%abs))
pipeline(dmat <-mat.invert((1/d.mi),invertor = solve,as.norm = F)%>%abs)
pipeline(dmat <-( mat.invert(exp(-d.mi),invertor = inv,as.norm = F)%>%abs))
pipeline(dmat <- (d.mi))
pipeline(dmat <- exp(-d.mi))
pipeline(dmat <- 1/(d.mi))

dmat <- exp(-d.mi)
diag(dmat)<-0
pmat <- apply(dmat,2,renorm)
pmat <- apply(pmat,2,function(x)x==max(x))
generator.mat <- apply(pmat,2,cumsum)

pmat <- dmat
install.packages.lazy('adaptMCMC')
library(adaptMCMC)
p.log <- function(g){pmat[,g]}

?MCMC
samp.1 <- MCMC(p.log, n=200, init=names(sample(exp(-H),1)),
               adapt=FALSE)
generator.mat <- apply(pmat,2,cumsum)
generator.mat <- apply(pmat,2,function(x)x==max(x))
init =sample(names(H),prob = exp(-H),1)


lookup<-function(r,cur,gen.mat){
  which(r<gen.mat[,cur])[1]
}
markov <- function(init.prob,generator.mat,nT=100){
  init =sample(seq_along(init.prob),prob = init.prob,1)
  rand <- runif(nT)
  out <- vector(length=nT)
  cur <- init
  for (i in seq_along(rand) ){
    out[i]=cur
    r = rand[i]
    cur <- lookup(r,cur,generator.mat)
  }
  out
}


samp <-lapply(1:100,function(x)markov(exp(-H),generator.mat))
mats <- lapply(samp,function(out)pair2adj(genes=genes,cbind(head(out,-1),tail(out,1))))
image(mats[[50]])
{
  sz= 1000
  dmat<-cor(expr.dat,method = 'spearman')
  dmat <- dmat
  dmat <- d.mi
  pmat <- apply(exp(-dmat),2,renorm)
  generator.mat <- apply(pmat,2,cumsum)
  
  # dmat <- exp(-d.mi)
  # diag(dmat)<-0
  # pmat <- apply(dmat,2,renorm)
  # pmat <- apply(pmat,2,function(x)x==max(x))
  # generator.mat <- apply(pmat,2,cumsum)
  
  samp <-lapply(1:sz,function(x)markov(exp(-H),generator.mat))
  mats <- lapply(samp,function(out)pair2adj(genes=genes,cbind(head(out,-1),tail(out,1))))
  
  M<-abind::abind(mats,along = 3)
  m0 <- apply(M,c(1,2),mean)
  mstd <- apply(M,c(1,2),sd)
  # m0 <- mats[[1]]
  # for (m in mats){
  #   m0 = m0 + m
  # }
  pipeline(m0)
  image(mstd)
  pipeline(mstd)
  pipeline(p<-pnorm(m0,sd=mstd))
}

{
  x <- samp[[1]]
  raiser <- (NN<-nrow(dmat))^(0:2)
  proj <- function(x) convolve(x-1,raiser,type = 'filter')
  out <- vector(length = 50^3)
  pjsamp<-sapply(samp,proj)
  pd <- table(as.integer(pjsamp))
  out[as.numeric(names(pd))+1] <- pd
  out  <- array(out,rep(50,3))
  
  adjs <- lapply(1:3,function(i)apply(out,setdiff(1:3,i),mean))
  adjs <- lapply(1:3,function(i)apply(out,setdiff(1:3,i),sd))
  
  pipeline(adjs[[1]]+adjs[[3]]-2*adjs[[2]])
  pipeline(adjs[[1]]*adjs[[3]]/adjs[[2]])
  pipeline(adjs[[1]]*adjs[[3]]/adjs[[2]])
}

image(apply(out,c(2,3),mean))
{
  pdn <-as.numeric(names(pd))
  x3 <- (pdn) %/% NN^2
  pdn <- pdn-x3*NN^2
  x2 <- (pdn) %/% NN^1 
  pdn <- pdn-x2*NN^1
  x1 <- (pdn) %/% 50^0
  # x2 <- x2-x3*50 
  # x1 <- x1-x2*50-x3*50^2 
  # x2<-x2 + 1
  # x1 <- x1 + 1
  df <- cbind(x1,x2,x3,pdn)
  df[,-ncol(df)] <- 1+  df[,-ncol(df)]
}
summary.table(df)
# summarise.t
table(df[,3])
df1 <- cbind(x1,x2,pdn)
df <- cbind(x1,x2,x3,pdn)
df[,1:3] <- df[,1:3] 

# df
groupColumns = c("x1","x2")
dataColumns = c("pdn")
# library(plyr)
# ?aggregate
sum1$x
df <- data.frame(df)
sum1 <-aggregate( df$pdn, df[,1:2], FUN = sum)
sum2 <-aggregate( df$pdn, df[,2:3], FUN = sum )
sum3 <-aggregate( df$pdn, df[,c(1,3)], FUN = sum )
sum1$pdn <- sum1$x
sum2$pdn <- sum2$x
sum3$pdn <- sum3$x
table(df$x2)
# res = plyr::ddply(df, groupColumns, function(x) colSums(x[dataColumns]))

DT <- data.table::data.table(df)
sum1 <- DT[, lapply(.SD,sum), by=list(x1,x2)]
sum2 <- DT[, lapply(.SD,sum), by=list(x2,x3)]
sum3 <- DT[, lapply(.SD,sum), by=list(x1,x3)]
with(sum1,cbind(x1,x2))
# pair2adj(genes=genes,as.matrix(sum1[,c('x1','x2')]),is.indexed = T,fill=sum1$pdn)
table(sum1$x1)
adj1 <- with(sum1,pair2adj(genes=genes,cbind(x1,x2),is.indexed = T,fill=pdn))
adj2 <- with(sum2,pair2adj(genes=genes,cbind(x2,x3),is.indexed = T,fill=pdn))
adj3 <- -with(sum3,pair2adj(genes=genes,cbind(x1,x3),is.indexed = T,fill=pdn))
image(adj1)
image(adj2)
image(adj3)
image(adj1 + adj2 + adj3)
xs <- bsxfun(cbind(pdn),t(cbind(raiser)),function(x,y){x%/%y/y})

# apply(-t(xs)*raiser,1,cumsum)
length(pd)/50^3
# cbind(pjsamp)

image(m0)
image(log(dmat))
diag(dmat) <-0
pnorm(m0,mstd)
ment <- apply(mats,c(1,2),entropise)
image(m0)
image(m0/mstd)
pipeline(m0/mstd)
pipeline(m0)
pipeline(dmat<-cor(expr.dat,method = 'spearman')%>%abs)
pipeline(dmat<-cor(expr.dat,method = 'pearson')%>%abs)

pipeline(minet::aracne())
image(m0)
dmat <- d.mi
diag(dmat) <- NaN
image((dmat))
out <- names(H)[out]

perfromance.pairs(genes,cbind(head(out,-1),tail(out,1)),true.pairs)


# markov.sampler <- function()
diag(dmat)<-NaN;pheatmap(dmat)
marg.eng<-function(x,y){
  # blanket <- neighbors(g.true,c(x,y),'total')
  blanket <-do.call(intersection,neighborhood(g.true,nodes=c(x,y)))
  # length(blanket)
  # print(blanket)
  if(length(blanket)==0){
    0
  }else{
    r <- d.e2d[blanket,x] + d.e2d[blanket,y] - H[blanket]
    mean(r)
  }
}
do.call(intersection,neighborhood(g.true,nodes=c(genes[1],genes[2])))
neighbors(g.true,genes[1],'in')
marg.eng(genes[1],genes[2])
marg.eng(genes[2],genes[1])
isSymmetric(d.e2d)
isSymmetric(d.eng)

df$d <- apply(df<-expand.grid(x=genes,y=genes),1,Fcompose(as.list,combine_args(marg.eng)))
d.eng <- matrix(df$d,dim(d.mi))
{
  par(mfrow=c(1,3))
  image(d.eng)
  dmat <- d.e2d
  diag(dmat)<-NaN
  image(dmat)
  # image(d.eng-d.e2d)
  image(-d.true)
}
d.mi <- make.mi(expr.dat)
pipeline(minet::minet(d.mi))
pipeline(minet::mrnet(d.mi))
pipeline(minet::aracne(d.mi))
pipeline(minet::mrnet(minet::build.mim(expr.dat)))
pipeline(minet::aracne(minet::build.mim(expr.dat)))
isSymmetric(d.eng)
image(d.eng-d.mi)
hist(d.eng)
hist(d.mi)
# marg.eng(genes[1],genes[9])
# pdist(rbind(genes[1:3]),marg.eng)
# marg.eng(genes[1],genes[2])
# pheatmap(d.mi)
plotnet(g.true)
dimnames(dmat) <- dimnames(d.mi)
e = eigen(M<-(dmat))
od <- order(apply(e$vectors*e$vectors,2,entropise),decreasing = T)
od <- rev(od)
e1 <- e$vectors[,gi<-od[1]]
e1 <- e$vectors[,od[2]]
e1 <- e$vectors[,od[3]]
e1 <- e$vectors[,od[4]]
e1 <- e$vectors[,od[5]]
e1sq <- e1^2
to = order(e1sq,decreasing = F)[1:5]
perfromance.pairs(genes,cbind(genes[gi],genes[to]),true.pairs)
# cumsum(order(log(e1sq),decreasing = T))sum(entropise(e1sq)/2
i<-which.max(e1^2);hist((e1^2)%>%log)
image((dmat<-outer(e1^2,e1^2,'*'))%>%log)
pipeline(dmat)
pipeline(dmat <- (outer(e1^2,e1^2,'*')))
pipeline(dmat <- (outer(e1,e1,'*')))

e1 <- e$vectors[,2]
pipeline(dmat <- (outer(e1,e1,'*')))

e1%*%M/e1
ev <-e$vectors
ev%*%
  norm(cbind(e$vectors[1,]),type = 'f')
diag(dmat)<-NaN
hist(d.mi)
dimnames(dmat) <- dimnames(d.mi);pheatmap(dmat)



plotnet <-function(g.true){
  par(mfrow=c(1,1))
  plot(g.true,vertex.size=0,edge.size=10,rescale=T,cex=1.25,
       edge.arrow.size=0,vertex.font.size=15
  )
}
plotnet(g.true )
routine.GENIE3(expr.dat,silent=0)
pipeline(log(inv(d.h/d.mi)%>%abs))
pipeline(log(inv(d.mi/d.e2d)%>%abs))
# diag(1/d.mi)
pipeline(chol2inv(chol((t(dmat)%*%dmat)))%>%abs%>%log)
image(chol2inv(chol((dmat)))%>%abs%>%log)
# ?is.positive.definite
image(dmat <- d.e2d)
image(1/d.mi)
image(dmat<-make.positive.definite(1/d.mi))
pipeline(abs(chol2inv(chol(dmat ))))
pipeline(log(abs(chol2inv(chol((C)) ))))
# image(-d.true)

chol2inv(chol(1/d.mi))
# ?chol2inv
pipeline(dmat <-f(1/((d.mi)-diag(H)))) ##### assert independence
pipeline(dmat <-f(1/d.mi)) ##### assert independence
load.assignment.data()
pipeline(f(d.h/d.mi))
pipeline(log(f(1/d.e2d)))
image(dmat);(d.true<-distances(g.true));image(-d.true)
routine.GENIE3(expr.dat,silent=0)
hist(dmat)
hist(qc.list$GENIE3$adjmat)
pipeline((dmat))
pipeline(log(resg+1E-3))
pipeline(log(exp(dmat)*resg+1E-3))
pipeline(resg <-qc.list$GENIE3$adjmat,dmat.prep = Fcompose(Fplus(1E-3),log))

pipeline(max.lvar)

image(distances(g.true))
# res <- (log(inv.cov)+max.lvar/2)
# res <- (max.lvar/2)
res <- (log(inv.cov)+mean.lvar+max.lvar)
res <- (log(inv.cov))
# res <- log(abs(C))+max.lvar/2
# res <- (abs(C))
# res <- ((C))
# res <- mean.lvar
{
  res <- l.inv.cov + mean.lvar/2
  res <-  mean.lvar
  res <-  mean.lvar
  res <- d.mi + outer(H,H,'+')
  pipeline(res)
}
routine.GENIE3(expr.dat,silent=0)
A = f(d.e2d)
# d.e2d <-
image(f(d.e2d))
image(inv.d.e2d <- mat.invert(d.e2d,as.norm = T,post = abs))
image(d.e2d)
image((d.e2d))
image(log(d.mi))
diag(d.mi) <- 0
d.mi
image(log(inv.d.e2d))
pipeline(f(exp(-d.mi)))

library(kernlab)
A<-d.e2d
# hclust(d = exp(-A))
# ?kernelMatrix
data(spam)
dt <- as.matrix(spam[c(10:20,3000:3010),-58])
image(log(dt+0.0000001))
# ?specc
clu <- hclust(d=as.dist(log(inv.d.e2d)),method = 'average')
{
  # dmat <- log(d.e2d)
  par(mfrow=c(1,3))
  # dmat <- d.h/d.mi
  # dmat <- d.maxh/d.mi
  dmat <- 1/d.e2d
  # dmat <- d.mi
  dmat <- log(f(dmat))
  dimnames(dmat)<-list(genes,genes)
  # dmat <- d.h/d.mi
  # dmat <- log(dmat)
  # dmat <- f(dmat)
  # dmat <- exp(-d.mi)
  # dmat <- d.mi/d.e2d
  # clu <- hclust(d=as.dist(dmat),method = 'average')
  # clu <- hclust(d=as.dist(dmat),method = 'complete')
  clu <- hclust(d=as.dist(dmat),method = 'ward.D2')
  # clu <- hclust(d=as.dist(dmat),method = 'median')
  # clu <- hclust(d=as.dist(dmat),method = 'mcquitty')
  # ?hclust
  plot(clu)
  diag(dmat)<-NaN
  image(dmat.od<-dmat[clu$order,clu$order])
  image(d.true[clu$order,clu$order])
}
# pipeline(dmat.od)
plot(g.true)
pipeline(dmat)

dimnames(dmat)<-dimnames(d.mi)
image(dmat.od)
{
  library(pheatmap)
  library(grid)
  library(gridExtra)
  library(gridGraphics)
  # install.packages.lazy('gridGraphics')
  par(mfrow=c(1,2))
  p1 <- pheatmap(dmat,Rowv = NA,Colv = NA,symm =T)
  # g <- grab_grob()
  
  par(mfrow=c(1,1))
  plot(g.true,vertex.size=0,edge.size=10,rescale=T,cex=1.25,
       edge.arrow.size=0,vertex.font.size=15
  )
  gridExtra::grid.arrange(p1,p2)
}
diag(dmat) <- 1
image(log(f(dmat)))
image(dmat.od[idx<-rSum > 0.6,idx])
image(log(f(dmat)))
pipeline(f(f(f(f(dmat)))))
pipeline((dmat))
plot(density(rSum<-rowSums(dmat[clu$order,clu$order],na.rm = T)/nrow(dmat)))
# plot(rSum)

od <- cutree(clu,k = 2)
# ?cutree
plot(clu)
{
  load.assignment.data()
  
  genes <- genes[od==3]
  
  subdata(genes)
}

load.assignment.data()
genes <- genes[od==2]

routine.GENIE3(expr.dat,silent = 0)

load.assignment.data()
image(dmat[od==1,od==1])
image(dmat[od==2,od==2])
image(d.true[od==1,od==2]==1)
image(d.true[od==1,od==1]==1)
image(d.true[od==2,od==2]==1)
image(d.true[od==3,od==3]==1)

pipeline(f((1/d.mi)))

pipeline(f(f((d.mi)))%>%log)
pipeline(f(d.e2d))
# pipeline(f(exp(-d.mi)))
# plot(exp(-d.mi),1/d.mi)
# cut(0.03,clu)
# image(d.true)
image(d.true[clu$order,clu$order])

ridx <- sample(1:nrow(d.true),replace = F)
image(d.true[ridx,ridx])
cutree(clu,h=0.03)
load.assignment.data()

pipeline(log(f(d.mi/d.e2d)))
pipeline(log(f(d.mi)))
plot(clu)
identify(clu)
image(log(d.e2d)[clu$order,clu$order])
# lea
# ?reorder
reorder(clu)
order.dendrogram(clu$labels)
pipeline(log(f((cov(expr.dat)))))
pipeline(log(inv.d.e2d))
hist(log(inv.d.e2d))
image(d.r<-array(runif(prod(dim(d.mi))),dim(d.mi)))
# image(exp(d.r));pipeline(exp(d.r))
pipeline(max.lvar)
plot(density(upper.tri.get(f(d.e2d))))
# pipeline(d.r)
diag(d.mi) <- H
pipeline(log(f(d.mi)))
# pipeline(f(cov(expr.dat)))
pipeline(qc.list$GENIE3$adjmat,dmat.prep = Fcompose(Fplus(1E-3),log))
plot(sort(log(inv.d.e2d+1E-3)))
plot(sort(qc.list$GENIE3$adjmat))
routine.GENIE3(expr.dat,silent=0,nCores = 12)
# hist(d.e2d)
mat.invert(d.e2d)
diag(d.e2d)
pipeline(d.mi)
# image(f(d.e2d))
pipeline(res)
pipeline(d.e2d)
pipeline(f(d.e2d))
pipeline(f(d.mi))
entropise(expr.dat[,1])
f(d.e2d)
entropy2d(expr.dat[,1],expr.dat[,2],eps=0)%>%print
d.e2d
eigen(f(d.e2d))
heatmap(f(d.e2d))

image(d.e2d)
pipeline(f(d.e2d))
pipeline(d.e2d)
pipeline(d.mi)

load.assignment.data()



routine.GENIE3(expr.dat,silent=0)
Fplus <- function(c){
  function(x){x+c}
}
pipeline(qc.list$GENIE3$adjmat,dmat.prep = Fcompose(Fplus(1E-3),log))
image(d.e2d)
image(res)

routine.GENIE3(expr.dat)


ld.var <- outer(lVAR,lVAR,'-')%>%abs
inv.cov <- mat.invert(cov(expr.dat))%>%abs
res <- (log(inv.cov))
pipeline(res)
pipeline(max.lvar)

res <- (log(inv.cov)+max.lvar/2)
pipeline(res)
res <- (log(inv.cov)+max.lvar)
pipeline(res)

load.assignment.data()
# setwd()
pipeline(qc.list$GENIE3$adjmat)
image(res)
'Variance dependent partial correlation'
pipeline(res)
image(log(inv.cov))

image(log(inv.cov)-d.var)
image(inv.cov)
image(d.var)
image(-g.dist)
plot(lVAR)

# eigen(d.kl)
# matplot(t(tail(expr.mat,3)),type='b')
matplot(((expr.dat[1:100,2:4])),type='b')
matplot(((expr.dat[1:100,ncol(expr.dat)+1-1:5])),type='b')
matplot(((expr.dat[1:100,ncol(expr.dat)+1-1:2])),type='b')


# image(-d.wt)
image(d.kl)
image(log(d.kl))
# image(log(d.mi)) ### low-mi = large distance
image(log(d.mi.sp)) ### low-mi = large distance
image(t(d.kl[1:10,1:10]))
image(1-exp(-d.mi)) ### hi-mi = low prob
image(log(d.mi)) ### hi-mi = low prob
image(-g.dist)


install.packages.lazy('transport')
# image(d.kl) ### hi-mi = low prob
# d.mi
image(g.dist)

image(-g.dist==-1)
X <- log(1-d.kl)
pc <- prcomp(X)
plot(pc)
s = summary(pc)
# image(s$x)
pca(log(d.kl))
# diag(X) <-
res <- f(X)
res <- (d.kl)
pipeline(res)
pipeline(qc.list$GENEIE3$adjmat)
par(mfrow=c(1,3))
image(log(d.kl))
image(-d.true)
image(-log(d.mi))
image(1-d.mi.sp%>%log)
plot(density(d.mi),ylim=c(0,1))
plot(density(d.mi.sp),ylim=c(0,1))

amat2igraph <- function(amat,weighted=T,...){
  graph_from_adjacency_matrix(amat,weighted = weighted,...)
}

par(mfrow=c(1,3))
image(d.true)
entropy::KL.empirical(1:10,10:1)
# image(distances(amat2igraph((d.mi)))
image(1-d.mi.sp%>%log)


