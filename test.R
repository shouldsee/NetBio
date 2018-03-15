# pair2adj()
# image(dmat)
image(pmat)
install.packages.lazy('expm')

colSums(inv(pmat%*%pmat))

ip <- inv(pmat)
# sqrtm
pa <- inv(ip%*%ip)
pb <-expm::sqrtm(pmat)


t(pb) %*% (pb)


# ?expm::sqrtm
pb

image(d.e2d)
entropy2d(X[,1],X[,1])

{
  d.e2d <- pdist(expr.dat,partial(entropy2d,nbin=c(50,50)))
  d.hx <- cbind(diag(d.e2d))
  # d.hx <- cbind(apply(expr.dat,2,entropise))
  d.H <- bsxfun(t(d.hx),d.hx,'+')
  d.Hmax <- bsxfun(t(d.hx),d.hx,pmax)
  d.CH <-bsxfun(d.e2d,d.hx,'-') 
  d.mi <- bsxfun(t(d.hx),d.CH,'-')
  image(d.CH)
  image(d.mi)
  d.mi.norm <- d.mi/d.Hmax
}

mat.view(d.mi.norm)
mat.view(d.mi/d.H)
hist(d.mi.norm)
which.max(d.mi.norm)
image(d.mi.norm)
image(d.H - d.e2d)
d.mi.norm <- 2*d.mi/d.H

# diag(d.H - d.e2d)
plot(cov(expr.dat)%>%abs,d.mi)

image(abs(cov(expr.dat)))
{
  d.e2d <- cov(expr.dat)
  d.hx <- sqrt(cbind(diag(d.e2d)))
  d.CH <-bsxfun(d.e2d,d.hx,'/')  %>% abs
  d.mi <- bsxfun(d.CH,t(d.hx),'/')
  image(log((d.CH)))
  image(log(d.mi))
}
mat.view(d.e2d)
mat.view(d.mi)

image(d.CH)
# d.H.Y_X <- bsxfun(d.e2d,d.hx,'-') 
diag(d.e2d) - d.hx
image(d.e2d)
image(d.H.Y_X)
# pa <-inv(ip%*%ip%*%ip%*%ip%*%ip%*%ip%*%ip%*%ip%*%ip%*%ip%*%ip)
# colSums(pa)
# image(pa %*% pa - pmat)
# pipeline(pmat)
# pa <- combine_args(Fcompose)((mat.sqrt,2))(pmat)
# make.mi
# image(d.mi)
mat.view<-function(mat){
  image(mat)
  bquote(Diagonal==.(diag(mat)))
}

mat.view(d.mi.norm)
mat.view(d.mi)
{

  dmat <- d.mi.norm
  # dmat <- d.mi
  # praw <- exp(-dmat)
  # praw <- exp(dmat)
  praw <- dmat
  pmat <- apply(praw,2,renorm)
  Z <- matrix(apply(praw,2,sum),nrow=1)
  pmat <- dmat
  # pmat <- exp(-dmat)
  generator.mat <- apply(pmat,2,cumsum)
}
image(pmat)
{
  Nit=4
  outlst <- Reduce(function(x,y)mat.sqrt(x),
                c(list(pmat),rep(1,Nit)),
                accumulate = T)
  # p <- pa + pmat
  # p <- apply(pa,2,renorm)
  # pa <- abs(pa)
  # pa <- pa ^ 2
  # pa <- pa/sum(abs(pa))
  # p <- pa + pmat
  out.raw <- abind::abind(outlst,along = 3)
  NORM <- apply(out.raw,c(3),function(m) norm(m,'F'))
  # lapply()
  # out <-apply(out.raw,c(1,2),function(m){m<-(m^2);m/sum(m)})
  out <- lapply(outlst,function(m){m<-(m^2);m/sum(m)})
  # out <- lapply(out.raw,function(m){m<-abs(m);m/sum(m)})
  out <- abind::abind(out,along = 3)
  # MEAN <- apply(out,c(1,2),median)
  # MEAN <- apply(out,c(1,2),mean)
  
  MEAN <- apply(out,c(1,2),min)
  # MEAN <- apply(out,c(1,2),mm)
  STD <- apply(out,c(1,2),sd)
  # p <- abs(pa)
  # p <- MEAN/STD
  # p <- MEAN
  p <- abs(outlst[[3]])
  # p <- out[,,5]
  # apply(out,)
  # p <- STD
  # p <- bsxfun(p,Z,'*')
  # p <- pa
  # p <- pa
  pipeline(p)
}

{
  NORM <- apply(out.raw,c(3),function(m) norm(m,'F'))
  NORM <- apply(out.raw,c(3),function(m) sqrt(sum(m^2)))
  # NORM.F <- apply(out.raw,c(3),function(m) norm(m,'F'))
  plot(log(NORM.F))
  NORM.1 <- apply(out.raw,c(3),function(m) norm(m,'1'))
  plot(log(NORM.1))
}

plot(NORM.F,NORM.1,log='xy')
abline(0,1)
image(out.raw[,,2])
image(out.raw[,,2]^2)
image(out[,,4])
# dim(adj.true)
hist(d.mi.norm[adj.true])
hist(d.mi)


{
  g.dist <- distances(g.true)
  dmat <- adj.true
  # dmat <- exp(-g.dist)
  praw <- dmat
  pmat <- array2array(dmat,as.numeric)
  # pmat <- pmat*d.mi.norm
  pmat <- pmat*(1-d.mi.norm)
  # pmat <- apply(praw,2,renorm)
  # Z <- matrix(apply(praw,2,sum),nrow=1)
  image(log(pmat%^%2))
  image(log(pmat%^%3))
  image(log(pmat%^%4))
  image(log(pmat%^%5))
  image(log(pmat%^%6))
  image(as.symmetric(log(pmat%^%20)))
  image(log(EM <- expm(pmat)))
  # colSums()
}
image(log(EM))

image((dmat*d.mi.norm))
image(-log(d.mi.norm))
mat.view(-log(d.mi.norm))
pipeline(log(d.mi.norm))
colSums(EM)
pmat
library(expm)


image(log(pmat%^%5))
image(log(as.matrix(expm(pmat))))
image(d.mi)
colSums(pmat)

d.mi.norm[1:3,]
{
  mat.view(d.mi)
  mat.view(1-d.mi.norm)
  pmat <- apply(exp(d.mi.norm),1,renorm)
  pmat <- d.mi.norm
  # pmat <- d.mi/d.H*2
  # pmat <- 1/d.e2d
  # pmat <- d.e2d/d.Hmax
  
  # pmat <- d.mi.norm
  # pmat <- 1/d.mi
  pmat <- 1-d.mi.norm
  # pmat <- 1-d.mi.norm
  # M <- logm(pmat)
  # M <- logm(pmat,'Eigen')
  M <- OpenMx::logm(pmat)
  M <- OpenMx::logm(M)
  # diag(M) <- max(M)
  diag(M) <- NaN
  # hist(M)
  pipeline(abs(M))
}

library(expm)
library(Rutil)
library(igraph)
g.true <- graph_from_adjacency_matrix(adj.true)
source('dream5.R')
{
  par(mfrow=c(1,3))
  image(log(expm(array2array(adj.true,as.numeric))))
  image(log(expm(exp(-distances(g.true)))))
  dmat <- array2array(adj.true,as.numeric)
  dmat[!adj.true] <- Inf
  image(dmat.t <- log(expm(exp(-dmat))))
}
image(distances(g.true))

image(expr.dat){
  M <- expm(1-d.mi.norm)
  diag(M) <- NaN
  image(log(M))
}
pipeline(-M)

image(-d.mi)

mean(d.mi.norm)
diag(-d.mi.norm)
image(-d.mi.norm)
image(log(abs(cov(expr.dat))))


pipeline(log(abs(mat.invert(cor(expr.dat),as.norm=F))))
pipeline(log(abs(mat.invert(cov(expr.dat),as.norm=T))))
pipeline(log(abs(cov(expr.dat))))
pipeline((abs(mat.invert(cov(expr.dat)))))

M <- abs(mat.invert(d.mi.norm,as.norm = T)) 
diag(M) <- NaN
image(log(M))

M <- cov(expr.dat)
# image(log(abs(M)))
pipeline(abs(M))
pipeline(log(abs(M)))
a<-(abs(mat.invert(M,as.norm = T)));pipeline(log(a))

{
  M <- cov(expr.dat)
  par(mfrow=c(1,3))
  image(dmat.t)
  # image(cor(expr.dat,method='spearman'))
  image(d.mi.norm)
  # image(log(d.mi))
  # image(cor(expr.dat,method='pearson'))
  image(log(abs(M)))
}
M <- d.mi.norm
image(M<-d.CH/d.Hmax)
image(log(abs(D<-M%^%.25)))
image(D)
diag(M) <-  0;
diag(D)<-NaN
pipeline(M)

plot(log(abs(M)),dmat.t)
image(log(abs(M)))
abs(mat.invert(M,as.norm = F)) %>%pipeline


image(log(abs(mat.invert(d.mi,as.norm=F))))
pipeline(log(abs(chol2inv(d.mi))))
(log(abs(mat.invert(d.mi.norm,as.norm = F))))%>%pipeline

# D <- dist(d.mi.norm,)%>%as.matrix
D <- dist(d.mi.norm,"maximum")%>%as.matrix
D <- dist(D,"maximum")%>%as.matrix
image(D)
# pmat <- exp(1-d.mi.norm) 
{
  # pmat <- exp(-d.mi.norm)
  # pmat <- exp(-d.mi)
  # pmat <- 1/d.mi.norm
  # pmat <- cor(expr.dat,method='spearman')
  pmat <- d.CH / d.Hmax
  pmat <- abs(pmat)
  
  # pmat <- cor(expr.dat)
  # pmat <- d.mi.norm
  # pmat <- 1/d.mi
  # pmat <- d.mi
  M <- pmat
  Z <- apply(M,2,sum)
  M <- apply(M,2,renorm)
  # image(OpenMx::logm(M) )
  par(mfrow=c(2,3))
  Niter = 10
  lst = vector('list',Niter)
  for (i in 1:Niter){
    # image(log(abs(M) ))
    image(log(abs(M*Z) ))
    M <- M/(mean(abs(M)))
    lst[[i]]=abs(M)
    title(i)
    M <- (mat.sqrt(M))
    # image(M)
    # M <- OpenMx::logm(M)
  }
  image(dmat.t)
  plot(sapply(lst,sd),type='b')
}

pmat <- 1-cor(expr.dat)
M <- pmat
Z <- apply(M,2,sum)
M <- apply(M,2,renorm)
pa <- mat.sqrt(M)
rowSums(M)
pb <- pa
abs(expm::sqrtm(pa))%>%colSums
{
  pa <-pb
  e <- eigen(pa)
  pb <- mat.sqrt(pa)
colSums(pb)
pa
e
pb
}
pb

pa%^%0.5
colSums(pb)
pb <- mat.sqrt(pb)
colSums(pb)

# diag(pmat)
# pipeline(mat.sqrt(mat.sqrt(mat.sqrt(cov(expr.dat)))))
pipeline(mat.sqrt(cov(expr.dat)))
pipeline(d.mi.norm)
# ou
pmin(lst[[1]],lst[[2]] , lst[[3]]  ,lst[[4]])%>%log%>%pipeline
pipeline((abs(lst[[1]]*Z)))

pipeline((-log((lst[[3]])) ))
array(runif(length(genes)^2),dim(d.mi))%>%pipeline
pipeline((log((lst[[5]])) ))
pipeline((-log((lst[[5]])) ))
pipeline((-log((lst[[4]])) ))
# sd(lst[[3]])
# sd(lst[[5]])
pipeline((-log(lst[[5]])))
pipeline((log(lst[[6]])))
pipeline((log(lst[[7]])))
pipeline((log(lst[[10]])))
pipeline(d.mi.norm)
pipeline(-(t(d.CH)+d.CH))

load.assignment.data()
{
  pipeline(minet::aracne(d.mi.norm))
  pipeline(minet::aracne(minet::build.mim(expr.dat)))
  pipeline(minet::mrnet(minet::build.mim(expr.dat)))
  pipeline(minet::minet(expr.dat))
}
plot(g.true)
pipeline(minet::mrnet(d.mi))
pipeline(minet::mrnet(d.mi.norm))
resg<-routine.GENIE3(expr.dat,silent=0)
pipeline((resg))
pheatmap::pheatmap(expr.dat)

install.packages.lazy('NMF')

f <- NMF::nmf(t(expr.dat),prod(dim(expr.dat))^0.5 %>%floor)
# image(f@fit@W)
library(pheatmap)
H <- rowSums(f@fit@H)
W <- t(f@fit@W)*H

D = dist(f@fit@W,'euclidean')%>%as.matrix
D = dist(f@fit@W,'euclidean')%>%as.matrix
VAR = apply(expr.dat,1,sd)
MEAN = apply(expr.dat,1,mean)
df <- reshape2::melt(expr.dat,varnames=c('obs','gene'),value.name='expr')
df$gene <- factor(df$gene)
df$obs <- factor(df$obs)

mdl = lm(expr~gene,df)
summary(mdl)
plot(g.true)
res <- GeneNet::ggm.estimate.pcor(expr.dat)
res <- GeneNet::network.test.edges(res)
D <-pair2adj(res[,c('node1','node2')],genes=genes,fill=res$pval)
# pipeline((D))
(cov(expr.dat)%>%abs%>%log)%>%image
image(-log(D))
image(distances(g.true))

install.packages.lazy('WGCNA')
image(dmat.t)
par(mfrow=c(1,3))
pheatmap(expr.dat)
load.assignment.data()
pipeline(cor(expr.dat)%>%abs)
pipeline(cor(expr.dat,method='spearman')%>%abs)
pipeline(cor(expr.dat,method='spearman')%>%abs%>%log)
image(C<-cov(expr.dat)%>%abs%>%log)
image(dmat.t)

pipeline(C)
plot(g.true)
clu <- hclust(dist(t(expr.dat),'euclidean'))
clu <- hclust(dist(t(expr.dat),'pearson'))
{
  par(mfrow=c(2,3))
  D <- cor(expr.dat,method='spearman')%>%abs%>% as.dist
  D <- 1-D
  clu <- hclust(D,method='ward.D2')
  # clu <- hclust(1-(cor(expr.dat,method='pearson') ))
  # clu <- hclust((cor(expr.dat,method='pearson') ))
  
  K = 6
  tr <- cutree(clu,h=1.0)
  
  # tr <- cutree(clu,k = K)
  m <- init.amat(genes)
  od <-order(tr)
  gps <- split(od,tr[od])
  for (k in 1:max(tr)){
    gi <- which(tr==k)
    gs <- genes[gi]
    # m[expand.grid(gi,gi)]=T
    m[gi,gi]=T
  }
  image(m[clu$order,clu$order])
  image(dmat.t[clu$order,clu$order])
  image(adj.true[clu$order,clu$order])
}


{
  par(mfrow=c(2,3))
  idx <- gps[[1]]
  dat <- expr.dat[,idx]
  gres<-GENIE3(t(dat),nCores = 6)
  plot(subgraph(g.true,idx))
}

{
  par(mfrow=c(2,3))
  image(adj.true[idx,idx])
image(dmat.t[idx,idx])
image(g.dist[idx,idx])
image(minet::minet(dat))
image(as.matrix(1-D)[idx,idx])
# image(gres)
image(as.symmetric(gres%>%log))
# image(as.symmetric(gres))
}
image(d.mi.norm)
pipeline(D<-cor(expr.dat,method='spearman')%>%abs)

C <- cor(expr.dat,method='spearman')
p <- apply(D,2,renorm)
image(mat.sqrt(p))
(m<-p) %>% image

(m*upper.tri(m))%>%log%>%image
mat.sqrt(p)%>%image
# mat.sqrt(p)%>%mat.sqrt
m <- expm::sqrtm(p)%>%abs%>%sqrtm
m <- expm::sqrtm(D)%>%abs
m <- expm::sqrtm(C)%>%abs
(m*upper.tri(m))%>%image
m <- expm::sqrtm(m)




# (m*upper.tri(m))%>%abs %>%log%>%image
pipeline(m*upper.tri(m))
# pipeline(mat.sqrt(D))
# pipeline(-d.mi.norm)
summary(d.mi.norm%>%as.vector)
image(
  (g.dist)
  (dmat.t[idx,idx])
)
{
  # pipeline(minet::aracne(d.mi.norm))
  pipeline(minet::aracne(minet::build.mim(dat)))
  pipeline(minet::mrnet(minet::build.mim(dat)))
  pipeline(minet::minet(dat))
}

image(dmat.t)

# image(dmat.t[clu$order,clu$order])
pipeline(m)
hist(dmat.t[m])
sum(m&adj.true) / sum(m)
sum(m&adj.true) / sum(adj.true)
# perfromance.pairs()

image(m)
pipeline(m)
image(adj.true)

# ?pheatmap
# ?heatmap
pheatmap(dmat.t)
image(dmat.t)
image(dmat.t[clu$order,clu$order])

# order(tr)
# for (k in 1:K){
#   which(tr==k)
# }

plot(clu)

clu$order
# ?pheatmap
# load.assignment.data()
ADJ <- WGCNA::adjacency(expr.dat)
tom <- WGCNA::TOMdist(ADJ)
pipeline(-tom)
pipeline(ADJ)
cov(expr.dat)%>%abs%>%pipeline
image(tom)
?WGCNA::TOMdist
WGCNA::TOMdist
# image(ADJ+tom)
WGCNA::plotNetworkHeatmap(df)
pipeline(-df)
expr.dat[,'G199']%>%hist
# lm(expr.dat)
plot(MEAN,VAR,xlim=c(0,11),ylim=c(0,2))
hist(expr.VAR)
pipeline(-D)
# pipeline(rma)
# ?NMF::rmatrix
D <- cov(W) 
D <- dist(t(W),'manhattan')%>%as.matrix
image(log(D))
# pipeline(cor(W)%>%abs)
# pipeline(log(D%>%abs))
# pipeline(abs(cor(expr.dat)))
mat.invert(D)%>%abs%>%pipeline
pipeline(D%>%abs)
pipeline(cov(expr.dat)%>%abs)
mat.invert(D)%>%abs%>%pipeline
image(abs(D)%>%log)
image(-distances(g.true))
image(dmat.t)

pdist(t(f@fit@W),'euclidean')
pheatmap(f@fit@H)
# f@fit
# f@fit@W
load.assignment.data()
# pipeline(-(lst[[1]]))
sdiag(dmat.t)
# colSums(expm(exp(-dmat)))
image(log(expm(dmat)))
# pipeline(log(expm(array2array(adj.true,as.numeric))))
image(log(mat.sqrt(d.mi.norm)%>%abs))
image(log(abs(cov(expr.dat))))
image(log(d.mi.norm))
image(d.mi)
image(log(d.mi/d.Hmax))
image(d.e2d/d.Hmax)
image(log(expm(array2array(adj.true,as.numeric))))
image(exp(-distances(g.true)))
image(OpenMx::logm(expm(array2array(adj.true,as.numeric))))
image(OpenMx::logm(expm(exp(-distances(g.true)))))
image(OpenMx::logm(exp(-distances(g.true))))
pipeline(OpenMx::logm(exp(-distances(g.true))))
image(pmat)


install.packages.lazy('OpenMx')


diag(M) <- NaN
image(M)
pipeline(M)
# ?logm

image((apply(d.mi.norm,2,renorm) ))
# ?logm
dim(out)
MEAN
# image(MEAN)

image(pa)
colSums(pa)


b = pa%*%pa
image(pa)

pa
image(pa)

image(pmat)

colSums(mat.sqrt(pa))
# mat.sqrt
# rowSums(pa)
colSums(pa)

load.assignment.data()
image(as.matrix(expr.dat))
pheatmap::pheatmap(as.matrix(expr.dat))
M = dist(expr.dat,)%>%as.matrix
diag(M) <- NaN
# pipeline(-M)
M%>%image

hist(unlist(expr.dat))
resg<-routine.GENIE3(expr.dat,silent=0)
image((resg))
image(expm(exp(-dmat)))

