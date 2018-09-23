
# source('header.R')
{
ADJ = !init.amat(genes)
C <- cor(expr.dat, method='spearman')
C.abs <- abs(C)
D <- 1 - C.abs
geneis <- seq_along(genes)

gd <- expand.grid(g1=geneis,g2=geneis,g3=geneis)%>%as.matrix
# gd <- expand.grid(g1=geneis,g2=geneis)%>%as.matrix
isValid <- ADJ[gd[,c(1,2)]] & ADJ[gd[,c(2,3)]] & ADJ[gd[,c(1,3)]]
gd <- gd[isValid,]
nrow(gd) %>% print


m <- C.abs
m <- D
m <- 1-m
# m <- -log(m)
# m <- 1/m
# dep <- abs(logm[gd[,c(1,3)]] - (logm[gd[,c(1,2)]] + logm[gd[,c(3,2)]])/2)
# dep <- logm[gd[,c(1,3)]] - (logm[gd[,c(1,2)]] + logm[gd[,c(3,2)]])/(logm[gd[,c(1,2)]] - logm[gd[,c(3,2)]])
# dep <- 1/m[gd[,c(1,3)]] - (1/m[gd[,c(1,2)]] + 1/m[gd[,c(3,2)]])/2
# dep <- m[gd[,c(1,3)]] - (m[gd[,c(1,2)]] + m[gd[,c(3,2)]])/2
s13 <- m[gd[,c(1,3)]]
s12 <- m[gd[,c(1,2)]]
s23 <- m[gd[,c(3,2)]]
dep <- m[gd[,c(1,3)]] - (m[gd[,c(1,2)]] + m[gd[,c(3,2)]])/2
nor <- m[gd[,c(1,3)]] + (m[gd[,c(1,2)]] + m[gd[,c(3,2)]])/2
dep <- dep / nor
dep <- abs(dep)

df.raw <- cbind(gd,s13=s13,s12=s12,s23=s23)
# pairs <- cbind(gd[,1],gd[,2],dep)
pairs <- cbind(gd,dep)
idx <- pairs[,1]!=pairs[,3] & pairs[,1]!=pairs[,2] & pairs[,3]!=pairs[,2]
pairs <- pairs[idx,]
df.raw <- df.raw[idx,]


hasEdge <- adj.true[pairs[,c(1,3)]]
pairs <- cbind(pairs,hasEdge)
pairs <- pairs[order(pairs[,'dep'],decreasing = T),]%>%as.matrix
hist(pairs[,'dep'],20)
DT <- data.table::data.table(pairs)
DT.raw <- data.table::data.table(df.raw)
}


{
  # BINS=seq(-8,8,length.out=30)
  BINS=seq(-1,1,length.out=30)
  par(mfrow=c(2,3))
  for(i in sample(geneis,6)){
    j = sample(geneis[-i],1)
    # hist(DT[g1==j & g3==i,]$dep,BINS,main='')
    df <- DT.raw[g1==j & g3==i,]
    plot(x=df$s12,y=df$s23,main='',xlim=c(0,1),ylim=c(0,1),xlab=i,ylab=j)
    title(bquote(.(as.character(adj.true[j,i]))~~.(df$s13)))
  }
}


image(-distances(g.true)[3:5,3:5])
image(-distances(g.true)[1:19,1:19])
image(distances(g.true))
g.true.big <- graph_from_adjacency_matrix(adj.true.big)
image(distances(g.true)[genes,genes]-distances(g.true.big)[genes,genes])
image(-distances(g.true.big)[genes,genes])

hist(DT[g1==1 & g3==8,]$dep,BINS)
hist(DT[g1==1 & g3==9,]$dep,breaks=BINS)
?hist
# hist(DT[g1==1 & g3==2,]$dep)
hist(DT[g1==1 & g3==3,]$dep)

image(t(lst))
load.assignment.data()
length(genes)
{
  # res <- pair2adj(sum1[,c(1,2)]%>%as.matrix,is.indexed = T,fill=sum1[,'dep'],genes=genes)
  ps <- tester(param,as.pair=T)%>%as.matrix
  res <- pair2adj(ps[,c(1,2)]%>%as.matrix,is.indexed = T,fill=ps[,3],genes=genes)
  # pipeline(res)
  # pipeline(res^2)
  pipeline(-res)
}


pipeline(C.abs)
pipeline(-C.abs)
image(C.abs)
image(-distances(g.true))
gres <- routine.GENIE3(expr.dat,silent=0)

param = c(0,1,0,0,0,.00000001)
param = c(1,0,0,0,0,.00000001)
param = c(1,0,0,0,0,0)
param = c(0,1,0,0,0,0)
param = c(0,0,1,0,0,0)
param = c(0,0,0,1,0,0)
param = c(0,0,0,1,1,1)
plot(param)


plot(g.true)
pipeline(gres)
pipeline(-gres)
param <- ores$par
param <- tail(lst,3)[1,]
image(lst%>%t)
pipeline(-res)
pipeline(C.abs)
pipeline(-C.abs)

tester<-function(w,l1 = 0.5,as.pred=F,as.pair=F){

  # w <- w/sum(abs(w))
  w <- w
  # agg <- function(x){ w[1]*max(x) + w[2]*min(x) + w[3]*median(x) + w[4]*mean(x) + w[5]*sd(x)}
  agg <- function(x){ (w[1]*max(x)+w[2]*min(x)+w[3]*median(x) + w[4]*mean(x))/(sd(x)*w[5]+1*(1-w[5])) }
  # agg <- function(x){ (w[1]*max(x)+w[3]*median(x) + w[4]*mean(x))/(sd(x)*w[5]+1*(1-w[5])) }
  # agg <- function(x){ w[1]*max(x) + w[2]*min(x) + w[3]*median(x) + w[4]*mean(x) + w[5]*sd(x)
  #   + w[6] * quantile(x,probs=0.85)}
  # agg <- function(x){ median(x)/sd(x)}
  # agg <- function(x){ abs(median(x)/sd(x))}
  
  # agg <- function(x){ max(x)/sd(x)}
  # agg <- function(x){ w[1]/-1.5*(-1.5*max(x) + 0.5*min(x))  + w[4]*mean(x) + w[5]*sd(x)}
  sum1 <- DT[, lapply(.SD, agg), by=list(g1,g3)]
  if(as.pair){ return(sum1[,c('g1','g3','dep')])}
  if(w[6]){sum1 <- abs(sum1)}
  hasEdge <- adj.true[sum1[,c(1,2)]%>%as.matrix]
  pd<-hasEdge[order(sum1[,'dep'])]
  if(as.pred){return(pd)}
  -get_AROC(pd) + l1*mean(abs(w))
}

ores <- optim(w0,tester,method = 'BFGS',control = list(maxit=10))
ores <- optim(ores$par,tester,method = 'BFGS',control = list(maxit=10))

image(lst%>%t)
param <- lst[nrow(lst),]

{
  # ores$par <- runif(5)
  # ores <- optim(ores$par,tester,method='BFGS',control = list(maxit=10))
  ores <- optim(ores$par,tester,control = list(maxit=50))
  # ores%>%print
  tester(ores$par,l1=0)%>%print
  lst <- rbind(lst,ores$par)
}

plot(g.true)

nudge <- function(x,r=0.1){
  x + r*(runif(length(x))-0.5)*mean(par0)
}

ores$par <- nudge(ores$par)
ores$par <- nudge(par0)
ores$par <- runif(5)
# parlst <- rbind()
parlst<- rbind(parlst,ores$par)
######## Finding best parameter
lst <- rbind()
# apply(parlst,1,tester)
for(i in 1:10){
  ores$par <- runif(5)
  # ores <- optim(ores$par,tester,method='BFGS',control = list(maxit=10))
  ores <- optim(ores$par,tester,control = list(maxit=100))
  # ores%>%print
  tester(ores$par,l1=0)%>%print
  lst <- rbind(lst,ores$par)
}

image(t(lst))
param <- tail(lst,5)[1,]
{
  # param <- parlst[13,]
  # param <- apply(lst,2,mean)
  # param <- apply(lst,2,median)
  pd <- tester(param,as.pred = T)
  print(param)

  # pd <- sum1[,'hasEdge']
  ys <- (cumsum(pd))/sum(pd)
  xs <- seq_along(pd)/length(pd)
  AROC <- abs(sum(ys-xs))*diff(xs[1:2])
  plot(xs,ys)
  title(bquote(.(AROC)))
  # plot(seq_along(),cumsum(pd)/sum(pd))
  abline(0,1,2)
}

image(lst%>%t)

par(mfrow=c(2,3))
pheatmap::pheatmap(lst)

parlst <- rbind(parlst,lst)
image(lst)
param <- lst[4,]
tester(param,l1=0)
{
  # sum1 <- DT[, lapply(.SD,min), by=list(g1,g3)]
  # sum1 <- DT[, lapply(.SD,max), by=list(g1,g3)]
  # sum1 <- DT[, lapply(.SD,median), by=list(g1,g3)]
  # sum1 <- DT[, lapply(.SD,function(x)(max(x))/sd(x) ), by=list(g1,g3)]
  # sum1 <- DT[, lapply(.SD,function(x)(max(x)-min(x))*sd(x) ), by=list(g1,g3)]
  # sum1 <- DT[, lapply(.SD,function(x)sd(x)), by=list(g1,g3)]
  # sum1 <- DT[, lapply(.SD,function(x)median(x)/sd(x)), by=list(g1,g3)]
  # sum1 <- DT[, lapply(.SD,sum), by=list(g1,g3)]
  # sum1 <- DT[, lapply(.SD,sum), by=list(g1,g3)]
  # w <- w/sum(abs(w))
  
  agg <- function(x){ w[1]*max(x) + w[2]*min(x) + w[3]*median(x) + w[4]*mean(x) + w[5]*sd(x)}
  sum1 <- DT[, lapply(.SD, agg), by=list(g1,g3)]
  hasEdge <- adj.true[sum1[,c(1,2)]%>%as.matrix]
  pd<-hasEdge[order(sum1[,'dep'])]
  # AROC <- abs(mean((cumsum(pd)-seq_along(pd))/sum(pd)/sum(pd)))
  # sum1 <- arrange(sum1,desc(dep))
  # sum1$hasEdge <- adj.true[sum1[,c(1,2)]%>%as.matrix]
  # sum1 <- arrange(sum1,desc(dep))
  # head(sum1,50)[,'hasEdge']%>%mean%>%print
  # tail(sum1,100)[,'hasEdge']%>%mean
  # head(sum1,100)[,'hasEdge']%>%mean
  # tail(sum1,800)[,'hasEdge']%>%mean
  # tail(sum1,nrow(sum1)%/%4)[,'hasEdge']%>%mean
  # tail(sum1,50)[,'hasEdge']%>%mean
  
  {
    # pd <- sum1[,'hasEdge']
    ys <- (cumsum(pd))/sum(pd)
    xs <- seq_along(pd)/length(pd)
    AROC <- abs(sum(ys-xs))*diff(xs[1:2])
    plot(xs,ys)
    title(bquote(.(AROC)))
    # plot(seq_along(),cumsum(pd)/sum(pd))
    abline(0,1,2)
  }
}
routine.GENIE3(expr.dat,silent=0)
pd <- sum1[,'hasEdge']
get_AROC <- function(pd){
  ys <- (cumsum(pd))/sum(pd)
  xs <- seq_along(pd)/length(pd)
  AROC <- abs(sum(ys-xs))*diff(xs[1:2])
}

par(mfrow=c(2,3))
pd <- sum1[,'hasEdge']
plot(cumsum(pd)/sum(pd))
head(sum1,800)[,'hasEdge']%>%mean%>%print
# ys <- sapply(xs,function(x)head(sum1,x)[,'hasEdge']%>%mean)
plot(1:20*20)
tail(sum1,700)[,'hasEdge']%>%mean%>%print

{
  trunc <- floor(nrow(sum1)*0.1)
  # trunc <- floor(nrow(sum1)*0.9)
  df <-head(sum1,trunc)
  # df <-tail(sum1,-trunc)
  idx <- df[,c(1,2)]%>%as.matrix
  ADJ[idx] = F
}

nrow(sum1)
idx <- tail(sum1,100)[,c(1,2)]%>%as.matrix
idx <- tail(sum1,100)[,c(1,2)]%>%as.matrix
ADJ[idx] = F
tail(sum1,100)$hasEdge%>%mean
# tail(sum1,50)$hasEdge%>%mean
hist(pairs[,'dep'])
hist(sum1$dep)



load.assignment.data()
# tail(sum1,200)[,'hasEdge']%>%mean

resg <- routine.GENIE3(expr.dat,silent=0)
# 
# image(res)
# image(g.dist)
# pipeline(-res)
# 
# head(sum1,50)[,'hasEdge']%>%mean
# hist(pairs[,'dep'])
# # order_by(sum1,dep)
# sum1 <- arrange(sum1,desc(dep))
# sum1$hasEdge <- sum1$hasEdge/length(genes)
# dat <- head(sum1,80)
# hist(pairs[,4])
# # tail(pairs,80)
# mean(dat$hasEdge)
# 
# 

# # ?order_by
# # ?arrange
# head(sum1)
# 
# pairs <- pairs[order(dep,decreasing = T),]%>%as.matrix
# perf <- partial(perfromance.pairs,genes= genes,true.pairs=true.pairs)
# ps <- head(pairs,2)
# 
# head(pairs,20)
# 
# adj.true[t(cbind(ps[1:2,1:2]))]
# 
# hist(dep)
