
load.assignment.data()



{
  require(plyr)
  func <- function(xx)
  {
    # return(data.frame(COR = cor(xx$s12, xx$s23,method='spearman')))
    # return(data.frame(COR = cor(xx$s12, xx$s23,method='spearman')))
    # return(data.frame(COR = max(abs(xx$s12-xx$s23))))
    # return(data.frame(COR = min(abs(xx$s12-xx$s23))))
    # return(data.frame(COR = max(abs(xx$s12*xx$s23))))
    # return(data.frame(COR = mean(xx$s12+xx$s23)))
    # return(data.frame(COR = median(xx$s12*xx$s23)))
    return(data.frame(COR = cor(xx$s12, xx$s23,method='pearson')))
    # return(data.frame(COR = cov(xx$s12, xx$s23)))
  }
  
  
  
}

routine.bnlearn.bootstrap <- function(expr.dat,iss=1,R=200,score='bde',...){
  require(bnlearn)
  
  # expr.dat.bin <- as.data.frame( binarise.aggF(expr.dat,median))
  expr.dat.bin <- bnlearn::discretize(as.data.frame(expr.dat),method = 'quantile',breaks=3)
  # expr.dat.bin <- bnlearn::discretize(as.data.frame(expr.dat),method = 'interval')
  # res.bnlearn <- boot.strength((expr.dat.bin),algorithm = 'hc',algorithm.args=list(score='bds',iss=iss))
  res.bnlearn <- boot.strength((expr.dat.bin),algorithm = 'hc',algorithm.args=list(
    score=score
    ,iss=iss),
  cpdag = T,R=R,...)
  # res.bnlearn <- boot.strength((expr.dat.bin),algorithm = 'hc',algorithm.args=list(score='bde'),
  #       cpdag = T,R=200)
  res.bnlearn <- arrange(res.bnlearn,desc(strength))
  # head(res.bnlearn)
  pipeline(res.bnlearn)
}
routine.bnlearn.bootstrap(expr.dat)
{
  # clu <- makeCluster(16)
  par(mfrow=c(3,3))
  routine.bnlearn.bootstrap(expr.dat,cluster=clu,score='bde')
# routine.bnlearn.bootstrap(expr.dat,cluster=clu,score='bds')
# routine.bnlearn.bootstrap(expr.dat,cluster=clu,score='bdla')
# stopCluster(clu)
}
pipeline(abs(cor(expr.dat,method='spearman')))
gres<-routine.GENIE3(expr.dat,silent=0)


pipeline((gres))

require(doParallel)
clu <- makeCluster(16)
close(clu)
doParallel::stopImplicitCluster()
# doParallel::registerDoParallel()

parallel::stopCluster(clu)
  {
  head(df.cor)
  df.cor <- arrange(df.cor,(COR))
  hasEdge <- adj.true[df.cor[,c(1,2)]%>%as.matrix]
  {
    par(mfrow=c(2,3))
    pd <- hasEdge
    pd <-rev(pd)
    # pd <- tester(param,as.pred = T)
    print(param)
    ys <- (cumsum(pd))/sum(pd)
    xs <- seq_along(pd)/length(pd)
    AROC <- abs(sum(ys-xs))*diff(xs[1:2])
    plot(xs,ys)
    title(bquote(.(AROC)))
    # plot(seq_along(),cumsum(pd)/sum(pd))
    abline(0,1,2)
  }
}

res.bnlearn <- arrange(res.bnlearn,desc(strength))
pd <- adj.true[res.bnlearn[,c(1,2)]%>%as.matrix]
# head(res.bnlearn)


{
  # ps <- df.cor
  ps <- res.bnlearn
  res <- pair2adj(ps[,c(1,2)]%>%as.matrix,is.indexed = T,fill=ps[,3],genes=genes)
  # pipeline(res)
  # pipeline(res^2)
  pipeline(res)
}

Ts  <- cumsum(pd)
Ps <- seq_along(pd)
PR <- Ts/Ps
RC <- Ts 
plot(RC,PR,type='l')
par(mfrow=c(1,1))
# plot()
# plot(Ps,Ts)
pheatmap::pheatmap(expr.dat%>%t)

library(ggplot2)

entropise <- entropy
{
  VAR <- apply(expr.dat,2,var)
  ENT <- apply(expr.dat,2,entropise)
  
  odvar = order(VAR,decreasing = T)
  od = order(ENT,decreasing = T)
  idx <-1:2
  # idx <-9:16
  print(ENT[od[idx]])
  expr.df <- as.data.frame(expr.dat[,od[idx]])
  expr.df <- reshape2::melt(expr.df,variable.name='gene',value.name='expr') 
  expr.df.bin <- reshape2::melt(expr.dat.bin[,od[idx]]%>%as.matrix,variable.name='gene',value.name='group') 
  expr.df <- cbind(expr.df, group=expr.df.bin$group)
  ggplot(expr.df) + geom_density(aes(x=expr,color=(gene),y=(..count..)) ) + xlab('expression') + 
    geom_density(aes(x=expr,y=(..count..),color=interaction(gene,group)),linetype=2) + 
    # geom_point(aes(x=expr,color=gene,y=group))
    NULL #scale_y_log10()
}
