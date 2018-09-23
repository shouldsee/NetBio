
routine.GENIE3<-function(dat,nCores=6,silent=1,...){
  require(GENIE3)
  require(doParallel)
  GENIE3(t(dat),nCores=nCores)
}

routine.bnlearn.bootstrap <- function(expr.dat,iss=1,R=200,score='bde',silent=0,...){
  require(bnlearn)
  
  expr.dat.bin <- bnlearn::discretize(as.data.frame(expr.dat),method = 'quantile',breaks=3)
  res.bnlearn <- boot.strength((expr.dat.bin),algorithm = 'hc',algorithm.args=list(
    score=score
    ,iss=iss),
    cpdag = T,R=R,...)
  res.bnlearn <- arrange(res.bnlearn,desc(strength))
}

routine.xgb <- function(expr.dat,...){
  
  source('xgb.R')
  df <- Gxgb.fit(expr.dat,model.dir = 'qc/',silent=1)
  pairs = as.matrix(df[,c('Feature','output')])
  VAR <- apply(expr.dat,2,var)
  res <- pair2adj(pairs,genes=genes,is.indexed = F,symmetric = F,
                  fill = df$Gain*VAR[df$output]*VAR[df$Feature])
  qc <- pipeline(res,...)
  qc$method <-'xgb'
}

routine.mrnet <- function(expr.dat){
  minet::mrnet(minet::build.mim(expr.dat))
}
# 
# if(interactive()){
#   PKGMethod='MI.inverted'
#   fig.cap=PKGMethod
#   PKG = strsplit(PKGMethod,'\\.')[[1]][1]
#   # library(PKG,character.only = T)
#   res <- mat.invert(make.mi(expr.dat),post=abs)
#   qc <- pipeline(res,silent=1)
#   qc$method =  PKGMethod
#   qc.list[[PKGMethod]] <- qc
#   
#   PKGMethod='MI.raw'
#   fig.cap=PKGMethod
#   PKG = strsplit(PKGMethod,'\\.')[[1]][1]
#   # library(PKG,character.only = T)
#   res <- (make.mi(expr.dat))
#   qc <- pipeline(res,silent=1)
#   qc$method =  PKGMethod
#   qc.list[[PKGMethod]] <- qc
#   
#   
#   PKGMethod='COV.inverted'
#   fig.cap=PKGMethod
#   PKG = strsplit(PKGMethod,'\\.')[[1]][1]
#   # library(PKG,character.only = T)
#   res <- mat.invert(cov(expr.dat),post = abs)
#   qc <- pipeline(res,silent=1)
#   qc$method =  PKGMethod
#   qc.list[[PKGMethod]] <- qc
#   
# }
# 

