library(Rutil)
library(dplyr)
if (interactive()){
  try(setwd(dirname(sys.frame(1)$ofile)))
}
source('routine.R')
source('subnet.R')
# source('minfo.R')


##### initialise an adjacency matrix
init.amat <- function(index,fill=F){
  m<- matrix(fill,ncol=length(index),nrow =length(index) )
  NAME <- names(index)
  colnames(m) <- NAME
  rownames(m) <- NAME
  m
}

#### Helper function to symmetrise matrixs
as.symmetric <- function(mat){
  tmat =  t(mat)
  if (is.logical(mat)){
    mat = mat | tmat
  }else if(is.numeric(mat)){
    mat = (mat + tmat)/2
  }
  mat
}
upper.tri.get <- function(mat,...){
  mat[upper.tri(mat,...)]
}

normalise <- function(x){
  (x-mean(x))/sd(x)
}
renorm <- function(x){
  x/sum(x)
}

dictify<-function(entity){
  index <- setNames(seq_along(entity),entity)
}
dict.get <- function(x,index){index[x]}

#### Backend for extracting ROC-statistics
confmat2list <- function(tb){
  if( 'TRUE' %in% rownames(tb)){
    TP = tb['TRUE','TRUE']
    FP = tb['TRUE','FALSE']
  }else{
    TP=0
    FP=0
  }
  if( 'FALSE' %in% rownames(tb)){
    TN = tb['FALSE','FALSE']
    FN = tb['FALSE','TRUE']
  }else{
    TN = 0
    FN = 0
  }
  list(
    TP = TP,
    FP = FP,
    TN = TN,
    FN = FN
  )
}

#### Take two binary vector and return ROC-statistics
confusion.matrix <- function(pred,true,as.list = T){
  res <- cbind(
    pred,
    true
  )
  colnames(res) <- c('pred','gdtruth')
  tb <- table(as.data.frame(res))
  if(as.list){
    # confmat2list(tb)
    res <- confmat2list(tb)
    res <- within(res,{PR <-TP/(TP+FP);res} )
    res <- within(res,{RC <-TP/(TP+FN);res} )
    res <- within(res,{SP <-TN/(TN+FP);res} )
    res <- within(res,{NPN<-TN/(TN+FN);res} )
    res <- within(res,{F1 <-2/(1/PR+1/RC) })
    # res$tab <- tb
    res
  }else{
    tb
  }
}

##### Take an ordered T/F prediction to calculate
##### ROC statistics at different thresholds 
confusion.matrix.from_pred <- function(tseq,as.list = T){
  TP <- cumsum(tseq)
  FP <- cumsum(1-tseq)
  P <- seq_along(tseq)
  P <- c(0,P)
  N <- rev(P)
  TP <- c(0,TP)
  FP <- c(0,FP)
  nT <- tail(TP,1)
  nF <- tail(FP,1)
  TN <- nF - FP
  data.frame(
  PR = TP/P,
  RC = TP/nT,
  SP = TN/nF,
  NPN= TN/N
  )
}




##### Convert a data.frame into an adjacency matrix
pair2adj <- function(pairs,as.list = F,genes=NULL,is.indexed=F,fill= T,symmetric = T){
  if(is.null(names(genes))) {
    index <- dictify(genes)
  }else{
    index <- genes
  }
  iadjmat <- init.amat(index)
  #### Maybe numeric avoids repeatedly hashing?
  indexer <- partial(dict.get,index=index)
  indexF <- function(x){
    x$index <- apply(x$pairs,c(2),indexer)
    x
  }
  v <- list(pairs=pairs)
  if (!is.indexed){
    v <- indexF(v)
  }else{
    v$index = v$pairs
  }
  v$adjmat <- iadjmat

  v$adjmat[v$index] <- fill
  if (symmetric){
    v$adjmat = as.symmetric(v$adjmat)
  }
  if(as.list){
    v
  }else{
    v$adjmat
  }
}

#### Take two data.frame and return the ROC-statistics
performance.pairs <- function(genes,
                              pred.pairs,
                              true.pairs){
  if(is.null(names(genes))) {
    index <- dictify(genes)
  }else{
    index <- genes
  }
  iadjmat <- init.amat(index)
  #### Maybe numeric avoids repeatedly hashing?
  indexer <- partial(dict.get,index=index)
  indexF <- function(x){
    x$index <- apply(x$mat,c(2),indexer)
    x
  }
  
  adjmats <- Map(partial(pair2adj,genes=genes),list(pred.pairs,true.pairs))

  
  #### play with the adjacency matrix
  res <- confusion.matrix(
    upper.tri.get(adjmats[[1]],diag= F)
    ,upper.tri.get(adjmats[[2]],diag= F)
    ,as.list= T
  )
  res
}




##### Benchmarking using an adjacency matrix 
##### or data.frame that specifies (from,to,score)
pipeline <- function(res,npt=100,nan.fill=0,
                     dmat.prep=identity, adj.true=NULL, genes=NULL,
                     ...){
  if(is.null(adj.true)){
    adj.true <- parent.frame()$adj.true
  }
  if(is.null(genes)){
    genes <- parent.frame()$genes
  }
  if(is.matrix(res)){
    
    res[is.nan(res)]<-nan.fill
    res<-dmat.prep(res)
    ps <- reshape2::melt(res,value.name='score')
    ps <- ps[ps[,1]!=ps[,2],]
  }else if(is.data.frame(res)){
    ps<-res
    if(!'score'%in%names(ps)){
      names(ps)[3] <- 'score'
    }
    res <- pair2adj(as.matrix(ps[,c(1,2)]),is.indexed = is.numeric(ps[,'score']),
                    fill=ps[,3],genes=genes)
  }
  
  ps <- arrange(ps,desc(score))
  
  pd <- adj.true[as.matrix(ps[,1:2])]
  qc <- confusion.matrix.from_pred(pd)
  qc.meta <- list(qc.dat=qc,adjmat=res)
  qc <- post.qc(qc.meta,dmat.prep=identity,...)
  qc
}


#### Compute AUPR and AROC after extracting ROC-statistics
post.qc<-function(qc.meta,silent=F,method='linear',...){
  qc <- qc.meta$qc.dat
  tryCatch(
    {
      funcs <- list(
        AUPR=approxfun(qc$RC,qc$PR,method=method,rule = 2,yright = 1E-4,yleft = 1E-4),
        AROC=approxfun(1-unlist(qc$SP),qc$RC,method=method,rule = 2)
      )
      #### Code looks complicated because old
      #### pipeline causes RC to not fully cover [0,1]
      funcs$F1 <- function(RC) { PR = funcs$AUPR(RC);2/(1/RC+1/PR)}
      globals = list()
      RG=rev(1-range(qc$SP))
      globals = lapply(funcs[c('AROC')], 
                       function(x) do.call(
                         integrate, c(f=x,as.list(RG),subdivisions=1000))$value)
      RG=range(qc$RC)
      globals = c(globals,
                  lapply(funcs[c('AUPR')],
                         function(x) do.call(
                           integrate, c(f=x,as.list(RG),subdivisions=1000))$value)
      )

      F1s <- 2/(1/qc$RC+1/qc$PR )
      Fidx <-which.max(F1s)
      # globals$F
      globals$F1 <- F1s[Fidx]
      coords = list(F1=c(qc$RC[Fidx], qc$PR[Fidx]) )
      
    },error=function(e){
      print("[ERR]fallback to basic plots")
      par(mfrow=c(1,2))
      within(qc,
             plot(RC,PR,type='b',ylim=c(0,1)))
      within(qc,
             plot(1-unlist(SP),RC,'b',cex=.8)
            
      )
      stop(e)
    }
  )

  qc.meta <- modifyList(qc.meta,
                        list(qc.dat=qc,globals=globals,coords=coords)
  )
  if (silent){
    # qc.meta
    ;
  }else{
    diagnose(qc.meta,...)
  }
  qc.meta
}


#### Components of Old Benchmarking pipeline
#### Now replaced with confusion.matrix.from_pred()
thresholder <- function(mat,thres){
  as.symmetric(mat >= thres)
}

select.cutoff<-function(res,N=100){
  RG = range(res)
  lst = combine_args(seq)(c(as.list(RG),length.out =N%/%2))#[-1]
  lst <- c(lst, quantile(unlist(res),seq(-0,1,length.out = N%/%2)))
  lst <- c(lst,RG[1]-.1,RG[2]+.1)
  lst <- sort(lst)
}


#### Produce diagnostic plots
diagnose <- function(qc.meta,add.dmat=T,dmat.prep=log){
  # qc <- qc.meta$qcdat
  within(qc.meta,
         {
      if(add.dmat){
        # par(mfrow = c(1,3))
        image(dmat.prep(adjmat))
        
      }else{
        # par(mfrow = c(1,2))
        ;
      }
           
      within(qc.dat,
             plot(RC,PR,type='b',ylim=c(0,1),
                  xlab='Recall',
                  ylab='Precision'))
      # lines(qc$RC,funcs$AUPR(qc$RC))
      combine_args(points)(c(as.list(coords$F1),pch=2,col=2))
      print(coords$F1)
      # text(bquote'')
      text(coords$F1[1],coords$F1[2]
           ,labels=bquote(F1==.(round(globals$F1,2) ))
           ,adj = c(-.25,-.25),pch=2,col=2) 
      title(bquote(atop(AUPR==.(globals$AUPR),
                 F1==.(globals$F1))))
      plot(1-unlist(qc.dat$SP),qc.dat$RC,'b',cex=.8
           ,ylab='True positive rate'
           ,xlab='False positive rate')
      abline(0,1,lty=2)
      title(bquote(AROC==.(globals$AROC)))
     }
  )
  
}



binarise.aggF <- function(df,aggF,MARGIN=2){
  #### Binarise a given network using some aggregation 
  apply(df,MARGIN,function(x){as.factor(x>aggF(x)) })
}


load.assignment.data <- function()
{
  with(parent.frame(),
       {
         ##### Define shared data
         datadir <-'assignment_1_files/'
         load(file.path(datadir,"medium_network_true.rda"),verbose = F)
         load(file.path(datadir,"medium_network.rda"),verbose = F)
         
         head(true.pairs)
         head(expr.dat[,1:5])
         
         genes = colnames(expr.dat)
         ngene = length(genes)
         expr.dat.bin <- as.data.frame( binarise.aggF(expr.dat,median))
         adj.true <- pair2adj(true.pairs,genes = genes)
         g.true <- igraph::graph_from_adjacency_matrix(adj.true)
       })
}

if (interactive()){
  # datadir <-'assignment_1_files/'
  # load(file.path(datadir,"medium_network_true.rda"),verbose = F)
  # load(file.path(datadir,"medium_network.rda"),verbose = F)
  load.assignment.data()
  
  {
    genes <- colnames(expr.dat)
    index <- dictify(genes)
    set.seed(0)
    N = 20
    test.pairs <- matrix(sample(genes,size = N*2,replace = F),ncol=2)
    res <-performance.pairs(genes,test.pairs,true.pairs)
    # print(res)
    res
  }
  
  routine.bnlearn.bootstrap(expr.dat,cluster=clu)
  pipeline(minet::mrnet(minet::build.mim(expr.dat)))
  # par(mfrow=c(2,3))
  # res <-routine.xgb(expr.dat)
  # pipeline(res)
  gres <- routine.GENIE3(expr.dat,nCores = 12)
  pipeline(gres)
  
  
  
  {
    A={function(){outer(index,index,function(x,y) x*0)}}
    B={function(){outer(index,index,function(x,y) rep(T,length(x)))}}
    C=function(){
      m<- matrix(F,ncol=length(index),nrow =length(index) )
      NAME <- names(index)
      colnames(m) <- NAME
      rownames(m) <- NAME
      m
    }
    microbenchmark::microbenchmark(A(),B(),C(),times=1000)
  }
  
}
