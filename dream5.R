library(Rutil)
library(dplyr)
if (interactive()){
  try(setwd(dirname(sys.frame(1)$ofile)))
}


# true.pairs


init.amat <- function(index,fill=F){
  m<- matrix(fill,ncol=length(index),nrow =length(index) )
  NAME <- names(index)
  colnames(m) <- NAME
  rownames(m) <- NAME
  m
}

as.symmetric <- function(mat){
  tmat =  t(mat)
  if (is.logical(mat)){
    mat = mat | tmat
  }else if(is.numeric(mat)){
    mat = (mat + tmat)/2
  }
  mat
}


dictify<-function(entity){
  index <- setNames(seq_along(entity),entity)
  # partial(dict.get,index=index)
}
dict.get <- function(x,index){index[x]}

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
df <- confusion.matrix.from_pred(pd)
plot(df$RC,df$PR)
plot(1-df$SP,df$RC)

upper.tri.get <- function(mat,...){
  mat[upper.tri(mat,...)]
}


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
  # lidx <- sub2ind(v$index[,1],v$index[,2],mat=iadjmat)
  # v$adjmat[lidx] <- T
  
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
  # pred <- list(mat=pred.pairs)
  # gdtruth<- list(mat=true.pairs)
  # for (vname in c('pred','gdtruth')){
  #   v <- indexF(get(vname))
  #   v$adjmat <- iadjmat
  #   # lidx <- sub2ind(v$index[,1],v$index[,2],mat=iadjmat)
  #   # v$adjmat[lidx] <- T
  #   v$adjmat[v$index] <- T
  #   v$adjmat <- t(v$adjmat)|v$adjmat
  #   assign(vname,v)
  # }
  
  
  
  #### play with the adjacency matrix
  res <- confusion.matrix(
    upper.tri.get(adjmats[[1]],diag= F)
    ,upper.tri.get(adjmats[[2]],diag= F)
    ,as.list= T
  )
  # res <- confmat2list(confmat)
  # res <- within(res,{PR <-TP/(TP+FP);res} )
  # res <- within(res,{RC <-TP/(TP+FN);res} )
  # res <- within(res,{SP <-TN/(TN+FP);res} )
  # res <- within(res,{NPN<-TN/(TN+FN);res} )
  # res <- within(res,{F1 <-2/(1/PR+1/RC) })
  # res$tab <- confmat
  res
}





##### Benchmarking pipeline
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


post.qc<-function(qc.meta,silent=F,method='linear',...){
  qc <- qc.meta$qc.dat
  # qc[,-ncol(qc)]<-unlist.df(qc[,-ncol(qc)])
  
  # qc= qc[complete.cases(qc[,-ncol(qc)]),]
  tryCatch(
    {
      funcs <- list(
        AUPR=approxfun(qc$RC,qc$PR,method=method,rule = 2,yright = 1E-4,yleft = 1E-4),
        AROC=approxfun(1-unlist(qc$SP),qc$RC,method=method,rule = 2)
        # AROC=approxfun(qc$RC,1-unlist(qc$SP),)
      )
      funcs$F1 <- function(RC) { PR = funcs$AUPR(RC);2/(1/RC+1/PR)}
      # F1 = )
      # AUPR <- integrate(,0,1)
      # AROC <- integrate(approxfun(1-unlist(qc$SP),qc$RC),0,1)
      globals = list()
      # within(globals,AUPR <- sum(convolve(qc$PR*qc$) 
      RG=rev(1-range(qc$SP))
      globals = lapply(funcs[c('AROC')], function(x) do.call(integrate, c(f=x,as.list(RG),subdivisions=1000))$value)
      RG=range(qc$RC)
      globals = c(globals,
                  lapply(funcs[c('AUPR')], function(x) do.call(integrate, c(f=x,as.list(RG),subdivisions=1000))$value)
      )
      # opres <- optim(c(0.5),Fcompose(funcs$F1, function(x){-x} ),method='BFGS')
      # globals$F1= -opres$value
      # coords = list(F1=c(opres$par,funcs$AUPR(opres$par) ) )
      
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


diagnose <- function(qc.meta,add.dmat=T,dmat.prep=log){
  # qc <- qc.meta$qcdat
  within(qc.meta,
         {
      if(add.dmat){
        par(mfrow = c(1,3))
        image(dmat.prep(adjmat))
        
      }else{
        par(mfrow = c(1,2))
      }
           
      within(qc.dat,
             plot(RC,PR,type='b',ylim=c(0,1)))
      # lines(qc$RC,funcs$AUPR(qc$RC))
      combine_args(points)(c(as.list(coords$F1),pch=2,col=2))
      print(coords$F1)
      # text(bquote'')
      text(coords$F1[1],coords$F1[2],labels=bquote(F1==.(round(globals$F1,2) )),adj = c(-.25,-.25),pch=2,col=2) 
      title(bquote(atop(AUPR==.(globals$AUPR),
                 F1==.(globals$F1))))
      plot(1-unlist(qc.dat$SP),qc.dat$RC,'b',cex=.8)
      abline(0,1,lty=2)
      title(bquote(AROC==.(globals$AROC)))
     }
  )
  
}


pipeline <- function(res,npt=100,nan.fill=0,dmat.prep=identity,...){
  if(is.matrix(res)){
    
  res[is.nan(res)]<-nan.fill
  res<-dmat.prep(res)
  # thres.lst <- select.cutoff(res,npt)
  # qc <- lapply(thres.lst,function(thres) {
  #   mat <- thresholder(res,thres)
  #   confusion.matrix( upper.tri.get(mat,diag=F),
  #                     upper.tri.get(adj.true,diag=F))} )
  # qc <- rbind_list(qc) 
  # qc <- cbind(thres=thres.lst,qc) %>% as.data.frame 
  
  ps <- reshape2::melt(res,value.name='score')
  }else if(is.data.frame(res)){
    # stop('Not implemented for data.frame')
    ps<-res
    if(!'score'%in%names(ps)){
      names(ps)[3] <- 'score'
    }
    
    res <- pair2adj(as.matrix(ps[,c(1,2)]),is.indexed = is.numeric(ps[,'score']),
                    fill=ps[,3],genes=genes)
    # ps <- res
  
  }
  
  ps <- arrange(ps,desc(score))
  
  pd <- adj.true[as.matrix(ps[,1:2])]
  qc <- confusion.matrix.from_pred(pd)
  qc.meta <- list(qc.dat=qc,adjmat=res)
  
  qc <- post.qc(qc.meta,dmat.prep=identity,...)
  # qc$adjmat <- res
  qc
}
# pipeline(res)


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

randomWalk <- function(links,size= 50,directed=F)
{
  ### A heuristic random walk on a graph
  if (is.matrix(links)){
    nodes <- colnames(links)
    amat <- links
  }else{
    nodes <- unique(unlist(links))
  }
  init <- sample(which(rowSums(amat)!=0,1),1)
  # cur<-nodes[init
  idxseq <- sample(seq_along(nodes),size=size-1,replace = T)
  idxseq <- c(idxseq)
  outseq <- c(init,idxseq)
  idx <- init
  for (i in seq_along(idxseq)){
    outseq[i+1] <- idx
    if(directed){
      pointer = amat[idx,] 
    }else{
      pointer = amat[idx,] | amat[,idx]
    }
    idx <- which(pointer)[ idxseq[i]%%sum(pointer)+1 ]
  }
  nodes[unique(outseq)]
}




#' Estimate entropy by KDE-fitted distribution 
#' @export
entropise <- function(x,eps=1E-3){
  ds <- density(x)
  dx <- diff(ds$x[1:2])
  lp <- ds$y*dx
  ent <- sum(lp*-log(lp+eps))
  # eps =
  # ent = -sum(ds$y*log(ds$y+eps))*dx
  # ent <- do.call(integrate,c(f=approxfun(ds$x,ds$y*log(ds$y+eps)),as.list(range(ds$x))))
}


#' Estimate total entropy of 2 r.v. by 2d-KDE
#' @export
entropy2d<-function(x,y,nbin=c(25,25),eps=1E-3){
  ds <- MASS::kde2d(x,y,n = nbin)
  # eps = 0
  dx <- diff(ds$x[1:2])
  dy <- diff(ds$y[1:2])
  lp <- ds$z*dx*dy
  -sum(lp*log(lp+eps))
}

wrong.entropy2d<-function(x,y,nbin=c(25,25),eps=1E-3){
  ds <- MASS::kde2d(x,y,n = nbin)
  # eps = 0
  dx <- diff(ds$x[1:2])
  dy <- diff(ds$y[1:2])
  ent <- -sum(ds$z*(log(ds$z+eps)))*dx*dy
}


#' Make pairwise MI matrix from a dataset
#' @export
make.mi <- function(expr.dat){
  d.e2d <- pdist(expr.dat,partial(entropy2d,nbin=c(30,30)))
  # d.e2d <- proxy::dist(expr.dat
  #                      ,method=partial(entropy2d,nbin=c(30,30))
  #                      ,diag = T,by_rows = F) %>%as.matrix
  expr.h <- apply(expr.dat,2,entropise)
  d.h <- outer(expr.h,expr.h,'+')
  d.mi = d.h - d.e2d 
}

#' Make pairwise distance matrix
#' @export
pdist<-function(X,distF){
  D <- proxy::dist(X,method =distF,by_rows = F,diag = T)%>%as.matrix
  diag(D)<-apply(X,2,function(x)distF(x,x))
  D
}


subdata <- function(genes){
  with(.GlobalEnv,
       {adj.true <- adj.true[genes,genes]
       expr.dat<-expr.dat[,genes]
       g.true <- subgraph(g.true,genes)
       })
}

routine.GENIE3<-function(dat,nCores=6,silent=1,...){
  require(GENIE3)
  require(doParallel)
  qc <- pipeline(GENIE3(t(dat),nCores=nCores,...),silent=silent)
  qc$method =  'GENIE3'
  .GlobalEnv$qc.list$GENIE3 <- qc
  qc$adjmat
}

normalise <- function(x){
  (x-mean(x))/sd(x)
}
renorm <- function(x){
  x/sum(x)
}


if (interactive()){
  datadir <-'assignment_1_files/'
  load(file.path(datadir,"medium_network_true.rda"),verbose = F)
  load(file.path(datadir,"medium_network.rda"),verbose = F)
  
  {
    genes <- colnames(expr.dat)
    index <- dictify(genes)
    set.seed(0)
    N = 20
    test.pairs <- matrix(sample(genes,size = N*2,replace = F),ncol=2)
    res <-perfromance.pairs(genes,test.pairs,true.pairs)
    # print(res)
    res
  }
  
  
  
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
