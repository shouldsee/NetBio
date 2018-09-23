

#### Generate a subnet 
#### from a list, return a list
subnet <- function(env.data,ntrial=10,ngene=10,...){
  if(is.null(env.data$gmat)){
    env.data$gmat = gmat_from_adjacency(env.data$adj.true)
  }
  with(env.data,
       {
         for (i in 1:ntrial){
           gis <- unique(markov(gmat,...))
           if(length(gis)>ngene){gis=head(gis,ngene);break}
         }
         
         list(
           genes = genes[gis],
           adj.true=adj.true[gis,gis],
           expr.dat = expr.dat[,gis],
           g.true = graph_from_adjacency_matrix(adj.true[gis,gis])
         )
       }
  )
}


#### Markov chain sampler
#### Perform random walk on a graph
gmat_from_adjacency <- function(amat){
  if (is(amat,'igraph')){
    amat <- as.symmetric(as.matrix(igraph::as_adjacency_matrix(amat)))
  }
  pmat <- apply(amat,2,renorm)
  generator.mat <- apply(pmat,2,cumsum)
}
lookup<-function(r,cur,gen.mat){
  which(r<gen.mat[,cur])[1]
}
markov <- function(generator.mat,init.prob=rep(1,nrow(generator.mat)),nT=100){
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



#### Runing bnlearn and GENIE3 side by side
contraster <- function(){
  env.sub <- subnet(env.data,nT=50)
  if(length(env.sub$genes)<5){return(list())}
  with(env.sub,
       {
         print(length(genes))
         t1 <- system.time({res <- routine.bnlearn.bootstrap(expr.dat,clu=clu)})
         r1 <- pipeline(res,adj.true = adj.true,silent=1)
         t2 <- system.time({res <- routine.GENIE3(expr.dat,nCores=16)})
         r2 <- pipeline(res,adj.true = adj.true,silent=1)
         r1$globals$runtime <- t1[3]
         r1$globals$method  <- 'bnlearn'
         r2$globals$runtime <- t2[3]
         r2$globals$method  <- 'GENIE3'
         
         d <-rbind(r1$globals,r2$globals)
         d <-cbind(d,ngene =length(genes))
         d
       }) ->out 
} 



plot.lst<-function(lst){
  require(ggplot2)
  df <- as.data.frame(lst)
  df <- unlist.df(df)
  df$i <- ceiling(1:nrow(df)/length(unique(df$method)))
  
  p1 <- ggplot(df) + geom_histogram(aes(x=AUPR,fill=method,y=..count..),position = "dodge",alpha=0.25) + 
    geom_density(aes(x=AUPR,color=method,y=..density..)) + xlim(0,1) + theme(legend.position = 'top')
  
  # ggplot(df) + geom_density(aes(x=AUPR,color=method,y=..count..))
  df.melt <- reshape2::melt(df,id.vars=c('i','method'),measure.vars=c('AUPR'),
                            value.name='AUPR')
  val1 <- filter(df,method=='bnlearn')
  val2 <- filter(df,method=='GENIE3')
  p2 <- ggplot()+geom_point(data=data.frame(bnlearn.AUPR=val1$AUPR,GENIE3.AUPR=val2$AUPR),
                            aes(x=bnlearn.AUPR,y=GENIE3.AUPR)) +
    geom_abline(intercept = 0,slope=1) +xlim(0,1)+ylim(0,1) + 
    ggtitle(bquote(atop(
      bnlearn.avg_runtime==.(mean(val1$runtime))~s,
      GENIE3.avg_runtime==.(mean(val2$runtime))~s
    )
    ))
  p <- gridExtra::grid.arrange(p1,p2,ncol=2)
}



### (Deprecated) Old heuristic random walk on a graph
### non-markovian
randomWalk <- function(links,size= 50,directed=F)
{
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
