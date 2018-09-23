nodes <- 1:3

Ln <- length(nodes)
# unlist(apply(combn(nodes,2),2,list),recursive = F)
unlist.once <- function(x,...){
  unlist(x,recursive = F,...)
}

nloop=2 ### number of nodes in the loop

rbind.list <- Rutil::combine_args(rbind)
make_itF<- function(nodes)
{ 
  itF<-function(nloop){ 
   arr = apply(combn(nodes,nloop),2,function(x) list(nodes=x,nloop=nloop))
   arr = rbind.list(arr)
  }
}


score.loop <- function(tb,loops,...){
  sapply( loops[,1],function(nodes){lp.nodes(nodes,tb=tb,...)})
}

score.wt <- function(wt,lpscore,nloop,check.valid = T){
  Ln = round(log2(length(nloop) + 1))
  if(check.valid){
    stopifnot( sum(nloop * wt ) == Ln)
  }
  sum(lpscore * wt)
}


make_tripgrid <- function(n){
  wts <- compose(combine_args(expand.grid),array2array)(rep(list(-1:1),n))
}
# lp.nodes()
aM=tail(all.graphs,1)[[1]]
amat2bn(aM)

plot(amat2bn(aM))


a = make_itF(nodes)(1)

loops <- rbind.list(lapply(1:Ln, make_itF(nodes)) )
nloop <-rbind.list(loops[,2])
loops[,1]

stopifnot(sum(nloop * wt ) == Ln)
wt <- c(1,1,0,-1,0,0,1)
loops[wt!=0,1]

# sess$dat=dat1
dat <- dat1
make_loops<-function(Ln){
  nodes = 1:Ln
  loops <- rbind.list(lapply(nodes, make_itF(nodes)) )
}
prepare_scorer<-function(dat,loops = NULL){
  tb <- table(dat)
  Ln = ncol(dat)
  if (is.null(loops)){
    loops<-make_loops(Ln)
  }
  nloop <-rbind.list(loops[,2])
  lpscore = score.loop(tb,loops,equiv=F)
  scorer <- partial(score.wt,lpscore=lpscore,nloop=nloop,check.valid=F)
}

Ln = 3
wts.valid <- seqF(Ln,return.wt = T)
loops <- make_loops(Ln)
{
  lst = matrix(0,ncol=28)
  # lst = matrix(ncol=28)
  
  for( i in 1:50) { 
    N= 300
    dat = matrix(runif(N*3)>0.5,ncol=3)  %>% as.data.frame
    # idx = sample(1:N,size=round(N/2),replace =F)
    # dat[idx,2] = (dat[idx,1]|dat[idx,3]) |  dat[idx,2]
    # idx = sample(1:N,size=round(N/2),replace =F)
    # dat[idx,3] = (dat[idx,1]|dat[idx,2]) |  dat[idx,3]
    # dat[idx,3] = dat[idx,2]
    # idx = sample(1:N,size=round(N/2),replace =F)
    # dat[idx,1] = dat[idx,3]
    # dat =dat2
    sc <- prepare_scorer(dat,loops)
    lls = apply(wts.valid,1,sc)
    lst = rbind(lst,lls)
  }
  lst = lst[-1,]
}

{
  par(mfcol=c(1,2))
  plot(MEAN<-apply(lst,2,mean))
  par(new=T)
  plot(apply(lst,2,sd),type='l',axes=F)
  axis(4)
  plot(lst[,16],lst[,28])
  abline(1,1)
}


bests <- order(MEAN,decreasing = T)[1:5]
cbind(wts.valid,MEAN)[bests,]
y = -wts.valid[28,]
y[7]=1
which(apply(wts.valid,1,function(x){all(x==y)}))
wts.valid
{
  bestgi <- which.max(MEAN)
  print(bestgi)
  print(lls[bestgi])
  wts.valid[bestgi,]
}
wts.valid
loops[,1]
# sc(aM)
# score.wt(wt,lpscore,nloop )
myalgo(aM,sess$dat,equiv=F,eta0=1)
bnlearn::score(amat2bn(aM),sess$dat,type='bde',iss=8)


wts.valid <- seqF(nrow(aM),return.wt = T)
engs <- apply(wts.valid,1,scorer)
max(engs)
wt = wts.valid[which.max(engs),]
print(loops[wt!=0,1])
print(wt)



seqF <- function(Ln,return.wt = F){
  itF = make_itF(1:Ln)
  loops <- rbind.list(lapply(1:Ln, itF))
  nloop <-rbind.list(loops[,2])
  D = diag(Ln)
  loop.mat <- rbind.list(lapply(loops[,1], function(x) colSums(D[x,,drop=F])))
  
  all(rowSums(loop.mat)==nloop)
  
  wts <- make_tripgrid(length(nloop))
  card.node <- apply(bsxfun(expand_dims(loop.mat,1),wts,'*'),c(1,3),sum)
  idx.valid = apply(card.node==1,1,all)
  #idx.valid = card==Ln
  # colSums()
  # card <- rowSums(bsxfun(array(nloop,dim=c(1,length(nloop))),wts,'*'))
  # all(rowSums(card.node)==card)
  # hist(card)
  
  wts.valid <- wts[idx.valid,]
  if (return.wt){
    return(wts.valid)
  }
  sum(idx.valid)
}
Map(seqF,1:3)

{
  Ln = 3
  itF = make_itF(1:Ln)
  loops <- rbind.list(lapply(1:Ln, itF))
  nloop <-rbind.list(loops[,2])
  D = diag(Ln)
  loop.mat <- rbind.list(lapply(loops[,1], function(x) colSums(D[x,,drop=F])))
  
  all(rowSums(loop.mat)==nloop)
  
  wts <- make_tripgrid(length(nloop))
  card.node <- apply(bsxfun(expand_dims(loop.mat,1),wts,'*'),c(1,3),sum)
  idx.valid = apply(card.node==1,1,all)
  #idx.valid = card==Ln
  # colSums()
  # card <- rowSums(bsxfun(array(nloop,dim=c(1,length(nloop))),wts,'*'))
  # all(rowSums(card.node)==card)
  # hist(card)
  
  wts.valid <- wts[idx.valid,]
}

# bnlearn::hist
wts.valid
sum(idx.valid)

