# with(.GlobalEnv,{
with(env.sub,{
    idx = which(adj.true,arr.ind = T)
  N = 36
  par(mfrow=c(sqrt(N),sqrt(N)),
      mar=c(0,0,2,0))
  # N = nrow(idx)
  for (i in 1:N ){
    dat <- expr.dat[,idx[i,]]
    ii <- idx[i,1]
    j <- idx[i,2]
    plot(dat)
    title(paste(ii,j,round(cor(dat[,1],dat[,2]),3),sep=','))
  }
}
)

which(adj.true[,1])
with(env.sub,{
  idx = which(adj.true,arr.ind = T)
  k = 5
  i =idx[k,2]
  j= idx[k,1]
  print(which(adj.true[,i]))
  print(which(adj.true[,j]))
  print(c(i,j))
  adj.true[5,27]
  par(mfrow=c(sqrt(N),sqrt(N)),
  mar=c(0,0,1,0))
  COL <- bnlearn::discretize(data.frame(expr.dat[,i,drop=F]))
  COL <- as.numeric(unlist(COL))
  for (ii in 1:nrow(adj.true) ){
  dat <- expr.dat[,c(ii,j)]
  plot(dat,col=COL)
  title(paste(ii,j,round(cor(dat[,1],dat[,2]),3),sep=','))
  }
  
  print(c(i,j))
}
)

env.sub <- subnet(env.data,ngene = 30)
runlm <- function(e){
with(e,{
# load.assignment.data()
# with(.GlobalEnv,{
  ##### The fitting using linear model is robust 
  #### across the observations
  
  # res <- lm(G3812~G333+G313+G268,data=data.frame(expr.dat))
  expr.df <- as.data.frame(expr.dat)
  gene <- genes[1]
  gene <- 'G3812'
  iter <- function(gene){
  form <- as.formula(sprintf('%s~.',gene)  )
  cvi <- cvFolds(nrow(expr.dat))
  res <- lm(form,data=expr.df[cvi$subsets[cvi$which!=1], ])
  # res <- cvFit(lmrob,form,data=expr.df)
  # res <- glm(form,data=data.frame(expr.dat))
  # res <- with(expr.df,glmnet(form,data))
  # plot(res$fitted.values,res$model$G3812)

  # plot(res$fitted.values,res$model[[gene]])
  test.data <- expr.df[cvi$subsets[cvi$which==1], ]
  pd <- cbind(test.data[[gene]],predict(res,test.data))
  plot(pd)
  title( paste(round(mean(abs(pd[,2]-pd[,1])),3),
               round(mean(abs(res$residuals)),3)
               ))
  # anova(res)
  # summary(res)$r.squared
  res
  }
  lapply(genes,iter)
  # iter(gene)
}) ->res #-> rsq
  rsq <- sapply(res,function(x){summary(x)$r.squared})
  summary(rsq)
  plot(density(rsq),xlim=c(0,1))
  res
}
N = (sqrt(length(genes))+1)^2
par(mfrow=c(sqrt(N),sqrt(N)),
    mar=c(0,0,1,0))
res <- runlm(subnet(env.data,ngene = 40))
res <- runlm(subnet(env.data,ngene = 40))

res <- runlm(env.sub)
res <- runlm(.GlobalEnv)


cvi <-cvFolds(100,K = 5)
# cvi$subsets
# data.frame(cvi)
lmrob
# ?cvFit
# a= cvSelect(res[[1]])

# ?cvFit
library(cvTools)
getcv <- function(i){
  form <- as.formula(sprintf('%s~.',parent.frame(2)$genes[i])  )
  cvLm(res[[i]])
}
# with(env.sub,lapply(seq_along(genes),getcv)) ->cvs
# cvs

# summary(o)
str(o)

image(ce <- rbind_list(lapply(res,coef))[,-1]%>%as.matrix%>%abs)
ce <- with(env.sub,lapply(seq_along(genes),function(i){
  x = coef(res[[i]])
  names(x)[1] <- genes[i]
  x[1] <- 0
  x
}))%>%rbind_list%>%as.matrix
dimnames(ce) <-NULL
# par(mf)
with(env.sub,
     {
       C <- cor(expr.dat)
       pipeline(abs(C))
     }
)
with(env.sub,pipeline(abs(ce)))
par(mfrow=c(3,3),mar=c(5,3,3,2))

with(env.sub,dim(adj.true))
# par(mfrow=c(2,2))
od <- order(a$coefficients[,4],decreasing = T)
a$cov.unscaled
a <- anova(res)
od <- order(anova(res)$F,decreasing = T)
with(env.sub,ggpairs(as.data.frame(expr.dat[,c(5,od[1:5])])))

library(glmnet)
summary(rsq)




