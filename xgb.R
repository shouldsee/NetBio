library(xgboost)
library(Rutil)
library(dplyr)

##### Predict gene expression (predictee) 
#### from all other genes (predictors)
##### 5-fold CV is ran to probe best parameter
####  to prevent overfitting
regGene <- function(expr.dat
  ,topred = 5
  ,nthread =1
  ,silent = 1
  ,objective ='reg:linear'
  # ,model.dir = 'model/'
  ,model.dir = NULL
  ,eta = 0.05
  ,max_depth = ncol(expr.dat)-1
  ,...
  ){
  genes<-colnames(expr.dat)
  gene<-genes[topred]
  dlabel = expr.dat[,topred]
  dtrain <- xgb.DMatrix(as.matrix(expr.dat[,-topred]),label=dlabel)

  param <- list(max_depth=max_depth,
                nthread = nthread
                ,objective=objective
                ,silent = silent
                ,eta = eta
            )
  param <- modifyList(param,list(...))
  res <- xgb.cv( param,dtrain,nrounds = 300
                 ,verbose = 1-silent
                 ,callbacks = list(cb.cv.predict(save_models = T),
                                   cb.early.stop(10,metric_name = 'test-rmse',
                                                 verbose = 1-silent)
                                   )
                 ,nfold= 5
                 )
  print(sprintf('Gene %s: \t best iteration:%d \t test_rmse_mean:%E'
                ,gene,res$best_iteration
        ,tail(res$evaluation_log$test_rmse_mean,1)))
  param$ntreelimit <- res$best_ntreelimit

  mdl = res$models[[1]]
  mdl = xgb.train(param,dtrain,nrounds=res$best_iteration)
  pd = predict(mdl,dtrain)
  xgb.attr(mdl,'gene')<-genes[topred]
  if (!is.null(model.dir)){
    fname = paste0(model.dir,genes[topred],'.xgb')
    xgb.save(fname = fname,model = mdl)
  }
  mdl
}


##### Convert a xgboost model to data.frame that 
##### outlines importance of each predictor
model2df <- function(mdl){
  gene <- xgb.attr(mdl,'gene')
  print(gene)
  Feature= setdiff(genes,gene)
  # xgb.plot.tree(mdl)
  tryCatch({
    imat <- xgb.importance(model=mdl,feature_names = Feature)
  },
  error=function(e){
    # print(e)
    if((grepl('Non-tree model detected!',e))) {
      print('[WARN] Constant tree encountered')
      imat <- data.frame(Feature,Gain=0,Cover=0,Frequency=0)
    }else{
      stop(e)
    }
  }) -> imat
  imat <- cbind(output=gene,imat)

}

### load all models in the cache folder
load.folder <- function(model.dir,callback=xgb.load,...){
  fs <- list.files(model.dir,full.names = T,...)
  print(bquote('[]Loading'~.(length(fs))~objects))
  lapply(fs, callback)
}


### The wrapper function that predict adjacency from 
### observed expression
Gxgb.fit <-function(dat,model.dir = NULL,post=T,silent=1,...){
  {
    # dat <- as.matrix(expr.dat)
    if(!is.null(model.dir)){
      dir.create(model.dir)
    }
    # browser()
    dat <- as.matrix(dat)
    genes<-colnames(dat)
    # expr.mat <- as.matrix(expr.dat)
    param <- list(
                 expr.dat=dat,
                 objective='reg:linear',
                 # objective='reg:logistic',
                 nthread=4,
                 eta=0.1,
                 max_depth= length(genes)-1,
                 alpha=0.25,lambda = 0.25,
                 model.dir = model.dir,
                 silent= silent
    )
    param <- modifyList(param,list(...))
    param$f = regGene
    f =  combine_args(partial)(param)
    mdls = lapply( seq_along(genes), f)
    if(post){
      df<-Gxgb.post(dat,mdls)
    }else{
      mdls
    }
  }
}


#### Post processing of the raw importance scores
Gxgb.post <- function(dat,mdls=NULL,model.dir=NULL){
  if (is.null(mdls)){
    mdls<-load.folder(model.dir,callback = xgb.load)
  }
  df = rbind_list(lapply(mdls,model2df))
  VAR = apply(dat,2,var)
  df$score <-  df$Gain * VAR[df$output] * VAR[df$Feature]
  df
}





##### Example  
# if (interactive()){
#   {
#     load('assignment_1_files/medium_network.rda',verbose = T)
#     pc <- prcomp(t(expr.dat))
#     par(mfrow=c(1,2))
#     biplot(pc)
#     plot(apply(expr.dat,2,var)%>%sort)
#     
#     # expr.dat <- as.matrix(expr.dat)
#     expr.dat.std <- apply(expr.dat,2,function(x)(x-mean(x))/sd(x) )
#     colnames(expr.dat.std) <- colnames(expr.dat)
#     # df <- Gxgb.fit(expr.dat.std)
#     df <- Gxgb.fit(expr.dat)
#     # df <- Gxgb.post(expr.mat,model.dir = 'model/')
#   }
#   
#   load('assignment_1_files/medium_network_true.rda')
#   genes <- colnames(expr.dat)
#   adj.true <- pair2adj(true.pairs,genes = genes)
#   {
#     pairs = as.matrix(df[,c('Feature','output')])
#     res <- pair2adj(pairs,genes=genes,is.indexed = F,symmetric = F
#                     ,fill = df$score
#                     # ,fill = df$Gain
#     )
#     qc <- pipeline(res)
#   }
# }

