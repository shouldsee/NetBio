source('fitting_beta.R')
source('bayes_infer.R')
library(Rutil)

stratifier <- function(.self,selfNode,eta0=1,
                       eta = array(eta0,dim(.self$tb)), #### potentially caching eta
                       parent_state=NULL,
                       ...){
  sess <- .self
  parNode = as.numeric(sess$F$get_parent(selfNode))
  #### ignore parent_state if parNode is not matched
  parent_state <- parent_state[seq_along(parNode)]
  #### ignore parNode if parent is not specified
  parNode <- parNode[seq_along(parent_state)]
  nodes <- c(parNode,selfNode)
  idx = cbind(rbind(parent_state,parent_state),child_state=1:2)
  list(
    tb.marg = margin.table(sess$tb,margin = nodes)[idx],
    eta.marg = margin.table(eta,margin = nodes)[idx])
}


logL_maker <- function(obsB,binom=T){
  #### Make a likelihood function based on a observation vector
  if(length(obsB)==2){
    N=sum(obsB)
    X=obsB[2]
  }else{
    N = length(obsB)
    X = sum(obsB)
  }
  if(binom){
    function(p1,...) {
      dbinom(x= X,size = N,p1,...)}
  }else{
    function(p1,...){
      p1^X*(1-p1)^(N-X)
    }
  }
}

net_posterior <- function(.self,
                          childNode,
                          parent_state=rep(1,length(.self$F$get_parent(childNode))) ,
                          eta0=1
                          ,debug=0
                          ,binom= F
                          ){
  sess <- .self
  out <- stratifier(sess,
                    childNode,
                    parent_state = parent_state)
  
  if(debug){print(out)}
  eta<- out$eta.marg
  
  ob <- .infer$new()
  ob$Fparam$prior <- partial(dbeta,shape1=eta[1],shape2=eta[2])
  L <-logL_maker(out$tb.marg,binom=binom)
  ob$Fparam$L <- L
  ob$Fparam$logL <- partial(L,log=T)

  {
    make_dpos(ob)
    marg_obs(ob)
    make_post(ob)
  }
  return(ob)
}


#### Examples
# 
# aM <- all.graphs[[3]]
# sess <- .PGM_binary$new(mdlgraph=igraph::graph_from_adjacency_matrix(aM))
# sess$preprocess()
# 
# {
#   sess$dat <-dat1
#   sess$make_table()
#   plot(sess$mdlgraph)
#   print(sess$tb)
#   ob <- sess$net_posterior(2,parent_state = 1,debug=1)
#   plot_bayes(ob)
#   title(bquote(beta(4,4)))
# }
# 
# {
#   sess$dat <-dat2
#   sess$make_table()
#   print(sess$tb)
#   ob <- sess$net_posterior(2,parent_state = c(1),debug=1)
#   plot_bayes(ob)
#   title(bquote(beta(4,4)))
# }
