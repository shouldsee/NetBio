Rutil::install.packages.lazy(c('dagitty','pcalg','jpeg'))
library(igraph)

require(igraph)
require(Rutil)
source('dirichlet.R')

datadir = './assignment_1_files/'
load(file.path(datadir,"small_network_1.rda"),verbose = 1)
load(file.path(datadir,"small_network_2.rda"),verbose = 1)


#' @export
lp.diri <- function(n,eta,equiv = T){
  Sn = sum(n)
  Se = sum(eta)
  logP = lgamma(Se) - lgamma(Se + Sn) + sum(lgamma(n+eta) - lgamma(eta))
  if (equiv){
    ### apply equivalency correction (optionial)
    logP = logP + ndchoose(n,log = T)
  }
  logP
}


lp.nodes <- function(nodes,tb,eta=array(1,dim(tb)),...){
  tb.marg = margin.table(tb,margin = nodes)
  eta.marg = margin.table(eta,margin = nodes)
  lp = lp.diri(tb.marg,eta.marg,...)
}
source('net_post.R')
.PGM_binary <- setRefClass('PGM_binary',
                      fields = list(
                        # dat=c('array','data.frame'),
                        dat='data.frame',
                        tb ='table',
                        mdlgraph='ANY',
                        F = 'list'
                      ),methods = list(
                        make_table = function(){tb <<- table(dat)},
                        preprocess = function(){
                          g <- mdlgraph
                          pars = adjacent_vertices(g, V(g), mode = c("in"))
                          F$get_parent <<-  function(x)pars[[x]]
                        },
                        lp.node= function(.self,selfNode,eta0=1,
                                           eta= array(eta0,dim(.self$tb)),...){
                          sess = .self
                          parNode = sess$F$get_parent(selfNode)
                          nodes = c(selfNode,parNode)
                          lp.nodes(nodes,sess$tb,eta,...) - lp.nodes(parNode,sess$tb,eta,...)
                        },
                        logL = function(.self,...){
                          sess = .self
                          lps = sapply( V(sess$mdlgraph), partial(sess$lp.node,...))
                          
                          sum(lps)
                        },
                        net_posterior = net_posterior
                      )
                    )


##### Testing the algorithm
source('graph_2edge.R')
sess = .PGM_binary$new(
  # dat=dat1,
  dat=dat2,
  mdlgraph=graph_from_adjacency_matrix(aM))
sess$make_table()
sess$preprocess()
sess$lp.node(1)%>%print
sess$logL()%>%print



##### My custom algorithm that computes the 
##### likelihood of a dag on a dataset
myalgo <- function(g,dat,...){
  if (is(g,'matrix')){
    g <- graph_from_adjacency_matrix(g)
  }
  sess = .PGM_binary$new(
    dat=dat,
    mdlgraph=g)
  sess$make_table()
  sess$preprocess()
  sess$logL(...)
}

source('graph_2edge.R')

library(bnlearn)
#' Convert an adjacency matrix to a "bnlearn" object
#' @export
amat2bn <- function(aM,nodes = rownames(aM)){
  # require(bn)
  if (is.null(nodes)){
    nodes = LETTERS[1:nrow(aM)]
  }
  g0 <- bnlearn::empty.graph(nodes)
  bnlearn::amat(g0) <- aM
  g0
}


###### Load all graphs (Stored as adjacency matrix)
lst = list()
data.scripts = c('graph_0edge.R',
                 'graph_1edge.R',
                 'graph_2edge.R'
)
for (f in data.scripts){
  source(f)
  lst = c(lst,aMs)
}
all.graphs <- lst



#### Helper functions to compute data.frames

aM2bn <- function(aM,dat){
  eta0 = 1
  bn <- amat2bn(aM)
  bn$iss = 2^nrow(aM) * eta0
  bn$dat <- dat
  bn
}


{bn2df <- function(bn, iss=bn$iss){
  bnscore <- score(bn,bn$dat,'bde',iss=iss)
  bnscore <- score(bn,bn$dat,'bde',iss=iss)
  # bn
  myalgo.bde<- myalgo( amat(bn),bn$dat,eta0 = iss/8,equiv=F)
  list(model=modelstring(bn),myalgo.bde.iss=myalgo.bde,bnlearn.bde.iss=bnscore
       ,bnlearn.bic=score(bn,bn$dat,'bic')
       ,iss=iss
  )
}}


main <- function(dat = dat1,...){
  # bn.list <- lapply(lst,aM2bn)
  f = partial(aM2bn,dat=dat)
  bn.list <- lapply(lst,f)
  # bn.list <- mapply(aM2bn,lst,dat)
  g <- partial(bn2df,...)
  # df <- sapply(bn.list,g)
  df <- Rutil::combine_args(rbind)(lapply(bn.list,g))
  df <- Rutil::unlist.df(data.frame(df))
  df <- cbind(df,dat= deparse(substitute(dat)))
}


