Rutil::install.packages.lazy(c('dagitty','pcalg','jpeg'))
library(igraph)
# 
# library(dagitty)
# library(pcalg)
# x <- dagitty('dag{
#   A -- B -- C
#              }')
# # attributes(x)
# # x@
# # x$dagitty
# # is(x)
# # classso
# x = graphLayout(x)
# plot(graphLayout(x))
# # plot(x)
# # equivalenceClass(x)
# eqs = equivalentDAGs(x)
# lapply( eqs,function(x)plot(graphLayout(x)))
# lapply( eqs,function(x)plot(graphLayout(x)))
# # ?equiv
# # dagitty::equ
# eqs[[1]]
# # ?equivalentDAGs
# # dagitty::a
# edges(x)


# setwd(dirname(sys.frame(1)$ofile))







source('dirichlet.R')

datadir = './assignment_1_files/'
load(file.path(datadir,"small_network_1.rda"),verbose = 1)
dat1[1:3,]
load(file.path(datadir,"small_network_2.rda"),verbose = 1)
# dat2[1:3,]

# # ?table
# 
# pars = adjacent_vertices(g, V(g), mode = c("out"))
# get_parent <- function(x)pars[[x]]
# get_parent <- function(x) adjacent_vertices(g, x, mode = c("out"))[[1]]
# g = graph_from_adjacency_matrix(aM)
# neighbors(g,V(g))
# edges(g)[[1]]
# igraph:::parent(g)
# neighbors(g,V(g))
# 

# vertex(g)
# as.data.frame(g)
dt = cbind(dat1,dat1)
dt = dat1


eta0 = 1
tb = table(dt)

# importIntoEnv
# igraph:::add_class
# is(g)
# findMethods('igraph')
# getClass('igraph')
require(igraph)

lp.nodes <- function(nodes,tb,eta=array(1,dim(tb)),...){
  tb.marg = margin.table(tb,margin = nodes)
  eta.marg = margin.table(eta,margin = nodes)
  # cat('marginal:','\n')
  # print(eta.marg)
  # print(tb.marg)
  lp = lp.diri(tb.marg,eta.marg,...)
}

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
                          # eta = array(eta0,dim(tb))
                          parNode = sess$F$get_parent(selfNode)
                          # parNode = c(1,3)
                          # selfNode = c(2)
                          nodes = c(selfNode,parNode)
                          # cat('nodes',nodes,'\n')
                          lp.nodes(nodes,sess$tb,eta,...) - lp.nodes(parNode,sess$tb,eta,...)
                        },
                        logL = function(.self,...){
                          sess = .self
                          lps = sapply( V(sess$mdlgraph), partial(sess$lp.node,...))
                          
                          sum(lps)
                        }
                      )
                    )

# adjacent_vertices(g, V(g), mode = c("in"))

source('graph_2edge.R')
sess = .PGM_binary$new(
  # dat=dat1,
  dat=dat2,
  mdlgraph=graph_from_adjacency_matrix(aM))
sess$make_table()
sess$preprocess()
sess$lp.node(1)%>%print
sess$logL()%>%print

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

# source('graph_1edge.R')
source('graph_2edge.R')
# for (aM in aMs){


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


###### Load all graphs 
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


# 
# main <- function(aM){
#   cat('\n','\n')
#   g = graph_from_adjacency_matrix(aM)
#   sess$mdlgraph = g
#   sess$preprocess()
#   # sess$lp.node(1)%>%print
#   eta0 = 1
#   lL = sess$logL(eta0=eta0,equiv=F)
#   lL %>%print
#   
#   bn <- amat2bn(aM)
#   bn$logL = lL
#   bn$iss = 2^nrow(aM) * eta0
#   bn$dat <- sess$dat
#   bn
# }
# 
# bn.list <- lapply(lst,main)
# 
# 
# {bn2df <- function(bn){
#   bnscore <- score(bn,bn$dat,'bde',iss=bn$iss)
#   bnscore <- score(bn,bn$dat,'bde',iss=8)
#   
#   list(model=modelstring(bn),myalgo.bde=bn$logL,bnlearn.bde=bnscore
#        ,bnlearn.bic=score(bn,bn$dat,'bic')
#        ,iss=bn$iss)
# }}(bn.list[[1]])
# 
# df <- Rutil::combine_args(rbind)(lapply(bn.list,bn2df))
# df

# deal::jointprior

# sess$logL

# source('plot_dag.R')
# boxed_equiv('graph_0edge.R')
# boxed_equiv('graph_1edge.R')
# boxed_equiv('graph_2edge.R')
