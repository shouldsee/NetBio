{
  library(Rutil)
  install.packages.lazy(c('deal','bnlearn'))
  library(deal)
  dat = dat1
  g.start = network(dat)
  g.prior = jointprior(g.start,N = 8)
  nets = networkfamily(dat,g.start,g.prior)[[1]]
  nets = nwfsort(nets)
  o = head(nets)
  class(o)<-'networkfamily'
  print(o)
  
  dat = dat2
  g.start = network(dat)
  g.prior = jointprior(g.start,N = 8)
  nets = networkfamily(dat,g.start,g.prior)[[1]]
  nets = nwfsort(nets)
  o = head(nets)
  class(o)<-'networkfamily'
  print(o)
}

library(bnlearn)
data("learning.test")
# hc
ntw.hc <- hc(learning.test,score='bde')
arcs(ntw.hc)
ntw.hc

ntw.tabu <- tabu(learning.test,score='bde')
# arcs(ntw.tabu)
ntw.tabu

dat <- dat2
res <- hc(dat,score='bde')

?hc

data(learning.test)
dag = hc(learning.test, score = "bic")

for (i in 1:3) {
  
  a = alpha.star(dag, learning.test)
  dag = hc(learning.test, score = "bde", iss = a)
  print(a)
}
res <- dag
dat <- learning.test
score( res, dat, type = "bde")

dag = hc(learning.test, score = "bde")
score( res, dat, type = "bde")

score( res, dat, type = "loglik",iss=8)
bnlearn:::check.label
bnlearn:::discrete.data.types
plot(res)


#' Convert an adjacency matrix to a "bnlearn" object
#' @export
amat2bn <- function(aM,nodes = rownames(aM)){
  if (is.null(nodes)){
    nodes = LETTERS[1:nrow(aM)]
  }
  g0 <- empty.graph(nodes)
  amat(g0) <- aM
  g0
}

bn = bnlearn::pc.stable(bn.list[[1]]$dat)
# modelstring(bn)
plot(bn)
# deal
# deal::
# bnlearn::h

# library(GeneNet,lib.loc = '/data/subliminal/RLib')
install.packages.lazy(c('GENIE3','GeneNet'))
BiocInstaller::biocLite('GENIE3')

# bnlearn::pri
# dg = pcalg::pcalg2dagitty(aM, V(g), type = "cpdag")
# dagitty::
# print(dg)
