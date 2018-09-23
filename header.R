library(expm)
library(Rutil)
library(igraph)
clu <- parallel::makeCluster(16)

source('dream5.R')

# load.assignment.data()
# g.true <- graph_from_adjacency_matrix(adj.true)
