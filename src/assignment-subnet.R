load.assignment.data()
env.data <- list(adj.true = adj.true,expr.dat=expr.dat,genes=genes,
                 gmat = gmat_from_adjacency(adj.true))
out <- replicate(50,contraster(),simplify = F)
lst <- combine_args(rbind)(out)
save(lst, env.data, file='assignment-subnet.rdf5')