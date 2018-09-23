{
  fname = 'dream5/training data/Network 3 - E. coli/net3_expression_data.tsv'
  f <- readLines(fname)
  dat <- read.csv(fname,sep='\t')
  expr.dat.big <- as.matrix(dat)
  
  
  fname <- 'dream5/test data/DREAM5_NetworkInference_GoldStandard_Network3 - E. coli.tsv'
  dat <- read.csv(fname,sep='\t')
  true.pairs.big<-dat[dat[,3]==1,1:2]
  genes.big <- colnames(expr.dat.big)
  adj.true.big <- pair2adj(true.pairs.big ,genes = genes.big)
  
  env.data <- list(adj.true = adj.true.big,expr.dat=expr.dat.big,genes=genes.big,
                   gmat = gmat_from_adjacency(adj.true.big))
}


out <- replicate(50,contraster(),simplify = F)
lst <- combine_args(rbind)(out)
save(lst, env.data, file='e-coli-subnet.rdf5')
#load('e-coli-subnet.rdf5')