tb = cbind(c(55,7),
           c(16,22))
eta = array(1,dim(tb))
lp.diri(tb,eta,equiv = F)
lp.nodes
L = 0
a = lp.nodes(c(1),tb,eta,equiv=F)
L = L+a
a = lp.nodes(c(2),tb,eta,equiv=F)
L = L+a
print(L)

L = 0
a = lp.nodes(c(1,2),tb,eta,equiv=F)
b = lp.nodes(c(1),tb,eta,equiv=F)
print(a-b)
b = lp.nodes(c(2),tb,eta,equiv=F)
print(a-b)
