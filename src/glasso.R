library(glasso)
?glasso
# 

# x<-matrix(rnorm(50*20),ncol=20)
# s<- var(x)
# a<-glasso(s, rho=.01)
# aa<-glasso(s,rho=.02, w.init=a$w, wi.init=a$wi)


s <- var(expr.dat)
image(a$wi)
a<-glasso(s, rho=.01)
aa<-glasso(s,rho=.02, w.init=a$w, wi.init=a$wi)
aaa<-glasso(s,rho=.001, w.init=a$w, wi.init=a$wi)
aaa<-glasso(s,rho=.01)
par(mfrow=c(3,3))
pipeline(a$wi%>%abs)
pipeline(aa$wi%>%abs)
pipeline(aaa$wi%>%abs)
# image


pipeline(a$w%>%abs)
pipeline(aa$w%>%abs)
