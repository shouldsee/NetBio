
#### Likelihood: 
logL_maker <- function(obsB){
  N = length(obsB)
  X = sum(obsB)
  function(p1,...) {
    dbinom(x= X,size = N,p1,...)}
  # return( )
  # d(data)
  # pa * obj$prior
}

.infer<- setRefClass('infer',fields = list(
  # prior='function',
  Fparam = 'list',
  # dpos = 'function', ### Un-normalised posterior
  marg = 'list'
))

make_dpos <- function(ob){
  ob$Fparam$dpos <- function(x){ob$Fparam$L(x)*ob$Fparam$prior(x)}
}

marg_obs <- function(ob){
  f = ob$Fparam$dpos
  val = integrate( f, 0,1)$value
  ob$marg = list(obs=val)
  val
}

bayes_infer <- function(prior,L,low=0,high=1){
  # cnorm = 
  post_raw <- function(PAR){prior(PAR)*L(PAR)}
  res = integrate(post_raw,low,high)
  cnorm = res$value
  post <- function(PAR){post_raw(PAR)/cnorm}
  
}
make_post <- function(ob){
  ob$Fparam$post <- bayes_infer(ob$Fparam$prior,ob$Fparam$L)
  # ob$Fparam$post <- function(PAR){ob$Fparam$L(PAR)*ob$Fparam$prior(PAR)/ob$marg$obs}
}


plot_bayes<-function(ob,ylab1=quote(P),ylab2='Likelihood') {
  # par(mar=c(4.2,4.1,1.1,4.1))
  par(lwd=2)
  i = 1
  preview(ob$Fparam$prior,type= 'l'
          ,xlab = quote(theta)
          ,ylab=''
          ,lty=i
          ,col=i
          ,ylim=c(0,6)
          ,pch = i
          # ,lwd = 2
  )
  mtext(ylab1,side=2,line=2)
  
  i = i+1
  # par(new = T)
  dat = preview(ob$Fparam$post,type= 'l',silent = T)
  do.call(lines,c(dat
                  ,lty=i
                  ,col=i
                  ,pch = i
                  # ,lwd = 2
  ))
  # text(labels = c("red line"))
  
  i = i+1
  legend(0, 2.25, 
         legend=c(
           # substitute(bquote(P(theta))),
           # deparse(substitute(P(theta))),
           expression(P(theta)),
           expression(P(theta~'|'~x)),
           expression(P(x~'|'~theta))),
         # "Prior",
         # 'Posterior',"Likelihood"),
         pch =c(NA,NA,3),
         col=1:i,
         lty=1:i, cex=0.6)
  
  par(new = T)
  preview(ob$Fparam$L,type= 'b',silent = F
          ,lty=i
          ,col=i
          ,pch=i
          # ,lwd = 2
          
          ,axes = F,xlab='',ylab='')
  
  mtext(side=4,line=2,ylab2);axis(side = 4)
  
  grid()

}

obsS = '011100101101'
obsB <- strsplit(obsS,'')[[1]] =='1'

out <-res2

ob <-list()
##### Prior: 
prior <- out$d
ob$Fparam$prior <- prior


L = logL_maker(obsB)
logL = partial(L,log = T)

ob$Fparam$L <- L
ob$Fparam$logL <- logL

ob = do.call(.infer$new,ob)


{
  make_dpos(ob)
  marg_obs(ob)
  make_post(ob)
  # ob
  
}
q2infer <- ob