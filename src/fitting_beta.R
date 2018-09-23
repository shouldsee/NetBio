library(Rutil)

preview<- function(p,xs =seq(0,1,length.out = 100),silent=F
                   ,xlab = 'x',ylab='y'
                   ,...){
  dots = list(...)
  # if (is.null(dots$xlab)){ dots$xlab='x'}
  # if (is.null(dots$ylab)){ dots$ylab='y'}
  ys= p(xs)
  xlab = (substitute(xlab))
  ylab = (substitute(ylab))
  # ys = deparse(substitute(y))
  # xs = deparse(substitute(xs))
  dat <- list(x=xs,y=ys
              ,xlab=xlab
              ,ylab=ylab
  )
  # dat <- list(y=ys)
  if (silent) {return(dat)}
  # with(data,plot())
  do.call(plot,args = c(dat,dots) )
}
pbeta.ma <- function(x,m,a){pbeta(x,a,m-a)}
dbeta.ma <- function(x,m,a){dbeta(x,a,m-a)}

loss <- function(param){
  m = param[1]
  a = param[2]
  # p<-function(x) {pbeta(x,a,m-a)}
  x = c(0.25,0.75)
  y_true = c(0.05,0.95) 
  ps <- pbeta(x,a,m-a)
  Metrics::mse(y_true,ps)
  # mean((ps - y_true)^2)
}


out = optim(c(2,1),loss,method = "BFGS",control = list(maxit=300))
out$par <- as.list(setNames(out$par,c('m','a')))
# outf = partial(dbeta.ma,m = out$par$m,a = out$par$a)
outf = partial(pbeta.ma,m = out$par$m,a = out$par$a)
out$p <- outf
out$d <- partial(dbeta.ma,m = out$par$m,a = out$par$a)
res1 <- out


loss <- function(param,debug = F){
  m = param[1]
  a = param[2]
  # a = max(a,1)
  # m = max(m,2)
  # p<-function(x) {pbeta(x,a,m-a)}
  # x = c(0.25,0.75)
  # y_true = c(0.05,0.95)
  
  p = partial(pbeta.ma,m=m,a=a)
  d = partial(dbeta.ma,m=m,a=a)
  # MODE = min(max((a-1)/(m-2),0),1)
  opout = optim( c(0.4),
                 # d,
                 {function(x){-d( min(max(x,1E-4),1-1E-4) )}},
                 method = 'BFGS', control = list(maxit=500))
  MODE= opout$par
  
  
  y_true = c( 0.4, 0.1)
  y_pred = c( MODE, p(0.3)) 
  if (debug){
    # print(MODE)
    print(p(0.3))
    cat('y_pred:',y_pred,'\n')
    cat('y_true:',y_true,'\n')
  }
  Metrics::mae(y_true, y_pred)
  # mean((ps - y_true)^2)
}

system.time({
  out = optim(c(10,5),loss,
              control = list(maxit=500)
              # ,method = "BFGS"
  )
  out$par <- as.list(setNames(out$par,c('m','a')))
  # outf = partial(dbeta.ma,m = out$par$m,a = out$par$a)
  outf = partial(pbeta.ma,m = out$par$m,a = out$par$a)
  out$p <- outf
  out$d <- partial(dbeta.ma,m = out$par$m,a = out$par$a)
  res2 <- out
  print(res2$value)
  loss(unlist(res2$par),debug=T)
  
})