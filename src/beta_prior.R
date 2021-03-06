## Finding priors for the beta distribution

Beta distribution is the natural prior for a binomial/bernoulli distribution. Considering distribution of a binomial variable $X\gvn\theta\sim Binom(\theta)$, in order to make its marginalisd distribution $P(X) = \int P(X\gvn\theta)P(\theta).d\theta$ analytically tractable, one of the choice is to assume $\theta\sim Beta(M,\alpha)$, so that:
  
  $$
  \begin{aligned}
P(\theta) &=\frac{x^{\alpha - 1}(1-x)^{M-\alpha-1}}{B(\alpha,M-\alpha)}
\\
E(\theta) &={\alpha\over M }
\end{aligned}
$$
  
  ### Constraint 1:
  
  $$
  \begin{aligned}
P(\theta\le 0.25) = P(\theta\ge 0.75) =  0.05 \\
P(\theta\le 0.75)=0.95
\end{aligned}
$$
  
  ### Constraint 2
  
  
  $$
  \begin{aligned}
\text{argmax}_\theta(P(\theta))=0.4 \\
P(\theta\le 0.3) = 0.1
\end{aligned}
$$
  
  ```{r}

partial <- function(f, ...) {
  'Source: https://stackoverflow.com/questions/32173901/how-to-efficiently-partially-apply-a-function-in-r  by josliber '
  l <- list(...)
  function(...) {
    do.call(f, c(l, list(...)))
  }
}

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
a=1;
m=2
p<-function(x) {pbeta(x,a,m-a)}
dat = preview(p,silent = F,xlab = bquote(theta))
# do.call(plot.default,args = dat)
# with(dat,plot(x,y))

pbeta.ma <- function(x,m,a){pbeta(x,a,m-a)}
dbeta.ma <- function(x,m,a){dbeta(x,a,m-a)}

source('fitting_beta.R')

```

```{r}
out <- res1

preview(out$p,type='l'
        ,xlab = bquote(theta)
        ,ylab = bquote('P'~(theta))
)
abline(h = c(0.95,0.05),col=2,lty=2)
x = c(.25,.75)
y = out$p(x)
points(x,y)
text(x,y,sprintf('(%.3f,%.3f)',x,y),adj = -.5*c(1,1))

options(digits = 3)
tl = with(data = out$par,bquote(M == .(m) ~","~ A == .(a) ~","~ loss==.(out$value) ))
title(tl)
grid()
# plot.default
```

```{r}
out <- res2
par(mar=c(4.2,4.1,1.1,4.1))
preview(out$p,type='l'
        ,xlab = quote(theta)
        ,ylab = quote(P(Theta<theta))
)
abline(v = 0.4,col=2,lty=2)

par(new=T)
x = c(.3)
y = out$p(x)
points(x,y)
text(x,y,sprintf('(%.3f,%.3f)',x,y),adj = -.25*c(-3,2))

abline(h = y,col=2,lty=2)
tl = with(data = out$par,bquote(M == .(m) ~","~ A == .(a) ~","~ loss==.(out$value) ))
title(tl)

par(new=T)

preview(out$d,type='l',ylab='',xlab='',axes = F,
        col=3)
axis(side = 4)
mtext(side= 4,line = 2,quote(P(Theta=theta)))
grid()
```


```{r}

obsS = '011100101101'
obsB <- strsplit(obsS,'')[[1]] =='1'

out <-res2

ob <-list()
##### Prior: 
prior <- out$d
ob$Fparam$prior <- prior

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

L = logL_maker(obsB)
logL = partial(L,log = T)

ob$Fparam$L <- L
ob$Fparam$logL <- logL

ob = do.call(.infer$new,ob)

# margL = logL
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


{
  make_dpos(ob)
  marg_obs(ob)
  make_post(ob)
  # ob
  
}
q2infer <- ob
# margL <- function(marg){
#   integrate(dpos,0,1)$value
#   }
# post <- function(p1){ 
#   L(p1)*prior(p1)/margL(p1) }

```

```{r bayes,fig.cap = cap}
cap = paste('Inference of posterior distribution on parameter $\\theta$ given mutatble sequence:', obsS)
par(lwd=2)
par(mar=c(4.2,4.1,1.1,4.1))
i = 1

ob <- q2infer

preview(ob$Fparam$prior,type= 'l'
        ,xlab = quote(theta)
        ,ylab = quote(P)
        ,lty=i
        ,col=i
        ,ylim=c(0,6)
        ,pch = i
        # ,lwd = 2
)

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

par(new = T)
preview(ob$Fparam$L,type= 'b',silent = F
        ,lty=i
        ,col=i
        ,pch=i
        # ,lwd = 2
        
        ,axes = F,xlab='',ylab='')

mtext(side=4,line=2,'Likelihood');axis(side = 4)

grid()
legend(0, .125, 
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
       lty=1:i, cex=0.8)
# do.call(lines,c(dat,col=2))
# preview(post,type='l')


# outcounts
```

```{r}
#
# ob$Fparam$post()
# preview(ob$Fparam$prior)

post <- ob$Fparam$post
L = ob$Fparam$L
preview(L,silent = F)
par(new = T)
preview(post,silent = F)
post <- bayes_infer(post,L)
combine_args(lines)(preview(post,silent = T))
post <- bayes_infer(post,L)
combine_args(lines)(preview(post,silent = T))
post <- bayes_infer(post,L)

# par(new = T)
# combine_args(plot)(c(preview(post,silent = T),list(log='y')))
# post <- bayes_infer(post,L)

```


```{r}
# print(ob$marg$obs)
integrate(ob$Fparam$dpos,0,1)
integrate(ob$Fparam$L,0,1)
integrate(ob$Fparam$prior,0,1)

tol = .Machine$double.eps^0.5
# pryr::e
a = integrate(function(x) ob$Fparam$L(x)*ob$Fparam$prior(x),0,1,rel.tol = tol) 
a
b = integrate(ob$Fparam$L,0,1,rel.tol = tol)
b
a$value/b$value

```

### Basics for Bayesian inference \label{sec:bayes-basics}

$$
  \begin{aligned}
\text{likelihood}&: f(\theta) = P(x \gvn \theta) \\
\text{prior}&: f(\theta) = P(\theta) \\
\text{posterior}&: f(\theta)  = P(\theta \gvn x) = \frac{  P(x\gvn \theta) P(\theta)  }{P(x) } \\
\text{marginal likelihood} &: P(x) =  \int P(x\gvn \theta) P(\theta).d\theta
\end{aligned}
$$
  
  The observed toss sequence is: `r obsS`

Assume each coin toss is independent from each other, the likelihood of an observed sequence is only dependent on the total number of heads and not the sequence it occurred in. Denoting the coin toss as a sequence $\{x_i\}$ where $x_i \in \{0,1\}$, we have

$$
  \begin{aligned}
\text{\#head} = \indicator\{x_i=1\} \\
  \text{\#head} \sim Binom(|\{x_i\}|,\theta)
    \end{aligned}
    $$
      Combining with the prior $\theta\sim$ Beta(`r res2$par$a`, `r with(res2$par,m-1)` ), I calculated the marginal likelihood numerically to be $P(x)=$ `r q2infer$marg$obs` and derived the posterior distribution accordingly (see figure \ref{fig:bayes}).
      
      