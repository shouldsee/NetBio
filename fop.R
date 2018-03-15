### reference: http://adv-r.had.co.nz/Function-operators.html


file.edit('~/Rutil/R/fop.R')
file.edit('~/Rutil/compile.R')
f = (function(x = 4) g(x) + h(x))

library(dplyr)
{
  formals(f) %>% print
  body(f) %>% print
  environment(f) %>% print
}
library(Rutil)
install.packages.lazy('pryr')
library(pryr)
pryr::compose

{
  f = function(x) x ^ 2
  g = function(x,y=1) x + 1
  # Reduce(compose())
  a=pryr::compose(f,g)(2)
  a%>%print
  b=Rutil::compose(f,g)(2) 
  b%>%print
  # Reduce(compose)
}

# pryr::compose

{
  xs = rep(c(f,g),c(1,500))
  x0=2
  #### Recursion size is beyond stack size easily in the first one
  system.time(Reduce(Rutil::compose,xs)(x0)%>%print)%>%print
  ####
  system.time({Rutil::combine_args(pryr::compose)(rev(xs))(x0)%>%print})%>%print
  #### 
  system.time({Rutil::combine_args(Rutil::compose)(c(xs,right=F) )(x0)%>%print})%>%print
}


{
  ###### Reduce is slower than for loop by a factor of 2
  xs = 1:1000000
  system.time({y=Reduce(g,xs,init=0)})%>%print
  print(y)
  system.time({
    y = 0
    for (x in xs){
      y = g(y,x)
    }
  })%>%print
  print(y)
}

