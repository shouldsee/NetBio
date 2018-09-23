knitr::opts_chunk$set(
  echo = F
  # ,results = 'hide' ##### Will hide tables
  ,cache = F
  ,eval = T
  ,eval.after = c('fig.cap','code')
)

options(stringsAsFactors = F)  #### Always disable this!
# options(warn=  -1)
read_wrap <-function(fname,wid = 100){
  s= readLines(fname)
  s = gsub('\t','  ',s)
  
  # sout = 
  idx = nchar(s) > wid
  
  slst <- as.list(s)
  # browser()
  slst[idx]<-lapply(slst[idx],function(x)strwrap(x,wid))
  unlist(slst)
  
}
