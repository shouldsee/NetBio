require(Rutil)
#' @export
lp.diri <- function(n,eta,equiv = T)
{
  # require(Rutil)
  Sn = sum(n)
  Se = sum(eta)
  # ?gamma
  # lgamma(c(2,3,4))
  logP = lgamma(Se) - lgamma(Se + Sn) + sum(lgamma(n+eta) - lgamma(eta))
  if (equiv){
    ### apply equivalency correction 
    logP = logP + ndchoose(n,log = T)
  }
  logP
}