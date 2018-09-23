aMs = list()
{
  aM = cbind(c(0,0,0),
             c(0,0,0),
             c(1,0,0))
  ess = pcalg::dag2essgraph(aM)
  # check_aM(aM)
  aMs = c(aMs,list(aM))
}

{
  aM = cbind(c(0,0,1),
             c(0,0,0),
             c(0,0,0))
  # check_aM(aM)
  aMs = c(aMs,list(aM))
}