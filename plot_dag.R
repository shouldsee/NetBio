
plot_jpeg = function(path, add=FALSE,xlef=1,ybot=1)
{
  require('jpeg')
  jpg = readJPEG(path, native=T) # read the file
  res = dim(jpg)[2:1] # get the resolution, [x, y]
  if (!add) # initialize an empty plot area if add==FALSE
  {plot(1,1,xlim=c(1,res[1]),ylim=c(1,res[2]),asp=1,type='n',xaxs='i',yaxs='i',xaxt='n',yaxt='n',xlab='',ylab='',bty='n')}
  rasterImage(jpg,xlef,ybot,xlef+res[1],ybot+res[2])
  return(c(xlef+res[1],ybot+res[2]))
}



boxed_equiv<- function (data.script) {
  # source('graph_1edge.R')
  # data.script='graph_1edge.R'
  # source('graph_2edge.R')
  source(data.script)
  {
    deco <- function(g,...){
      vertex_attr(g,'name') <- c('A','B','C')
      vertex_attr(g,'color') <- c(0,0,0)
      Lv =length(V(g))
      l = igraph::layout_on_grid(g,width = 3,height = 1.5)
      # l = igraph::layout_with_kk(g,coords = cbind(1:Lv,rep(0,Lv)),miny=rep(0,Lv),maxy=rep(0.1,Lv))
      # l = igraph::layout_with_kk(g,miny=rep(0,Lv),maxy=rep(0.1,Lv))
      # l[,2]=0.1
      l[,1]=c(1,3,2)/3
      # l[,2]=0.1
      # print(l)
      # l[,2]=1:3
      # print(l)
      plot(g,
           layout=l,
           # edge.width=E(g)$width,
           vertex.size=15,
           ...)
      g
    }
    deco(dag,rescale=T)
  }
  check_aM <- function(aM){
    require(igraph)
    # par(mfrow = c(1,2))
    par(mfrow = c(1,1))
    g = igraph::graph_from_adjacency_matrix(aM)
    # g0 = deco(g)
    ess = pcalg::dag2essgraph(aM)
    # g = igraph::graph_from_adjacency_matrix(ess)
    # g1 = deco(g)
    return(list(g,ess))
  }
  dags = Rutil::combine_args(rbind)(lapply(aMs,check_aM))
  uniq = unique(dags[,2])
  for (i in 1:length(uniq)){
    f = function(x) all(x[[2]]==uniq[[i]])
    idx = apply(dags,MARGIN=1,f)
    idx = which(idx)
    # par(mfrow=c(1,length(idx)))
    wd = 600
    # asp = 0.25
    asp = 0.125
    ht = wd*asp
    print(i)
    jpeg(paste0(i,'.jpg'),width = wd,height = ht*length(idx)+50)
    par(mfrow=c(length(idx),1))
    par(mar=c(1,0,0,0),
        cex=1.,
        oma=c(1.5,.5,.5,.5))
    for (id in (idx) ){
      dag=dags[[id,1]]
      deco(dag,asp=ht/wd/3
           ,margin=c(1.80,0,0,0)
      )
    }
    box("inner")
    dev.off()
  }
  
  {
    ht = ((wd*asp)+20)*length(aMs)+30
    fname= paste0(tools::file_path_sans_ext(data.script),'.jpg')
    # pdf( fname,width = wd/96, height = ht/96)
    jpeg( fname,width = wd, height = ht)
    par(
      mar=c(.0,0,2,0),
      cex=1.,
      oma=c(.5,.5,.5,.5))
    par(mfrow = c(1,1))
    plot(c(0,wd),c(0,ht),'n',axes = F,xlab = '',ylab='')
    res = c(0,0)
    for (i in 1:length(uniq)){
      res = plot_jpeg(path=paste0(i,'.jpg'),xlef = 0,ybot = res[2], add=T)
      # res = plot_jpeg(path=paste0(i,'.jpg'),xlef = 0,ybot = res[2],add=T)
    }
    title(data.script,line=0)
    box('inner',lty=2)
    dev.off()
  }
}