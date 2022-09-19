# heatmap of correlation matrix
heatmap.cor <- function(Cmat,main=NULL) {
  require(cluster)
  require(grid)
  require(ggplot2)
  dist <- 1 - Cmat
  clust  <- agnes(x = dist, diss = TRUE)
  orderLab <- clust$order.lab
  Cmat <- Cmat[rev(orderLab), rev(orderLab)]
  
  Cmat[!lower.tri(Cmat)] <- NA
  hCols <- rev(rainbow(256, start = 0, end = 2/3))
  n <- nrow(Cmat)
  hh <- rev(rainbow(256, start = 0, end = 2/3))
  cor.df <- data.frame(cor=as.vector(Cmat),x=rep(1:n,each=n),y=rep(1:n,n))
  oo <- ggplot(data=cor.df,aes(x=x,y=rev(y),fill=cor)) + geom_tile() + 
    scale_fill_gradientn(limits=c(-1,1),colours=hh,na.value='white') 
  oo <- oo + scale_x_continuous(breaks=1:n,labels=rev(orderLab)) + 
    scale_y_continuous(breaks=1:n,labels=orderLab) + xlab("ENV") + ylab("ENV") +
    theme(aspect.ratio=1,axis.text=element_text(size=10,face='bold'),axis.text.x = element_text(angle=90, vjust=1),
      axis.title=element_text(size=10), legend.key.height=unit(2,'cm'), legend.key.width=unit(0.8,'cm'))
  if(!is.null(main))
    oo <- oo + ggtitle(main) 
  return(oo)
}
