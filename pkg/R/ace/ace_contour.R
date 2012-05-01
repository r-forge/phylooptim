rm(list=ls())

name_all <- c("geo","aqui","mono","mam")
well <- c("spg", "Rcgmin", "Rvmmin", "bobyqa","L-BFGS-B","nlminb","ucminf","Nelder-Mead","nlm","CG","BFGS","newuoa")

for (j in c(1:length(name_all))){
  j <- 1
  name <- name_all[[j]]
  load(paste("/home/michels/Hallowed/repository/phylooptim/pkg/R/ace/",name,"acects.RData",sep=""))
    for (i in c(1:length(well))){
      xx <- l[[i]]$foostore$sigma2
      yy <- l[[i]]$foostore$ace
      lx <- max(xx)-min(xx)
      ly <- max(yy)-min(yy)
      x1range <- seq(min(xx)+0.0000001,max(xx),length=1000)
      x2range <- seq(min(yy),max(yy),length=1000)

      x <- as.vector(dv$data)
      phy <- dv$phy
      type <- 'continuous'
      numcol<-50 #number of colors to use
      if (!inherits(phy, "phylo"))
          stop('object "phy" is not of class "phylo".')
      if (is.null(phy$edge.length))
          stop("tree has no branch lengths")
      type <- match.arg(type, c("continuous", "discrete"))
      nb.tip <- length(phy$tip.label)
      nb.node <- phy$Nnode
      if (nb.node != nb.tip - 1)
          stop('"phy" is not rooted AND fully dichotomous.')
      if (length(x) != nb.tip)
          stop("length of phenotypic and of phylogenetic data do not match.")
      if (!is.null(names(x))) {
          if(all(names(x) %in% phy$tip.label))
            x <- x[phy$tip.label]
          else warning("the names of 'x' and the tip labels of the tree do not match: the former wfere ignored in the analysis.")
      }
      tip <- phy$edge[, 2] <= nb.tip
      dev.BM <- function(p) {
        if (p[1] < 0) return(1e100) # in case sigma² is negative
        x1 <- p[-1][phy$edge[, 1] - nb.tip]
        x2 <- numeric(length(x1))
        x2[tip] <- x[phy$edge[tip, 2]]
        x2[!tip] <- p[-1][phy$edge[!tip, 2] - nb.tip]
        -2 * (-sum((x1 - x2)^2/phy$edge.length)/(2*p[1]) -
        nb.node * log(p[1]))
      }

      res<-t(sapply(x1range,function(x,y=x2range) sapply(y,function(z) dev.BM(c(x,rep(z, length.out = nb.node))) )))

      if(any(res[!is.na(res)]<0)) {
      	res[!is.na(res) & res>=0 & res<1]<-0
      	res[!is.na(res) & res<0]<--1*log(-1*res[!is.na(res) & res<0])
      }
      res[!is.na(res) & res>0]<-log(res[!is.na(res) & res>0])

      s <- seq(length(xx)-1)

      nlevels=25
      levels=pretty(range(res,na.rm=TRUE),nlevels,finite=TRUE)
      k<-pretty(trunc(range(res,na.rm=T)),nlevels)
      wr.pal<-colorRampPalette(c("grey","black"))
      wr<-wr.pal(round(numcol*sum(k>=0)/length(k)))
      bw.pal<-colorRampPalette(c("white","grey"))
      bw<-bw.pal(round(numcol*sum(k<=0)/length(k)))
      cols<-c(bw,wr)
      rcols<-cols[round((1:(length(levels)))*length(cols)/length(levels))]

      kk <- floor(length(l[[i]]$foostore$loglik)/10)     
      newseq <- seq(0,kk*10,by=10)
      xxsig <- rep(NA,length=kk)
      yyace <- rep(NA,length=kk)
      for (k in (1:kk)){
          lnLTries <- na.omit(l[[i]]$foostore$loglik[(newseq[k]+1):newseq[k+1]])
          if (sum(attr(lnLTries,"na.action"))<sum(1:10)){
          sigma2x <- na.omit(l[[i]]$foostore$sigma2[(newseq[k]+1):newseq[k+1]])
          acey <- na.omit(l[[i]]$foostore$ace[(newseq[k]+1):newseq[k+1]])
          bestIndex<-which.min(lnLTries[1:20])
          xxsig[k] <- sigma2x[bestIndex]
          yyace[k] <- acey[bestIndex]}else{next}
      }
      
      png(file=paste("ace_cts_",name,"_",well[i],"_contour.png",sep=""))
      filled.contour(x1range,x2range,res,levels=levels,col=rcols,  
      		plot.axes={
      			axis(1)
      			axis(2)
            plot.title=title(main=paste("Ace Continuous ",name,"",well[i]),xlab="sigma2", ylab="ace")
      			contour(x1range,x2range,res,add=TRUE,nlevels=length(levels))
      			arrows(xxsig[s], yyace[s], xxsig[s+1], yyace[s+1],length=.1,col="yellow",lwd=2)
      			points(l[[i]]$sigma2,l[[i]]$ace[1],col=2,pch=15,cex=3)
            points(xx[1],yy[1],col="green",pch=19,cex=2)
      		}
      )
      dev.off()
    }
rm(list=ls()[-1*c(which(ls()=="well"),which(ls()=="name_all"),which(ls()=="j"))])
}

#Where the true min is approximately
xmin <- x1range[which(res==min(res),arr.ind=TRUE)[1]]
ymin <- x2range[which(res==min(res),arr.ind=TRUE)[2]]

#Likelihoods of "the best"
good <- NULL
for (i in c(1:length(well))){
  if(!is.na(l[[i]]$sigma2)|!is.na(l[[i]]$ace[1])){if(round(xmin,4)+.05 > round(l[[i]]$sigma2,4) & round(ymin,4)-.25 < round(l[[i]]$ace[1],4)){good <- c(good,i)}}}
llh <- data.frame(loglik=rep(NA,length(good)),sigma2=rep(NA,length(good)),ace=rep(NA,length(good)))
for (i in c(1:length(good))){llh[i,] <- c(round(l[[good[i]]]$loglik,5),round(l[[good[i]]]$sigma2,5),round(l[[good[i]]]$ace[1],5))}
rownames(llh) <- well[good]
llh
