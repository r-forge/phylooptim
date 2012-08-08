#Note: This is for the Aquilegia data set

require(ouch)
require(geiger)

load("ouch.contour.sv.RData")

#sqrt.alpha
x1range1 <- seq(0,2,length=100)
x1range <- seq(2.05,5,length=100)
x1range <- c(x1range1,x1range)

#sigma
x2range1 <- seq(0,2,length=100)
x2range <- seq(2.05,5,length=100)
x2range <- c(x2range1,x2range)

#The response
  res<-t(sapply(x1range,function(x,y=x2range) sapply(y,function(z) ou.lik.fn(
                               tree=tree,
                               alpha=x^2,
                               sigma=z^2,
                               beta=beta,
                               dat=dat
                               )$deviance/-2)))

#Needed for a filled contour plot, in grayscale
      numcol<-50
      wr.pal<-colorRampPalette(c("white","grey"))
      wr <- wr.pal(numcol)
      bw.pal<-colorRampPalette(c("grey","black"))
      bw <- bw.pal(numcol)
      cols<-c(wr,bw)

x.labels <- c(0:5)
y.labels <- c(0:5)

lev <- c(seq(-30-83*.5,-30.5,by=.5),seq(-30,-22,length.out=17))

for (i in c(1:length(well))){
  for (j in c(1:length(l[[i]][[2]]))){
    l[[i]][[2]][[j]][[1]]$sqrt.alpha <- abs(l[[i]][[2]][[j]][[1]]$sqrt.alpha)
    l[[i]][[2]][[j]][[1]]$sigma <- abs(l[[i]][[2]][[j]][[1]]$sigma)
  }
}

sv.correct.list <- list()
sv.correct.list[[1]] <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25)
sv.correct.list[[2]] <- c(1, 2, 6, 7, 8, 9, 10, 12, 13, 14, 15, 17, 18, 19, 20, 23, 24, 25)
sv.correct.list[[3]] <- c(1, 2, 3, 4, 7, 8, 9, 10, 12, 13, 14, 15, 17, 18, 19, 20, 23, 24, 25)
sv.correct.list[[4]] <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25)
sv.correct.list[[5]] <- c(1, 2, 3, 4, 5, 7, 8, 9, 13, 14, 18, 19)
sv.correct.list[[6]] <- c(1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23)
sv.correct.list[[7]] <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25)
sv.correct.list[[8]] <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24)
sv.correct.list[[9]] <- c(1, 5, 7, 8, 9, 10, 13, 14, 15, 18, 19, 20, 24, 25)
sv.correct.list[[10]] <- c(1, 5, 7, 8, 9, 10, 13, 14, 15, 18, 19)
sv.correct.list[[11]] <- c(1, 5, 7, 8, 9, 10, 13, 14, 15, 19, 20, 24)
sv.correct.list[[12]] <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 23, 24, 25)

for (j in c(1:length(well))){
  png(file=paste(well[j],".contour.png",sep=""))
  filled.contour2(x1range,
                  x2range,
                  res,
                  levels=lev,
                  col=cols,
                  axes=FALSE,
        	  plot.axes = {
                              axis(side = 1, at = x.labels,labels=x.labels)
                              axis(side = 2, at = y.labels,labels=y.labels)
                              plot.title = title(xlab="sqrt.alpha",
                                ylab="sigma", main=paste(well[j]),cex.main=1.5)

                  #First start value
                  for (i in c(1:length(l[[j]][[2]]))[c(1:25)]){
                    points(l[[j]][[2]][[i]][[1]]$sqrt.alpha[1],
                           l[[j]][[2]][[i]][[1]]$sigma[1],
                           col="black",
                           pch=20,
                           cex=2,
                           lwd=2)
                  }

                  #End start value
                  for (i in c(1:length(l[[j]][[2]]))[sv.correct.list[[j]]]){
                    points(l[[j]][[2]][[i]][[1]]$sqrt.alpha[1],
                           l[[j]][[2]][[i]][[1]]$sigma[1],
                           col="white",
                           pch=20,
                           cex=1.25,
                           lwd=2)
                  }
                              }
                  )
  dev.off()
}
