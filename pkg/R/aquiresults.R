
R version 2.14.2 (2012-02-29)
Copyright (C) 2012 The R Foundation for Statistical Computing
ISBN 3-900051-07-0
Platform: x86_64-pc-linux-gnu (64-bit)

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

  Natural language support but running in an English locale

R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.

Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

> .help.ESS <- help
> options(STERM='iESS', editor='emacsclient')
> rm(list = ls())
> 
> require(gplots)
Loading required package: gplots
Loading required package: gtools
Loading required package: gdata
gdata: read.xls support for 'XLS' (Excel 97-2004) files ENABLED.

gdata: read.xls support for 'XLSX' (Excel 2007+) files ENABLED.

Attaching package: ‘gdata’

The following object(s) are masked from ‘package:stats’:

    nobs

The following object(s) are masked from ‘package:utils’:

    object.size

Loading required package: caTools
Loading required package: bitops
Loading required package: grid
Loading required package: KernSmooth
KernSmooth 2.23 loaded
Copyright M. P. Wand 1997-2009

Attaching package: ‘gplots’

The following object(s) are masked from ‘package:stats’:

    lowess

> 
> p <- function(x, digits=4, prefix="", cex.cor, ...)
+ {
+     r <- x$MLE
+     txt <- format(c(r, 0.123456789), digits=digits)
+     txt <- paste(prefix, txt, sep="")[1:12]
+ #    if(missing(cex.cor)) cex.cor <- 0.5/strwidth(txt)
+ #    text(0.5, 0.5, txt, cex = cex.cor )
+     return(txt)
+ }
> 
> #Results for aquilegia data
> load("/home/michels/Hallowed/repository/phylooptim/pkg/R/geiger/aquigeigercts.RData")
> geiger_cts <- l
> l <- NULL
> geiger_cts$mean_lik <- round(c(mean(as.numeric(geiger_cts[[1]][,1])),mean(as.numeric(geiger_cts[[2]][,1])),mean(as.numeric(geiger_cts[[3]][,1])),mean(as.numeric(geiger_cts[[4]][,1])),mean(as.numeric(geiger_cts[[5]][,1])),mean(as.numeric(geiger_cts[[6]][,1])),mean(as.numeric(geiger_cts[[7]][,1])),mean(as.numeric(geiger_cts[[8]][,1])),mean(as.numeric(geiger_cts[[9]][,1])),mean(as.numeric(geiger_cts[[10]][,1])),mean(as.numeric(geiger_cts[[11]][,1])),mean(as.numeric(geiger_cts[[12]][,1]))),4)
> names(geiger_cts$mean_lik) <- names(geiger_cts)[1:12]
> geiger_cts$MLE <- round(c(max(as.numeric(geiger_cts[[1]][,1])),max(as.numeric(geiger_cts[[2]][,1])),max(as.numeric(geiger_cts[[3]][,1])),max(as.numeric(geiger_cts[[4]][,1])),max(as.numeric(geiger_cts[[5]][,1])),max(as.numeric(geiger_cts[[6]][,1])),max(as.numeric(geiger_cts[[7]][,1])),max(as.numeric(geiger_cts[[8]][,1])),max(as.numeric(geiger_cts[[9]][,1])),max(as.numeric(geiger_cts[[10]][,1])),max(as.numeric(geiger_cts[[11]][,1])),max(as.numeric(geiger_cts[[12]][,1]))),4)
> names(geiger_cts$MLE) <- names(geiger_cts)[1:12]
> 
> load("/home/michels/Hallowed/repository/phylooptim/pkg/R/geiger/aquigeigerdisc.RData")
> geiger_disc <- l
> l <- NULL
> geiger_disc$mean_lik <- round(c(mean(as.numeric(geiger_disc[[1]][,7])),mean(as.numeric(geiger_disc[[2]][,7])),mean(as.numeric(geiger_disc[[3]][,7])),mean(as.numeric(geiger_disc[[4]][,7])),mean(as.numeric(geiger_disc[[5]][,7])),mean(as.numeric(geiger_disc[[6]][,7])),mean(as.numeric(geiger_disc[[7]][,7])),mean(as.numeric(geiger_disc[[8]][,7])),mean(as.numeric(geiger_disc[[9]][,7])),mean(as.numeric(geiger_disc[[10]][,7])),mean(as.numeric(geiger_disc[[11]][,7])),mean(as.numeric(geiger_disc[[12]][,7]))),4)
> names(geiger_disc$mean_lik) <- names(geiger_disc)[1:12]
> geiger_disc$MLE <- round(c(max(as.numeric(geiger_disc[[1]][,7])),max(as.numeric(geiger_disc[[2]][,7])),max(as.numeric(geiger_disc[[3]][,7])),max(as.numeric(geiger_disc[[4]][,7])),max(as.numeric(geiger_disc[[5]][,7])),max(as.numeric(geiger_disc[[6]][,7])),max(as.numeric(geiger_disc[[7]][,7])),max(as.numeric(geiger_disc[[8]][,7])),max(as.numeric(geiger_disc[[9]][,7])),max(as.numeric(geiger_disc[[10]][,7])),max(as.numeric(geiger_disc[[11]][,7])),max(as.numeric(geiger_disc[[12]][,7]))),4)
> names(geiger_disc$MLE) <- names(geiger_disc)[1:12]
> 
> load("/home/michels/Hallowed/repository/phylooptim/pkg/R/ace/aquiacects.RData")
> ace_cts <- l
> l <- NULL
> ace_cts$mean_lik <- round(c(mean(ace_cts[[1]][,5]),mean(ace_cts[[2]][,5]),mean(ace_cts[[3]][,5]),mean(ace_cts[[4]][,5]),mean(ace_cts[[5]][,5]),mean(ace_cts[[6]][,5]),mean(ace_cts[[7]][,5]),mean(ace_cts[[8]][,5]),mean(ace_cts[[9]][,5]),mean(ace_cts[[10]][,5]),mean(ace_cts[[11]][,5]),mean(ace_cts[[12]][,5])),4)
> names(ace_cts$mean_lik) <- names(ace_cts)[1:12]
> ace_cts$MLE <- round(c(max(ace_cts[[1]][,5]),max(ace_cts[[2]][,5]),max(ace_cts[[3]][,5]),max(ace_cts[[4]][,5]),max(ace_cts[[5]][,5]),max(ace_cts[[6]][,5]),max(ace_cts[[7]][,5]),max(ace_cts[[8]][,5]),max(ace_cts[[9]][,5]),max(ace_cts[[10]][,5]),max(ace_cts[[11]][,5]),max(ace_cts[[12]][,5])),4)
> names(ace_cts$MLE) <- names(ace_cts)[1:12]
> 
> load("/home/michels/Hallowed/repository/phylooptim/pkg/R/ace/aquiacedisc.RData")
> ace_disc <- l
> l <- NULL
> ace_disc$mean_lik <- round(c(mean(ace_disc[[1]][,5]),mean(ace_disc[[2]][,5]),mean(ace_disc[[3]][,5]),mean(ace_disc[[4]][,5]),mean(ace_disc[[5]][,5]),mean(ace_disc[[6]][,5]),mean(ace_disc[[7]][,5]),mean(ace_disc[[8]][,5]),mean(ace_disc[[9]][,5]),mean(ace_disc[[10]][,5]),mean(ace_disc[[11]][,5]),mean(ace_disc[[12]][,5])),4)
> names(ace_disc$mean_lik) <- names(ace_disc)[1:12]
> ace_disc$MLE <- round(c(max(ace_disc[[1]][,5]),max(ace_disc[[2]][,5]),max(ace_disc[[3]][,5]),max(ace_disc[[4]][,5]),max(ace_disc[[5]][,5]),max(ace_disc[[6]][,5]),max(ace_disc[[7]][,5]),max(ace_disc[[8]][,5]),max(ace_disc[[9]][,5]),max(ace_disc[[10]][,5]),max(ace_disc[[11]][,5]),max(ace_disc[[12]][,5])),4)
> names(ace_disc$MLE) <- names(ace_disc)[1:12]
> 
> load("/home/michels/Hallowed/repository/phylooptim/pkg/R/ouch/aquiouchcts.RData")
> ouch <- l
> l <- NULL
> ouch$mean_lik <- round(c(mean(ouch[[1]][,5]),mean(ouch[[2]][,5]),mean(ouch[[3]][,5]),mean(ouch[[4]][,5]),mean(ouch[[5]][,5]),mean(ouch[[6]][,5]),mean(ouch[[7]][,5]),mean(ouch[[8]][,5]),mean(ouch[[9]][,5]),mean(ouch[[10]][,5]),mean(ouch[[11]][,5]),mean(ouch[[12]][,5])),4)
> names(ouch$mean_lik) <- names(ouch)[1:12]
> ouch$MLE <- round(c(max(ouch[[1]][,5]),max(ouch[[2]][,5]),max(ouch[[3]][,5]),max(ouch[[4]][,5]),max(ouch[[5]][,5]),max(ouch[[6]][,5]),max(ouch[[7]][,5]),max(ouch[[8]][,5]),max(ouch[[9]][,5]),max(ouch[[10]][,5]),max(ouch[[11]][,5]),max(ouch[[12]][,5])),4)
> names(ouch$MLE) <- names(ouch)[1:12]
> 
> #True Liklihood
> tglik_cts <- round(max(na.omit(as.numeric(geiger_cts$MLE))),4)
> gliktf_cts <- rep(NA,length(geiger_cts$mean_lik))
> for (i in c(1:length(gliktf_cts))){if (is.na(geiger_cts$mean_lik[i])){gliktf_cts[i] <- "NA"}else{if (as.numeric(geiger_cts$mean_lik[i])==tglik_cts){gliktf_cts[i] <- "Y"}else{if (tglik_cts-.05 < as.numeric(geiger_cts$mean_lik[i]) && tglik_cts+.05 > as.numeric(geiger_cts$mean_lik[i])){gliktf_cts[i] <- "O"}else{gliktf_cts[i] <- "--"}}}}
> 
> for (i in c(1:length(well))){
+   prptn <- rep(NA,length(geiger_cts[[i]][,1]))
+     for (j in c(1:length(prptn))){
+       if (is.na(geiger_cts[[i]][,1][j])){prptn[j] <- 0}else{if (round(as.numeric(geiger_cts[[i]][,1]),4)[j]==tglik_cts){prptn[j] <- 1}else{prptn[j] <- 0}}}
+   geiger_cts[[i]] <- cbind(geiger_cts[[i]],prptn)
+ }
> 
> gprptn_cts <- c(round(tglik_cts,3),round(mean(geiger_cts[[1]][,10]),3),round(mean(geiger_cts[[2]][,10]),3),round(mean(geiger_cts[[3]][,10]),3),round(mean(geiger_cts[[4]][,10]),3),round(mean(geiger_cts[[5]][,10]),3),round(mean(geiger_cts[[6]][,10]),3),round(mean(geiger_cts[[7]][,10]),3),round(mean(geiger_cts[[8]][,10]),3),round(mean(geiger_cts[[9]][,10]),3),round(mean(geiger_cts[[10]][,10]),3),round(mean(geiger_cts[[11]][,10]),3),round(mean(geiger_cts[[12]][,10]),3))
> 
> tglik_disc <- round(max(na.omit(as.numeric(geiger_disc$MLE))),4)
> gliktf_disc <- rep(NA,length(geiger_disc$mean_lik))
> for (i in c(1:length(gliktf_disc))){if (is.na(geiger_disc$mean_lik[i])){gliktf_disc[i] <- "NA"}else{if (as.numeric(geiger_disc$mean_lik[i])==tglik_disc){gliktf_disc[i] <- "Y"}else{if (tglik_disc-.05 < as.numeric(geiger_disc$mean_lik[i]) && tglik_disc+.05 > as.numeric(geiger_disc$mean_lik[i])){gliktf_disc[i] <- "O"}else{gliktf_disc[i] <- "--"}}}}
> 
> for (i in c(1:length(well))){
+   prptn <- rep(NA,length(geiger_disc[[i]][,1]))
+     for (j in c(1:length(prptn))){
+       if (is.na(geiger_disc[[i]][,7][j])){prptn[j] <- 0}else{if (round(as.numeric(geiger_disc[[i]][,7]),4)[j]==tglik_disc){prptn[j] <- 1}else{prptn[j] <- 0}}}
+   geiger_disc[[i]] <- cbind(geiger_disc[[i]],prptn)
+ }
> 
> gprptn_disc <- c(round(tglik_disc,3),round(mean(geiger_disc[[1]][,10]),3),round(mean(geiger_disc[[2]][,10]),3),round(mean(geiger_disc[[3]][,10]),3),round(mean(geiger_disc[[4]][,10]),3),round(mean(geiger_disc[[5]][,10]),3),round(mean(geiger_disc[[6]][,10]),3),round(mean(geiger_disc[[7]][,10]),3),round(mean(geiger_disc[[8]][,10]),3),round(mean(geiger_disc[[9]][,10]),3),round(mean(geiger_disc[[10]][,10]),3),round(mean(geiger_disc[[11]][,10]),3),round(mean(geiger_disc[[12]][,10]),3))
>             
> talik_cts <- round(max(na.omit(as.numeric(ace_cts$MLE))),4)
> aliktf_cts <- rep(NA,length(ace_cts$mean_lik))
> for (i in c(1:length(aliktf_cts))){if (is.na(ace_cts$mean_lik[i])){aliktf_cts[i] <- "NA"}else{if (as.numeric(ace_cts$mean_lik[i])==talik_cts){aliktf_cts[i] <- "Y"}else{if (talik_cts-.05 < as.numeric(ace_cts$mean_lik[i]) && talik_cts+.05 > as.numeric(ace_cts$mean_lik[i])){aliktf_cts[i] <- "O"}else{aliktf_cts[i] <- "--"}}}}
> 
> for (i in c(1:length(well))){
+   prptn <- rep(NA,length(ace_cts[[i]][,1]))
+     for (j in c(1:length(prptn))){if (is.na(ace_cts[[i]][,5][j])){prptn[j] <- 0}else{if (round(ace_cts[[i]][,5],4)[j]==talik_cts){prptn[j] <- 1}else{prptn[j] <- 0}}}
+   ace_cts[[i]] <- cbind(ace_cts[[i]],prptn)
+ }
> aprptn_cts <- c(round(talik_cts,3),round(mean(ace_cts[[1]][,8]),3),round(mean(ace_cts[[2]][,8]),3),round(mean(ace_cts[[3]][,8]),3),round(mean(ace_cts[[4]][,8]),3),round(mean(ace_cts[[5]][,8]),3),round(mean(ace_cts[[6]][,8]),3),round(mean(ace_cts[[7]][,8]),3),round(mean(ace_cts[[8]][,8]),3),round(mean(ace_cts[[9]][,8]),3),round(mean(ace_cts[[10]][,8]),3),round(mean(ace_cts[[11]][,8]),3),round(mean(ace_cts[[12]][,8]),3))
> 
> talik_disc <- round(max(na.omit(as.numeric(ace_disc$MLE))),4)
> aliktf_disc <- rep(NA,length(ace_disc$mean_lik))
> for (i in c(1:length(aliktf_disc))){if (is.na(ace_disc$mean_lik[i])){aliktf_disc[i] <- "NA"}else{if (as.numeric(ace_disc$mean_lik[i])==talik_disc){aliktf_disc[i] <- "Y"}else{if (talik_disc-.05 < as.numeric(ace_disc$mean_lik[i]) && talik_disc+.05 > as.numeric(ace_disc$mean_lik[i])){aliktf_disc[i] <- "O"}else{aliktf_disc[i] <- "--"}}}}
> 
> for (i in c(1:length(well))){
+   prptn <- rep(NA,length(ace_disc[[i]][,1]))
+     for (j in c(1:length(prptn))){if (is.na(ace_disc[[i]][,5][j])){prptn[j] <- 0}else{if (round(ace_disc[[i]][,5],4)[j]==talik_disc){prptn[j] <- 1}else{prptn[j] <- 0}}}
+   ace_disc[[i]] <- cbind(ace_disc[[i]],prptn)
+ }
> aprptn_disc <- c(round(talik_disc,3),round(mean(ace_disc[[1]][,8]),3),round(mean(ace_disc[[2]][,8]),3),round(mean(ace_disc[[3]][,8]),3),round(mean(ace_disc[[4]][,8]),3),round(mean(ace_disc[[5]][,8]),3),round(mean(ace_disc[[6]][,8]),3),round(mean(ace_disc[[7]][,8]),3),round(mean(ace_disc[[8]][,8]),3),round(mean(ace_disc[[9]][,8]),3),round(mean(ace_disc[[10]][,8]),3),round(mean(ace_disc[[11]][,8]),3),round(mean(ace_disc[[12]][,8]),3))
>             
> tolik <- round(max(na.omit(as.numeric(ouch$MLE))),4)
> oliktf <- rep(NA,length(ouch$mean_lik))
> for (i in c(1:length(oliktf))){if (is.na(ouch$mean_lik[i])){oliktf[i] <- "NA"}else{if (as.numeric(ouch$mean_lik[i])==tolik){oliktf[i] <- "Y"}else{if (tolik-.05 < as.numeric(ouch$mean_lik[i]) && tolik+.05 > as.numeric(ouch$mean_lik[i])){oliktf[i] <- "O"}else{oliktf[i] <- "--"}}}}
> 
> for (i in c(1:length(well))){
+   prptn <- rep(NA,length(ouch[[i]][,1]))
+     for (j in c(1:length(prptn))){if (is.na(ouch[[i]][,5][j])){prptn[j] <- 0}else{if (round(ouch[[i]][,5],4)[j]==tolik){prptn[j] <- 1}else{prptn[j] <- 0}}}
+   ouch[[i]] <- cbind(ouch[[i]],prptn)
+ }
> oprptn <- c(round(tolik,3),round(mean(ouch[[1]][,10]),3),round(mean(ouch[[2]][,10]),3),round(mean(ouch[[3]][,10]),3),round(mean(ouch[[4]][,10]),3),round(mean(ouch[[5]][,10]),3),round(mean(ouch[[6]][,10]),3),round(mean(ouch[[7]][,10],3)),round(mean(ouch[[8]][,10]),3),round(mean(ouch[[9]][,10]),3),round(mean(ouch[[10]][,10]),3),round(mean(ouch[[11]][,10]),3),round(mean(ouch[[12]][,10]),3))
> smoke <- matrix(c(gliktf_cts,aliktf_cts,oliktf,gliktf_disc,aliktf_disc),ncol=12,byrow=TRUE)
> colnames(smoke) <- names(geiger_cts)[1:12]
> rownames(smoke) <- c("geiger_cts","ace_cts","ouch","geiger_disc","ace_disc")
> aquilegia <- as.table(smoke)
> aquilegia
            spg Rcgmin Rvmmin bobyqa L-BFGS-B nlminb ucminf Nelder-Mead nlm CG
geiger_cts  Y   Y      Y      Y      Y        Y      Y      Y           Y   Y 
ace_cts     --  --     NA     Y      --       Y      Y      --          Y   O 
ouch        --  --     --     --     --       --     --     --          --  --
geiger_disc --  --     --     O      --       --     --     --          --  --
ace_disc    --  --     NA     --     --       --     --     --          --  --
            BFGS newuoa
geiger_cts  Y    Y     
ace_cts     Y    Y     
ouch        --   --    
geiger_disc --   --    
ace_disc    --   --    
> 
> smoke1 <- matrix(c(gprptn_cts,aprptn_cts,oprptn,gprptn_disc,aprptn_disc),ncol=13,byrow=TRUE)
> colnames(smoke1) <- c("Liklihood",names(ace)[1:12])
Error in `colnames<-`(`*tmp*`, value = "Liklihood") : 
  length of 'dimnames' [2] not equal to array extent
> rownames(smoke1) <- c("geiger_cts","ace_cts","ouch","geiger_disc","ace_disc")
> aquiprptn <- as.table(smoke1)
> aquiprptn
                  A       B       C       D       E       F       G       H
geiger_cts  -22.666   1.000   1.000   1.000   1.000   1.000   1.000   1.000
ace_cts      48.754   0.000   0.000   0.000   1.000   0.000   1.000   1.000
ouch        -26.097   0.455   0.455   0.500   0.500   0.409   0.455   0.000
geiger_disc -24.518   0.000   0.000   0.000   0.000   0.000   0.727   0.818
ace_disc    -24.769   0.318   0.273   0.000   0.045   0.000   0.273   0.364
                  I       J       K       L       M
geiger_cts    1.000   1.000   1.000   1.000   1.000
ace_cts       0.000   1.000   0.000   1.000   1.000
ouch          0.364   0.409   0.318   0.364   0.545
geiger_disc   0.864   0.545   0.000   0.000   0.955
ace_disc      0.318   0.182   0.182   0.182   0.045
> 
> smoke2 <- matrix(c(p(geiger_cts),p(ace_cts),p(ouch),p(geiger_disc),p(ace_disc),ncol=12,byrow=TRUE)
+ colnames(smoke2) <- names(geiger_cts)[1:12]
Error: unexpected symbol in:
"smoke2 <- matrix(c(p(geiger_cts),p(ace_cts),p(ouch),p(geiger_disc),p(ace_disc),ncol=12,byrow=TRUE)
colnames"
> rownames(smoke2) <- c("geiger_cts","ace_cts","ouch","geiger_disc","ace_disc")
Error in rownames(smoke2) <- c("geiger_cts", "ace_cts", "ouch", "geiger_disc",  : 
  object 'smoke2' not found
> aqui_lik <- as.table(smoke2)
Error in as.table(smoke2) : object 'smoke2' not found
> aqui_lik
Error: object 'aqui_lik' not found
> aprptn_disc
 [1] -24.769   0.318   0.273   0.000   0.045   0.000   0.273   0.364   0.318
[10]   0.182   0.182   0.182   0.045
> gprptn_disc
 [1] -24.518   0.000   0.000   0.000   0.000   0.000   0.727   0.818   0.864
[10]   0.545   0.000   0.000   0.955
> oprptn
 [1] -26.097   0.455   0.455   0.500   0.500   0.409   0.455   0.000   0.364
[10]   0.409   0.318   0.364   0.545
> aprptn_cts
 [1] 48.754  0.000  0.000  0.000  1.000  0.000  1.000  1.000  0.000  1.000
[11]  0.000  1.000  1.000
> gprptn_cts
 [1] -22.666   1.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000
[10]   1.000   1.000   1.000   1.000
> smoke <- matrix(c(gliktf_cts,aliktf_cts,oliktf,gliktf_disc,aliktf_disc),ncol=12,byrow=TRUE)
> colnames(smoke) <- names(geiger_cts)[1:12]
> rownames(smoke) <- c("geiger_cts","ace_cts","ouch","geiger_disc","ace_disc")
> aquilegia <- as.table(smoke)
> aquilegia
            spg Rcgmin Rvmmin bobyqa L-BFGS-B nlminb ucminf Nelder-Mead nlm CG
geiger_cts  Y   Y      Y      Y      Y        Y      Y      Y           Y   Y 
ace_cts     --  --     NA     Y      --       Y      Y      --          Y   O 
ouch        --  --     --     --     --       --     --     --          --  --
geiger_disc --  --     --     O      --       --     --     --          --  --
ace_disc    --  --     NA     --     --       --     --     --          --  --
            BFGS newuoa
geiger_cts  Y    Y     
ace_cts     Y    Y     
ouch        --   --    
geiger_disc --   --    
ace_disc    --   --    
> smoke1 <- matrix(c(gprptn_cts,aprptn_cts,oprptn,gprptn_disc,aprptn_disc),ncol=13,byrow=TRUE)
> colnames(smoke1) <- c("Liklihood",names(ace_cts)[1:12])
> rownames(smoke1) <- c("geiger_cts","ace_cts","ouch","geiger_disc","ace_disc")
> aquiprptn <- as.table(smoke1)
> aquiprptn
            Liklihood     spg  Rcgmin  Rvmmin  bobyqa L-BFGS-B  nlminb  ucminf
geiger_cts    -22.666   1.000   1.000   1.000   1.000    1.000   1.000   1.000
ace_cts        48.754   0.000   0.000   0.000   1.000    0.000   1.000   1.000
ouch          -26.097   0.455   0.455   0.500   0.500    0.409   0.455   0.000
geiger_disc   -24.518   0.000   0.000   0.000   0.000    0.000   0.727   0.818
ace_disc      -24.769   0.318   0.273   0.000   0.045    0.000   0.273   0.364
            Nelder-Mead     nlm      CG    BFGS  newuoa
geiger_cts        1.000   1.000   1.000   1.000   1.000
ace_cts           0.000   1.000   0.000   1.000   1.000
ouch              0.364   0.409   0.318   0.364   0.545
geiger_disc       0.864   0.545   0.000   0.000   0.955
ace_disc          0.318   0.182   0.182   0.182   0.045
> ace_disc
$spg
               I lb    ub          P         L conv        time prptn
 [1,]  0.0100000  0 1e+50  0.1086145 -24.76883    0 0.012566667     1
 [2,]  0.4857143  0 1e+50  0.1086145 -24.76883    0 0.010150000     1
 [3,]  0.9614286  0 1e+50  0.1086145 -24.76883    0 0.010700000     1
 [4,]  1.4371429  0 1e+50  0.1086145 -24.76883    0 0.010616667     1
 [5,]  1.9128571  0 1e+50  0.1086145 -24.76883    0 0.010700000     1
 [6,]  2.3885714  0 1e+50  0.1086146 -24.76883    0 0.011883333     1
 [7,]  2.8642857  0 1e+50  0.1086145 -24.76883    0 0.010800000     1
 [8,]  3.3400000  0 1e+50  3.3400000 -31.85976    0 0.009050000     0
 [9,]  3.8157143  0 1e+50  3.8157143 -31.85976    0 0.004366667     0
[10,]  4.2914286  0 1e+50  4.2914286 -31.85976    0 0.004283333     0
[11,]  4.7671429  0 1e+50  4.7671429 -31.85976    0 0.004133333     0
[12,]  5.2428571  0 1e+50  5.2428571 -31.85976    0 0.004083333     0
[13,]  5.7185714  0 1e+50  5.7185714 -31.85976    0 0.002383333     0
[14,]  6.1942857  0 1e+50  6.1942857 -31.85976    0 0.002150000     0
[15,]  6.6700000  0 1e+50  6.6700000 -31.85976    0 0.002200000     0
[16,]  7.1457143  0 1e+50  7.1457143 -31.85976    0 0.002200000     0
[17,]  7.6214286  0 1e+50  7.6214286 -31.85976    0 0.004050000     0
[18,]  8.0971429  0 1e+50  8.0971429 -31.85976    0 0.003216667     0
[19,]  8.5728571  0 1e+50  8.5728571 -31.85976    0 0.002950000     0
[20,]  9.0485714  0 1e+50  9.0485714 -31.85976    0 0.002900000     0
[21,]  9.5242857  0 1e+50  9.5242857 -31.85976    0 0.001933333     0
[22,] 10.0000000  0 1e+50 10.0000000 -31.85976    0 0.001766667     0

$Rcgmin
               I lb    ub          P         L conv        time prptn
 [1,]  0.0100000  0 1e+50  0.0100000 -31.85976    0 0.003333333     0
 [2,]  0.4857143  0 1e+50  0.1086146 -24.76883    0 0.011633333     1
 [3,]  0.9614286  0 1e+50  0.1086146 -24.76883    0 0.008850000     1
 [4,]  1.4371429  0 1e+50  0.1086146 -24.76883    0 0.015466667     1
 [5,]  1.9128571  0 1e+50  0.1086146 -24.76883    0 0.028366667     1
 [6,]  2.3885714  0 1e+50  0.1086146 -24.76883    0 0.224783333     1
 [7,]  2.8642857  0 1e+50  0.1086146 -24.76883    0 0.080816667     1
 [8,]  3.3400000  0 1e+50  3.2364969 -31.85976    1 0.214866667     0
 [9,]  3.8157143  0 1e+50  3.8157143 -31.85976    0 0.004583333     0
[10,]  4.2914286  0 1e+50  4.2914286 -31.85976    0 0.004916667     0
[11,]  4.7671429  0 1e+50  4.7671429 -31.85976    0 0.004266667     0
[12,]  5.2428571  0 1e+50  5.2428571 -31.85976    0 0.004283333     0
[13,]  5.7185714  0 1e+50  5.7185714 -31.85976    0 0.002483333     0
[14,]  6.1942857  0 1e+50  6.1942857 -31.85976    0 0.002283333     0
[15,]  6.6700000  0 1e+50  6.6700000 -31.85976    0 0.002233333     0
[16,]  7.1457143  0 1e+50  7.1457143 -31.85976    0 0.002333333     0
[17,]  7.6214286  0 1e+50  7.6214286 -31.85976    0 0.003716667     0
[18,]  8.0971429  0 1e+50  8.0971429 -31.85976    0 0.003300000     0
[19,]  8.5728571  0 1e+50  8.5728571 -31.85976    0 0.003083333     0
[20,]  9.0485714  0 1e+50  9.0485714 -31.85976    0 0.003050000     0
[21,]  9.5242857  0 1e+50  9.5242857 -31.85976    0 0.002033333     0
[22,] 10.0000000  0 1e+50 10.0000000 -31.85976    0 0.001850000     0

$Rvmmin
               I    lb           ub          P         L        conv
 [1,]  0.0100000 0e+00 1.000000e+50   0.010000 -34.21554 0.000000000
 [2,]  0.4857143 0e+00 1.000000e+50         NA        NA 0.001633333
 [3,]  0.9614286 0e+00 1.000000e+50         NA        NA 0.001466667
 [4,]  1.4371429 0e+00 1.000000e+50         NA        NA 0.003100000
 [5,]  1.9128571 0e+00 1.000000e+50         NA        NA 0.024233333
 [6,]  2.3885714 0e+00 1.000000e+50         NA        NA 0.146733333
 [7,]  2.8642857 0e+00 1.000000e+50         NA        NA 0.139383333
 [8,]  3.3400000 0e+00 1.000000e+50         NA        NA 0.148183333
 [9,]  3.8157143 0e+00 1.000000e+50   3.815712 -31.85976 0.000000000
[10,]  4.2914286 1e+50 4.291429e+00 -31.859756   0.00000 0.009050000
[11,]  4.7671429 1e+50 4.767143e+00 -31.859756   0.00000 0.008300000
[12,]  5.2428571 1e+50 5.242857e+00 -31.859756   0.00000 0.009633333
[13,]  5.7185714 1e+50 5.718571e+00 -31.859756   0.00000 0.004966667
[14,]  6.1942857 1e+50 6.194286e+00 -31.859756   0.00000 0.003833333
[15,]  6.6700000 1e+50 6.670000e+00 -31.859756   0.00000 0.003833333
[16,]  7.1457143 1e+50 7.145714e+00 -31.859756   0.00000 0.005483333
[17,]  7.6214286 0e+00 1.000000e+50   7.621429 -31.85976 0.000000000
[18,]  8.0971429 1e+50 8.097143e+00 -31.859756   0.00000 0.006783333
[19,]  8.5728571 1e+50 8.572857e+00 -31.859756   0.00000 0.004500000
[20,]  9.0485714 1e+50 9.048571e+00 -31.859756   0.00000 0.005350000
[21,]  9.5242857 1e+50 9.524285e+00 -31.859756   0.00000 0.004733333
[22,] 10.0000000 1e+50 1.000000e+01 -31.859756   0.00000 0.003000000
             time prptn
 [1,] 0.003216667     0
 [2,] 0.000000000     0
 [3,] 0.000000000     0
 [4,] 0.000000000     0
 [5,] 0.000000000     0
 [6,] 0.000000000     0
 [7,] 0.000000000     0
 [8,] 0.000000000     0
 [9,] 0.008100000     0
[10,] 0.000000000     0
[11,] 0.000000000     0
[12,] 0.000000000     0
[13,] 0.000000000     0
[14,] 0.000000000     0
[15,] 0.000000000     0
[16,] 0.000000000     0
[17,] 0.006200000     0
[18,] 0.000000000     0
[19,] 0.000000000     0
[20,] 0.000000000     0
[21,] 0.000000000     0
[22,] 0.000000000     0

$bobyqa
               I lb    ub          P         L conv        time prptn
 [1,]  0.0100000  0 1e+50  0.1086146 -24.76883    0 0.008300000     1
 [2,]  0.4857143  0 1e+50  0.1069759 -24.76941    0 0.005750000     0
 [3,]  0.9614286  0 1e+50  0.1880127 -25.57516    0 0.006433333     0
 [4,]  1.4371429  0 1e+50  0.2810413 -27.23021    0 0.007750000     0
 [5,]  1.9128571  0 1e+50  0.3740699 -28.78720    0 0.007216667     0
 [6,]  2.3885714  0 1e+50  0.4670985 -29.96885    0 0.004350000     0
 [7,]  2.8642857  0 1e+50  0.5601271 -30.75704    0 0.004616667     0
 [8,]  3.3400000  0 1e+50  0.6531557 -31.23827    0 0.004416667     0
 [9,]  3.8157143  0 1e+50  0.7461843 -31.51546    0 0.006800000     0
[10,]  4.2914286  0 1e+50  0.8392129 -31.66996    0 0.006666667     0
[11,]  4.7671429  0 1e+50  0.8510320 -31.68378    0 0.007183333     0
[12,]  5.2428571  0 1e+50  1.3267462 -31.85064    0 0.008783333     0
[13,]  5.7185714  0 1e+50  1.8024605 -31.85920    0 0.004400000     0
[14,]  6.1942857  0 1e+50  2.2781748 -31.85972    0 0.004216667     0
[15,]  6.6700000  0 1e+50  2.7538891 -31.85975    0 0.004550000     0
[16,]  7.1457143  0 1e+50  7.1362143 -31.85976    0 0.003216667     0
[17,]  7.6214286  0 1e+50  1.0136510 -31.79707    0 0.006033333     0
[18,]  8.0971429  0 1e+50  9.9961463 -31.85976    0 0.005466667     0
[19,]  8.5728571  0 1e+50  8.6030960 -31.85976    0 0.005133333     0
[20,]  9.0485714  0 1e+50  8.0985714 -31.85976    0 0.004883333     0
[21,]  9.5242857  0 1e+50  9.6260714 -31.85976    0 0.002933333     0
[22,] 10.0000000  0 1e+50 10.0000000 -31.85976    0 0.002650000     0

$`L-BFGS-B`
               I lb    ub         P              L conv        time prptn
 [1,]  0.0100000  0 1e+50        NA -4.494233e+307 9999 0.008116667     0
 [2,]  0.4857143  0 1e+50        NA -4.494233e+307 9999 0.001800000     0
 [3,]  0.9614286  0 1e+50        NA -4.494233e+307 9999 0.001566667     0
 [4,]  1.4371429  0 1e+50        NA -4.494233e+307 9999 0.002916667     0
 [5,]  1.9128571  0 1e+50        NA -4.494233e+307 9999 0.004316667     0
 [6,]  2.3885714  0 1e+50  2.388348  -3.185974e+01    0 0.003050000     0
 [7,]  2.8642857  0 1e+50  2.864270  -3.185975e+01    0 0.003466667     0
 [8,]  3.3400000  0 1e+50  3.339999  -3.185976e+01    0 0.002933333     0
 [9,]  3.8157143  0 1e+50  3.815714  -3.185976e+01    0 0.008733333     0
[10,]  4.2914286  0 1e+50  4.291429  -3.185976e+01    0 0.015350000     0
[11,]  4.7671429  0 1e+50  4.767143  -3.185976e+01    0 0.004583333     0
[12,]  5.2428571  0 1e+50  5.242857  -3.185976e+01   52 0.012766667     0
[13,]  5.7185714  0 1e+50  5.718571  -3.185976e+01    0 0.002866667     0
[14,]  6.1942857  0 1e+50  6.194286  -3.185976e+01    0 0.002900000     0
[15,]  6.6700000  0 1e+50  6.670000  -3.185976e+01    0 0.005800000     0
[16,]  7.1457143  0 1e+50  7.145714  -3.185976e+01    0 0.002400000     0
[17,]  7.6214286  0 1e+50  7.621429  -3.185976e+01    0 0.005166667     0
[18,]  8.0971429  0 1e+50  8.097143  -3.185976e+01    0 0.003400000     0
[19,]  8.5728571  0 1e+50  8.572857  -3.185976e+01   52 0.009850000     0
[20,]  9.0485714  0 1e+50  9.048571  -3.185976e+01    0 0.003333333     0
[21,]  9.5242857  0 1e+50  9.524286  -3.185976e+01    0 0.004783333     0
[22,] 10.0000000  0 1e+50 10.000000  -3.185976e+01   52 0.004883333     0

$nlminb
               I lb ub          P         L conv        time prptn
 [1,]  0.0100000 NA NA  0.1086146 -24.76883    0 0.008400000     1
 [2,]  0.4857143 NA NA  0.1086146 -24.76883    0 0.006233333     1
 [3,]  0.9614286 NA NA  0.1086146 -24.76883    0 0.006333333     1
 [4,]  1.4371429 NA NA  0.1086146 -24.76883    0 0.007233333     1
 [5,]  1.9128571 NA NA  0.1086146 -24.76883    0 0.009800000     1
 [6,]  2.3885714 NA NA  0.1086147 -24.76883    0 0.005400000     1
 [7,]  2.8642857 NA NA  2.8642699 -31.85975    0 0.003083333     0
 [8,]  3.3400000 NA NA  3.3399987 -31.85976    0 0.002850000     0
 [9,]  3.8157143 NA NA  3.8157137 -31.85976    0 0.005183333     0
[10,]  4.2914286 NA NA  4.2914278 -31.85976    0 0.005200000     0
[11,]  4.7671429 NA NA  4.7671429 -31.85976    0 0.004383333     0
[12,]  5.2428571 NA NA  5.2428571 -31.85976    0 0.004133333     0
[13,]  5.7185714 NA NA  5.7185714 -31.85976    0 0.002533333     0
[14,]  6.1942857 NA NA  6.1942857 -31.85976    0 0.002300000     0
[15,]  6.6700000 NA NA  6.6700001 -31.85976    0 0.002433333     0
[16,]  7.1457143 NA NA  7.1457143 -31.85976    0 0.002266667     0
[17,]  7.6214286 NA NA  7.6214286 -31.85976    0 0.003633333     0
[18,]  8.0971429 NA NA  8.0971429 -31.85976    0 0.003366667     0
[19,]  8.5728571 NA NA  8.5728571 -31.85976    0 0.003750000     0
[20,]  9.0485714 NA NA  9.0485713 -31.85976    0 0.003416667     0
[21,]  9.5242857 NA NA  9.5242857 -31.85976    0 0.002066667     0
[22,] 10.0000000 NA NA 10.0000000 -31.85976    0 0.001933333     0

$ucminf
               I lb ub          P         L conv        time prptn
 [1,]  0.0100000 NA NA  0.1086145 -24.76883    0 0.006533333     1
 [2,]  0.4857143 NA NA  0.1086145 -24.76883    0 0.006566667     1
 [3,]  0.9614286 NA NA  0.1086146 -24.76883    0 0.010833333     1
 [4,]  1.4371429 NA NA  0.1086145 -24.76883    0 0.006750000     1
 [5,]  1.9128571 NA NA  0.1086145 -24.76883    0 0.007383333     1
 [6,]  2.3885714 NA NA  0.1086145 -24.76883    0 0.003833333     1
 [7,]  2.8642857 NA NA  0.1086146 -24.76883    0 0.006266667     1
 [8,]  3.3400000 NA NA  0.1086145 -24.76883    0 0.004383333     1
 [9,]  3.8157143 NA NA  3.8157143 -31.85976    0 0.004366667     0
[10,]  4.2914286 NA NA  4.2914286 -31.85976    0 0.004883333     0
[11,]  4.7671429 NA NA  4.7671429 -31.85976    0 0.004416667     0
[12,]  5.2428571 NA NA  5.2428571 -31.85976    0 0.003950000     0
[13,]  5.7185714 NA NA  5.7185714 -31.85976    0 0.002416667     0
[14,]  6.1942857 NA NA  6.1942857 -31.85976    0 0.002233333     0
[15,]  6.6700000 NA NA  6.6700000 -31.85976    0 0.002166667     0
[16,]  7.1457143 NA NA  7.1457143 -31.85976    0 0.002116667     0
[17,]  7.6214286 NA NA  7.6214286 -31.85976    0 0.003466667     0
[18,]  8.0971429 NA NA  8.0971429 -31.85976    0 0.003166667     0
[19,]  8.5728571 NA NA  8.5728571 -31.85976    0 0.003383333     0
[20,]  9.0485714 NA NA  9.0485714 -31.85976    0 0.002783333     0
[21,]  9.5242857 NA NA  9.5242857 -31.85976    0 0.002033333     0
[22,] 10.0000000 NA NA 10.0000000 -31.85976    0 0.001783333     0

$`Nelder-Mead`
               I lb ub          P         L conv        time prptn
 [1,]  0.0100000 NA NA  0.1086250 -24.76883    0 0.006916667     1
 [2,]  0.4857143 NA NA  0.1086217 -24.76883    0 0.007416667     1
 [3,]  0.9614286 NA NA  0.1086302 -24.76883    0 0.007583333     1
 [4,]  1.4371429 NA NA  0.1086278 -24.76883    0 0.007866667     1
 [5,]  1.9128571 NA NA  0.1086256 -24.76883    0 0.007816667     1
 [6,]  2.3885714 NA NA  0.1086403 -24.76884    0 0.004883333     1
 [7,]  2.8642857 NA NA  0.1085995 -24.76883    0 0.005133333     1
 [8,]  3.3400000 NA NA  3.3400000 -31.85976    0 0.002800000     0
 [9,]  3.8157143 NA NA  3.8157143 -31.85976    0 0.004450000     0
[10,]  4.2914286 NA NA  4.2914286 -31.85976    0 0.004266667     0
[11,]  4.7671429 NA NA  4.7671429 -31.85976    0 0.005316667     0
[12,]  5.2428571 NA NA  5.2428571 -31.85976    0 0.003850000     0
[13,]  5.7185714 NA NA  5.7185714 -31.85976    0 0.002666667     0
[14,]  6.1942857 NA NA  6.1942857 -31.85976    0 0.002166667     0
[15,]  6.6700000 NA NA  6.6700000 -31.85976    0 0.002183333     0
[16,]  7.1457143 NA NA  7.8602857 -31.85976    0 0.002016667     0
[17,]  7.6214286 NA NA  8.3835714 -31.85976    0 0.003300000     0
[18,]  8.0971429 NA NA  8.9068571 -31.85976    0 0.002716667     0
[19,]  8.5728571 NA NA  8.5728571 -31.85976    0 0.003266667     0
[20,]  9.0485714 NA NA  9.0485714 -31.85976    0 0.002933333     0
[21,]  9.5242857 NA NA 10.4767143 -31.85976    0 0.002133333     0
[22,] 10.0000000 NA NA 10.0000000 -31.85976    0 0.001783333     0

$nlm
               I lb ub          P         L conv        time prptn
 [1,]  0.0100000 NA NA 17.3885063 -31.85976    0 0.003800000     0
 [2,]  0.4857143 NA NA  0.1086141 -24.76883    0 0.006816667     1
 [3,]  0.9614286 NA NA  0.1086141 -24.76883    0 0.007316667     1
 [4,]  1.4371429 NA NA  0.1086141 -24.76883    0 0.007783333     1
 [5,]  1.9128571 NA NA  0.1086146 -24.76883    0 0.025100000     1
 [6,]  2.3885714 NA NA  2.3647688 -31.85973    1 0.020600000     0
 [7,]  2.8642857 NA NA  2.8642857 -31.85975    0 0.002683333     0
 [8,]  3.3400000 NA NA  3.3400000 -31.85976    0 0.002850000     0
 [9,]  3.8157143 NA NA  3.8157143 -31.85976    0 0.004333333     0
[10,]  4.2914286 NA NA  4.2914286 -31.85976    0 0.004883333     0
[11,]  4.7671429 NA NA  4.7671429 -31.85976    0 0.005733333     0
[12,]  5.2428571 NA NA  5.2428571 -31.85976    0 0.004083333     0
[13,]  5.7185714 NA NA  5.7185714 -31.85976    0 0.002900000     0
[14,]  6.1942857 NA NA  6.1942857 -31.85976    0 0.002233333     0
[15,]  6.6700000 NA NA  6.6700000 -31.85976    0 0.002300000     0
[16,]  7.1457143 NA NA  7.1457143 -31.85976    0 0.002216667     0
[17,]  7.6214286 NA NA  7.6214286 -31.85976    0 0.003566667     0
[18,]  8.0971429 NA NA  8.0971429 -31.85976    0 0.002283333     0
[19,]  8.5728571 NA NA  8.5728571 -31.85976    0 0.003450000     0
[20,]  9.0485714 NA NA  9.0485714 -31.85976    0 0.002983333     0
[21,]  9.5242857 NA NA  9.5242857 -31.85976    0 0.002150000     0
[22,] 10.0000000 NA NA 10.0000000 -31.85976    0 0.001816667     0

$CG
               I lb ub          P         L conv        time prptn
 [1,]  0.0100000 NA NA  4.9680254 -31.85976    0 0.004733333     0
 [2,]  0.4857143 NA NA  0.1086146 -24.76883    0 0.026666667     1
 [3,]  0.9614286 NA NA  0.1086146 -24.76883    0 0.024483333     1
 [4,]  1.4371429 NA NA  0.1086149 -24.76883    0 0.032416667     1
 [5,]  1.9128571 NA NA  0.1086146 -24.76883    0 0.052683333     1
 [6,]  2.3885714 NA NA  2.3650231 -31.85973    1 0.038233333     0
 [7,]  2.8642857 NA NA  2.8627171 -31.85975    1 0.040900000     0
 [8,]  3.3400000 NA NA  3.3400000 -31.85976    0 0.002783333     0
 [9,]  3.8157143 NA NA  3.8157143 -31.85976    0 0.004333333     0
[10,]  4.2914286 NA NA  4.2914286 -31.85976    0 0.004916667     0
[11,]  4.7671429 NA NA  4.7671429 -31.85976    0 0.004083333     0
[12,]  5.2428571 NA NA  5.2428571 -31.85976    0 0.003983333     0
[13,]  5.7185714 NA NA  5.7185714 -31.85976    0 0.002483333     0
[14,]  6.1942857 NA NA  6.1942857 -31.85976    0 0.002216667     0
[15,]  6.6700000 NA NA  6.6700000 -31.85976    0 0.002250000     0
[16,]  7.1457143 NA NA  7.1457143 -31.85976    0 0.002200000     0
[17,]  7.6214286 NA NA  7.6214286 -31.85976    0 0.003683333     0
[18,]  8.0971429 NA NA  8.0971429 -31.85976    0 0.003150000     0
[19,]  8.5728571 NA NA  8.5728571 -31.85976    0 0.003383333     0
[20,]  9.0485714 NA NA  9.0485714 -31.85976    0 0.002966667     0
[21,]  9.5242857 NA NA  9.5242857 -31.85976    0 0.002000000     0
[22,] 10.0000000 NA NA 10.0000000 -31.85976    0 0.001866667     0

$BFGS
               I lb ub          P         L conv        time prptn
 [1,]  0.0100000 NA NA  9.9260509 -31.85976    0 0.004616667     0
 [2,]  0.4857143 NA NA  0.1086173 -24.76883    0 0.009316667     1
 [3,]  0.9614286 NA NA  0.1086146 -24.76883    0 0.008350000     1
 [4,]  1.4371429 NA NA  0.1086134 -24.76883    0 0.012666667     1
 [5,]  1.9128571 NA NA  0.1086096 -24.76883    0 0.032350000     1
 [6,]  2.3885714 NA NA  2.3883485 -31.85974    0 0.002733333     0
 [7,]  2.8642857 NA NA  2.8642699 -31.85975    0 0.002966667     0
 [8,]  3.3400000 NA NA  3.3399988 -31.85976    0 0.002883333     0
 [9,]  3.8157143 NA NA  3.8157143 -31.85976    0 0.004783333     0
[10,]  4.2914286 NA NA  4.2914286 -31.85976    0 0.005416667     0
[11,]  4.7671429 NA NA  4.7671429 -31.85976    0 0.004066667     0
[12,]  5.2428571 NA NA  5.2428571 -31.85976    0 0.004966667     0
[13,]  5.7185714 NA NA  5.7185714 -31.85976    0 0.002716667     0
[14,]  6.1942857 NA NA  6.1942857 -31.85976    0 0.002483333     0
[15,]  6.6700000 NA NA  6.6700000 -31.85976    0 0.002433333     0
[16,]  7.1457143 NA NA  7.1457143 -31.85976    0 0.002183333     0
[17,]  7.6214286 NA NA  7.6214286 -31.85976    0 0.003516667     0
[18,]  8.0971429 NA NA  8.0971429 -31.85976    0 0.003466667     0
[19,]  8.5728571 NA NA  8.5728571 -31.85976    0 0.003933333     0
[20,]  9.0485714 NA NA  9.0485714 -31.85976    0 0.002983333     0
[21,]  9.5242857 NA NA  9.5242857 -31.85976    0 0.002283333     0
[22,] 10.0000000 NA NA 10.0000000 -31.85976    0 0.002150000     0

$newuoa
               I lb ub          P              L conv        time prptn
 [1,]  0.0100000 NA NA  0.1086146  -2.476883e+01    0 0.006683333     1
 [2,]  0.4857143 NA NA         NA -4.494233e+307 9999 0.002533333     0
 [3,]  0.9614286 NA NA         NA -4.494233e+307 9999 0.002516667     0
 [4,]  1.4371429 NA NA         NA -4.494233e+307 9999 0.002333333     0
 [5,]  1.9128571 NA NA         NA -4.494233e+307 9999 0.001366667     0
 [6,]  2.3885714 NA NA         NA -4.494233e+307 9999 0.001283333     0
 [7,]  2.8642857 NA NA         NA -4.494233e+307 9999 0.001433333     0
 [8,]  3.3400000 NA NA         NA -4.494233e+307 9999 0.001383333     0
 [9,]  3.8157143 NA NA         NA -4.494233e+307 9999 0.002050000     0
[10,]  4.2914286 NA NA         NA -4.494233e+307 9999 0.001983333     0
[11,]  4.7671429 NA NA         NA -4.494233e+307 9999 0.001883333     0
[12,]  5.2428571 NA NA         NA -4.494233e+307 9999 0.001916667     0
[13,]  5.7185714 NA NA         NA -4.494233e+307 9999 0.001383333     0
[14,]  6.1942857 NA NA         NA -4.494233e+307 9999 0.001200000     0
[15,]  6.6700000 NA NA         NA -4.494233e+307 9999 0.001316667     0
[16,]  7.1457143 NA NA  7.1298810  -3.185976e+01    0 0.003233333     0
[17,]  7.6214286 NA NA         NA -4.494233e+307 9999 0.001966667     0
[18,]  8.0971429 NA NA 10.0921429  -3.185976e+01    0 0.004633333     0
[19,]  8.5728571 NA NA  8.6203571  -3.185976e+01    0 0.005716667     0
[20,]  9.0485714 NA NA  8.0985714  -3.185976e+01    0 0.004566667     0
[21,]  9.5242857 NA NA  9.4283357  -3.185976e+01    0 0.003150000     0
[22,] 10.0000000 NA NA 10.0000000  -3.185976e+01    0 0.002733333     0

$mean_lik
           spg         Rcgmin         Rvmmin         bobyqa       L-BFGS-B 
 -2.960360e+01  -2.992590e+01             NA  -3.037950e+01 -1.021417e+307 
        nlminb         ucminf    Nelder-Mead            nlm             CG 
 -2.992590e+01  -2.928120e+01  -2.960360e+01  -3.057050e+01  -3.057050e+01 
          BFGS         newuoa 
 -3.057050e+01 -3.064250e+307 

$MLE
        spg      Rcgmin      Rvmmin      bobyqa    L-BFGS-B      nlminb 
   -24.7688    -24.7688          NA    -24.7688    -31.8597    -24.7688 
     ucminf Nelder-Mead         nlm          CG        BFGS      newuoa 
   -24.7688    -24.7688    -24.7688    -24.7688    -24.7688    -24.7688 

> smoke <- matrix(c(gliktf_cts,aliktf_cts,oliktf,gliktf_disc,aliktf_disc),ncol=12,byrow=TRUE)
> colnames(smoke) <- names(geiger_cts)[1:12]
> rownames(smoke) <- c("geiger_cts","ace_cts","ouch","geiger_disc","ace_disc")
> aquilegia <- as.table(smoke)
> aquilegia
            spg Rcgmin Rvmmin bobyqa L-BFGS-B nlminb ucminf Nelder-Mead nlm CG
geiger_cts  Y   Y      Y      Y      Y        Y      Y      Y           Y   Y 
ace_cts     --  --     NA     Y      --       Y      Y      --          Y   O 
ouch        --  --     --     --     --       --     --     --          --  --
geiger_disc --  --     --     O      --       --     --     --          --  --
ace_disc    --  --     NA     --     --       --     --     --          --  --
            BFGS newuoa
geiger_cts  Y    Y     
ace_cts     Y    Y     
ouch        --   --    
geiger_disc --   --    
ace_disc    --   --    
> 
> smoke1 <- matrix(c(gprptn_cts,aprptn_cts,oprptn,gprptn_disc,aprptn_disc),ncol=13,byrow=TRUE)
> colnames(smoke1) <- c("Liklihood",names(ace_cts)[1:12])
> rownames(smoke1) <- c("geiger_cts","ace_cts","ouch","geiger_disc","ace_disc")
> aquiprptn <- as.table(smoke1)
> aquiprptn
            Liklihood     spg  Rcgmin  Rvmmin  bobyqa L-BFGS-B  nlminb  ucminf
geiger_cts    -22.666   1.000   1.000   1.000   1.000    1.000   1.000   1.000
ace_cts        48.754   0.000   0.000   0.000   1.000    0.000   1.000   1.000
ouch          -26.097   0.455   0.455   0.500   0.500    0.409   0.455   0.000
geiger_disc   -24.518   0.000   0.000   0.000   0.000    0.000   0.727   0.818
ace_disc      -24.769   0.318   0.273   0.000   0.045    0.000   0.273   0.364
            Nelder-Mead     nlm      CG    BFGS  newuoa
geiger_cts        1.000   1.000   1.000   1.000   1.000
ace_cts           0.000   1.000   0.000   1.000   1.000
ouch              0.364   0.409   0.318   0.364   0.545
geiger_disc       0.864   0.545   0.000   0.000   0.955
ace_disc          0.318   0.182   0.182   0.182   0.045
> 
> smoke2 <- matrix(c(p(geiger_cts),p(ace_cts),p(ouch),p(geiger_disc),p(ace_disc),ncol=12,byrow=TRUE)
+ colnames(smoke2) <- names(geiger_cts)[1:12]
Error: unexpected symbol in:
"smoke2 <- matrix(c(p(geiger_cts),p(ace_cts),p(ouch),p(geiger_disc),p(ace_disc),ncol=12,byrow=TRUE)
colnames"
> rownames(smoke2) <- c("geiger_cts","ace_cts","ouch","geiger_disc","ace_disc")
Error in rownames(smoke2) <- c("geiger_cts", "ace_cts", "ouch", "geiger_disc",  : 
  object 'smoke2' not found
> aqui_lik <- as.table(smoke2)
Error in as.table(smoke2) : object 'smoke2' not found
> aqui_lik
Error: object 'aqui_lik' not found
> p(ace_disc)
 [1] "-24.7688" "-24.7688" "      NA" "-24.7688" "-31.8597" "-24.7688"
 [7] "-24.7688" "-24.7688" "-24.7688" "-24.7688" "-24.7688" "-24.7688"
> p(geiger_disc)
 [1] "-24.5249" "-24.5249" "-24.5249" "-24.5249" "-24.5249" "-24.5176"
 [7] "-24.5176" "-24.5176" "-24.5176" "-24.5458" "-24.5178" "-24.5176"
> smoke2 <- matrix(c(p(geiger_cts),p(ace_cts),p(ouch),p(geiger_disc),p(ace_disc),ncol=12,byrow=TRUE)
+ )
> smoke2 <- matrix(c(p(geiger_cts),p(ace_cts),p(ouch),p(geiger_disc),p(ace_disc)),ncol=12,byrow=TRUE)
> colnames(smoke2) <- names(geiger_cts)[1:12]
> rownames(smoke2) <- c("geiger_cts","ace_cts","ouch","geiger_disc","ace_disc")
> aqui_lik <- as.table(smoke2)
> aqui_lik
            spg         Rcgmin      Rvmmin      bobyqa      L-BFGS-B   
geiger_cts  -22.6662    -22.6662    -22.6662    -22.6662    -22.6662   
ace_cts      -4.507e+00 -4.494e+307          NA   4.875e+01 -4.494e+307
ouch        -26.0968    -26.0968    -26.0968    -26.0968    -26.0968   
geiger_disc -24.5249    -24.5249    -24.5249    -24.5249    -24.5249   
ace_disc    -24.7688    -24.7688          NA    -24.7688    -31.8597   
            nlminb      ucminf      Nelder-Mead nlm         CG         
geiger_cts  -22.6662    -22.6662    -22.6662    -22.6662    -22.6662   
ace_cts       4.875e+01   4.875e+01   2.876e+01   4.875e+01   4.874e+01
ouch        -26.0968    -26.0968    -26.0968    -26.0968    -26.0968   
geiger_disc -24.5176    -24.5176    -24.5176    -24.5176    -24.5458   
ace_disc    -24.7688    -24.7688    -24.7688    -24.7688    -24.7688   
            BFGS        newuoa     
geiger_cts  -22.6662    -22.6662   
ace_cts       4.875e+01   4.875e+01
ouch        -26.0968    -26.0968   
geiger_disc -24.5178    -24.5176   
ace_disc    -24.7688    -24.7688   
> smoke <- matrix(c(gliktf_cts,aliktf_cts,oliktf,gliktf_disc,aliktf_disc),ncol=12,byrow=TRUE)
> colnames(smoke) <- names(geiger_cts)[1:12]
> rownames(smoke) <- c("geiger_cts","ace_cts","ouch","geiger_disc","ace_disc")
> aquilegia <- as.table(smoke)
> aquilegia
            spg Rcgmin Rvmmin bobyqa L-BFGS-B nlminb ucminf Nelder-Mead nlm CG
geiger_cts  Y   Y      Y      Y      Y        Y      Y      Y           Y   Y 
ace_cts     --  --     NA     Y      --       Y      Y      --          Y   O 
ouch        --  --     --     --     --       --     --     --          --  --
geiger_disc --  --     --     O      --       --     --     --          --  --
ace_disc    --  --     NA     --     --       --     --     --          --  --
            BFGS newuoa
geiger_cts  Y    Y     
ace_cts     Y    Y     
ouch        --   --    
geiger_disc --   --    
ace_disc    --   --    
> 
> smoke1 <- matrix(c(gprptn_cts,aprptn_cts,oprptn,gprptn_disc,aprptn_disc),ncol=13,byrow=TRUE)
> colnames(smoke1) <- c("Liklihood",names(ace_cts)[1:12])
> rownames(smoke1) <- c("geiger_cts","ace_cts","ouch","geiger_disc","ace_disc")
> aquiprptn <- as.table(smoke1)
> aquiprptn
            Liklihood     spg  Rcgmin  Rvmmin  bobyqa L-BFGS-B  nlminb  ucminf
geiger_cts    -22.666   1.000   1.000   1.000   1.000    1.000   1.000   1.000
ace_cts        48.754   0.000   0.000   0.000   1.000    0.000   1.000   1.000
ouch          -26.097   0.455   0.455   0.500   0.500    0.409   0.455   0.000
geiger_disc   -24.518   0.000   0.000   0.000   0.000    0.000   0.727   0.818
ace_disc      -24.769   0.318   0.273   0.000   0.045    0.000   0.273   0.364
            Nelder-Mead     nlm      CG    BFGS  newuoa
geiger_cts        1.000   1.000   1.000   1.000   1.000
ace_cts           0.000   1.000   0.000   1.000   1.000
ouch              0.364   0.409   0.318   0.364   0.545
geiger_disc       0.864   0.545   0.000   0.000   0.955
ace_disc          0.318   0.182   0.182   0.182   0.045
> 
> smoke2 <- matrix(c(p(geiger_cts),p(ace_cts),p(ouch),p(geiger_disc),p(ace_disc)),ncol=12,byrow=TRUE)
> colnames(smoke2) <- names(geiger_cts)[1:12]
> rownames(smoke2) <- c("geiger_cts","ace_cts","ouch","geiger_disc","ace_disc")
> aqui_lik <- as.table(smoke2)
> aqui_lik
            spg         Rcgmin      Rvmmin      bobyqa      L-BFGS-B   
geiger_cts  -22.6662    -22.6662    -22.6662    -22.6662    -22.6662   
ace_cts      -4.507e+00 -4.494e+307          NA   4.875e+01 -4.494e+307
ouch        -26.0968    -26.0968    -26.0968    -26.0968    -26.0968   
geiger_disc -24.5249    -24.5249    -24.5249    -24.5249    -24.5249   
ace_disc    -24.7688    -24.7688          NA    -24.7688    -31.8597   
            nlminb      ucminf      Nelder-Mead nlm         CG         
geiger_cts  -22.6662    -22.6662    -22.6662    -22.6662    -22.6662   
ace_cts       4.875e+01   4.875e+01   2.876e+01   4.875e+01   4.874e+01
ouch        -26.0968    -26.0968    -26.0968    -26.0968    -26.0968   
geiger_disc -24.5176    -24.5176    -24.5176    -24.5176    -24.5458   
ace_disc    -24.7688    -24.7688    -24.7688    -24.7688    -24.7688   
            BFGS        newuoa     
geiger_cts  -22.6662    -22.6662   
ace_cts       4.875e+01   4.875e+01
ouch        -26.0968    -26.0968   
geiger_disc -24.5178    -24.5176   
ace_disc    -24.7688    -24.7688   
> 