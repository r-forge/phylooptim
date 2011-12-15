#Results for geospiza

load("/home/michels/repository/phylooptim/pkg/R/geiger/geigererror.RData")
geiger <- l
geiger$mean_lik <- round(c(mean(geiger[[1]][,2]),mean(geiger[[2]][,2]),mean(geiger[[3]][,2]),mean(geiger[[4]][,2]),mean(geiger[[5]][,2]),mean(geiger[[6]][,2]),mean(geiger[[7]][,2]),mean(geiger[[8]][,2]),mean(geiger[[9]][,2]),mean(geiger[[10]][,2]),mean(geiger[[11]][,2]),mean(geiger[[12]][,2])),5)
names(geiger$mean_lik) <- names(geiger)[1:12]
                 
load("/home/michels/repository/phylooptim/pkg/R/ace/aceerror.RData")
ace <- l
ace$mean_lik <- round(c(mean(ace[[1]][,5]),mean(ace[[2]][,5]),mean(ace[[3]][,5]),mean(ace[[4]][,5]),mean(ace[[5]][,5]),mean(ace[[6]][,5]),mean(ace[[7]][,5]),mean(ace[[8]][,5]),mean(ace[[9]][,5]),mean(ace[[10]][,5]),mean(ace[[11]][,5]),mean(ace[[12]][,5])),5)
names(ace$mean_lik) <- names(ace)[1:12]

load("/home/michels/repository/phylooptim/pkg/R/ouch/oucherror.RData")
ouch <- l
ouch$mean_lik <- round(c(mean(ouch[[1]][,5]),mean(ouch[[2]][,5]),mean(ouch[[3]][,5]),mean(ouch[[4]][,5]),mean(ouch[[5]][,5]),mean(ouch[[6]][,5]),mean(ouch[[7]][,5]),mean(ouch[[8]][,5]),mean(ouch[[9]][,5]),mean(ouch[[10]][,5]),mean(ouch[[11]][,5]),mean(ouch[[12]][,5])),5)
names(ouch$mean_lik) <- names(ouch)[1:12]

#Find the mode
f <- function(x){
xt <- table(x)
t <- as.numeric(names(xt[xt == max(xt)]))  
return(t)}

#True Liklihood
tglik <- f(as.numeric(geiger$mean_lik))
gliktf <- rep(NA,length(geiger$mean_lik))
for (i in c(1:length(gliktf))){if (as.numeric(geiger$mean_lik[i])==tglik){gliktf[i] <- "Y"}else{if (tglik-.05 < as.numeric(geiger$mean_lik[i]) && tglik+.05 > as.numeric(geiger$mean_lik[i])){gliktf[i] <- "O"}else{gliktf[i] <- "--"}}}

talik <- f(as.numeric(ace$mean_lik))
aliktf <- rep(NA,length(ace$mean_lik))
for (i in c(1:length(aliktf))){if (as.numeric(ace$mean_lik[i])==talik){aliktf[i] <- "Y"}else{if (talik-.05 < as.numeric(ace$mean_lik[i]) && talik+.05 > as.numeric(ace$mean_lik[i])){aliktf[i] <- "O"}else{aliktf[i] <- "--"}}}

tolik <- f(as.numeric(ouch$mean_lik))
oliktf <- rep(NA,length(ouch$mean_lik))
for (i in c(1:length(oliktf))){if (as.numeric(ouch$mean_lik[i])==tolik){oliktf[i] <- "Y"}else{if (tolik-.05 < as.numeric(ouch$mean_lik[i]) && tolik+.05 > as.numeric(ouch$mean_lik[i])){oliktf[i] <- "O"}else{oliktf[i] <- "--"}}}

smoke <- matrix(c(gliktf,aliktf,oliktf),ncol=12,byrow=TRUE)
colnames(smoke) <- names(geiger)[1:12]
rownames(smoke) <- c("geiger","ace","ouch")
geospiza <- as.table(smoke)

#Results for aquilegia data
load("/home/michels/repository/phylooptim/pkg/R/geiger/aquigeigererror.RData")
geiger <- l
geiger$mean_lik <- round(c(mean(geiger[[1]][,2]),mean(geiger[[2]][,2]),mean(geiger[[3]][,2]),mean(geiger[[4]][,2]),mean(geiger[[5]][,2]),mean(geiger[[6]][,2]),mean(geiger[[7]][,2]),mean(geiger[[8]][,2]),mean(geiger[[9]][,2]),mean(geiger[[10]][,2]),mean(geiger[[11]][,2]),mean(geiger[[12]][,2])),5)
names(geiger$mean_lik) <- names(geiger)[1:12]
                 
load("/home/michels/repository/phylooptim/pkg/R/ace/aquiaceerror.RData")
ace <- l
ace$mean_lik <- round(c(mean(ace[[1]][,5]),mean(ace[[2]][,5]),mean(ace[[3]][,5]),mean(ace[[4]][,5]),mean(ace[[5]][,5]),mean(ace[[6]][,5]),mean(ace[[7]][,5]),mean(ace[[8]][,5]),mean(ace[[9]][,5]),mean(ace[[10]][,5]),mean(ace[[11]][,5]),mean(ace[[12]][,5])),5)
names(ace$mean_lik) <- names(ace)[1:12]

load("/home/michels/repository/phylooptim/pkg/R/ouch/aquioucherror.RData")
ouch <- l
ouch$mean_lik <- round(c(mean(ouch[[1]][,5]),mean(ouch[[2]][,5]),mean(ouch[[3]][,5]),mean(ouch[[4]][,5]),mean(ouch[[5]][,5]),mean(ouch[[6]][,5]),mean(ouch[[7]][,5]),mean(ouch[[8]][,5]),mean(ouch[[9]][,5]),mean(ouch[[10]][,5]),mean(ouch[[11]][,5]),mean(ouch[[12]][,5])),5)
names(ouch$mean_lik) <- names(ouch)[1:12]

#Find the mode
f <- function(x){
xt <- table(x)
t <- as.numeric(names(xt[xt == max(xt)]))  
return(t)}

#True Liklihood
tglik <- f(as.numeric(geiger$mean_lik))
gliktf <- rep(NA,length(geiger$mean_lik))
for (i in c(1:length(gliktf))){if (as.numeric(geiger$mean_lik[i])==tglik){gliktf[i] <- "Y"}else{if (tglik-.05 < as.numeric(geiger$mean_lik[i]) && tglik+.05 > as.numeric(geiger$mean_lik[i])){gliktf[i] <- "O"}else{gliktf[i] <- "--"}}}

talik <- f(as.numeric(ace$mean_lik))
aliktf <- rep(NA,length(ace$mean_lik))
for (i in c(1:length(aliktf))){if (as.numeric(ace$mean_lik[i])==talik){aliktf[i] <- "Y"}else{if (talik-.05 < as.numeric(ace$mean_lik[i]) && talik+.05 > as.numeric(ace$mean_lik[i])){aliktf[i] <- "O"}else{aliktf[i] <- "--"}}}

tolik <- f(as.numeric(ouch$mean_lik))
oliktf <- rep(NA,length(ouch$mean_lik))
for (i in c(1:length(oliktf))){if (as.numeric(ouch$mean_lik[i])==tolik){oliktf[i] <- "Y"}else{if (tolik-.05 < as.numeric(ouch$mean_lik[i]) && tolik+.05 > as.numeric(ouch$mean_lik[i])){oliktf[i] <- "O"}else{oliktf[i] <- "--"}}}

smoke <- matrix(c(gliktf,aliktf,oliktf),ncol=12,byrow=TRUE)
colnames(smoke) <- names(geiger)[1:12]
rownames(smoke) <- c("geiger","ace","ouch")
aquilegia <- as.table(smoke)
