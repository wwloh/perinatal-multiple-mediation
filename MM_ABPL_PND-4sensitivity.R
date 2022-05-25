rm(list=ls())
source("MM_ABPL_PND-1prep.R")

# fitted regression models for sensitivity analysis using MC estimation #######
# logistic regression model for binary Y, given A, M1, M2
fitY <- glm(as.formula(paste0("Y~",paste(paste0("A*M1*M2*",l_names),collapse="+"))),
            family = binomial("logit"), data = Data_agg, weights= WGT)

# logistic regression model for binary M2, given M1
fitM2 <- glm(as.formula(paste0("M2~",paste(paste0("A*M1*",l_names),collapse="+"))),
             family = binomial("logit"), data = Data_agg, weights= WGT)

# logistic regression model for binary M1
fitM1 <- glm(as.formula(paste0("M1~",paste(paste0("A*",l_names),collapse="+"))),
             family = binomial("logit"), data = Data_agg, weights= WGT)

# sensitivity analysis for cross-world correlation ############################
new.Data_agg <- list()
for (a in 0:1) {
  new.Data_agg.a <- data.table(Data_agg)
  new.Data_agg.a[,"A" := a] # set hypothetical exposure level
  new.Data_agg.a[,paste0("pM1a",a,"_eq0") := 1-predict.glm(
    fitM1,newdata=new.Data_agg.a,type="response")] # Pr(M1(a)=0|L)
  new.Data_agg.a[,"A" := NULL]
  setkey(new.Data_agg.a)
  new.Data_agg[[paste0("a",a)]] <- new.Data_agg.a
  rm(new.Data_agg.a)
}
new.Data_agg <- cbind(new.Data_agg$a0,"pM1a1_eq0"=new.Data_agg$a1[,pM1a1_eq0])
setkey(new.Data_agg)

if(FALSE) {
  # check predicted counterfactual mediators
  ## prob. of SGA=0 should be greater for PA=0
  plot(new.Data_agg[,list(pM1a0_eq0,pM1a1_eq0)],
       cex=new.Data_agg[,WGT]/new.Data_agg[,max(WGT)],
       xlab="Pr(M1(a=0)=0|L)",ylab="Pr(M1(a=1)=0|L)",
       xlim=c(0,1),ylim=c(0,1))
  abline(a=0,b=1,lty=2)
}


if(!file.exists("MM_ABPL_PND-4sensitivity.Rdata")) {
  
  # only needs to be calculated once for the observed dataset
  source("cross-world-correlation-sensitivity.R")
  (rho.grid <- (0:20)/20)
  prM1a3_eq1.M1a1.rho <- list()
  for (rr in 1:length(rho.grid)) {
    ptm=proc.time()[3]
    
    # pr(M1(a3=1)=1|M1(a1=0)=m)
    prM1a3_eq1.M1a1.0 <- t(apply(new.Data_agg, 1, function(x) {
      pM1.i <- as.numeric(x[c("pM1a0_eq0","pM1a1_eq0")])
      names(pM1.i) <- c("a1","a3")
      CrossWorldCorrelation(pZ1_0=pM1.i,rho=rho.grid[rr])
    }))
    
    # pr(M1(a3=0)=1|M1(a1=1)=m)
    prM1a3_eq1.M1a1.1 <- t(apply(new.Data_agg, 1, function(x) {
      pM1.i <- as.numeric(x[c("pM1a1_eq0","pM1a0_eq0")])
      names(pM1.i) <- c("a1","a3")
      CrossWorldCorrelation(pZ1_0=pM1.i,rho=rho.grid[rr])
    }))
    
    prM1a3_eq1.M1a1.rho[[rr]] <- list("a1.0"=prM1a3_eq1.M1a1.0,
                                      "a1.1"=prM1a3_eq1.M1a1.1)
    rm(prM1a3_eq1.M1a1.0,prM1a3_eq1.M1a1.1)
    
    cat("rho=",rho.grid[rr],"time (mins)=",round((proc.time()[3]-ptm)/60),"\n")
  }
  save.image("MM_ABPL_PND-4sensitivity.Rdata")
  q()
}

load("MM_ABPL_PND-4sensitivity.Rdata")
library("data.table")

OneEst <- function(a0,a1,a2,a3,rho) {
  # sample M1(a1)
  new.Data_agg.a <- data.table(Data_agg)
  new.Data_agg.a[,"A" := a1] # set hypothetical exposure level
  pM1a1 <- predict.glm(fitM1,newdata=new.Data_agg.a,type="response")
  new.Data_agg.a[,"M1a1" := rbinom(length(pM1a1),1,pM1a1)]
  rm(pM1a1)
  
  # sample M1(a3) given M1(a1)
  if (a3 == a1) {
    new.Data_agg.a[,"M1a3" := M1a1]
  } else {
    prM1a3_eq1.M1a1 <- prM1a3_eq1.M1a1.rho[[which(rho==rho.grid)]]
    prM1a3_eq1.M1a1 <- prM1a3_eq1.M1a1[[paste0("a1.",a1)]]
    new.Data_agg.a[,"M1a3" := sapply(1:nrow(new.Data_agg.a), function(i) {
      m1a1 <- new.Data_agg.a$M1a1[i]
      rbinom(1,1,prM1a3_eq1.M1a1[i,m1a1+1])
    })]
    rm(prM1a3_eq1.M1a1)
  }
  
  # sample M2(a2,M1(a3))
  new.Data_agg.a[,"A" := a2] # set hypothetical exposure level
  new.Data_agg.a[,"M1" := M1a3] # set counterfactual mediator
  pM2a2 <- predict.glm(fitM2,newdata=new.Data_agg.a,type="response")
  new.Data_agg.a[,"M2a2" := rbinom(length(pM2a2),1,pM2a2)]
  rm(pM2a2)
  
  # sample Y(a0,M1(a1),M2(a2,M1(a3)))
  new.Data_agg.a[,"A" := a0] # set hypothetical exposure level
  new.Data_agg.a[,"M1" := M1a1] # set counterfactual mediator
  new.Data_agg.a[,"M2" := M2a2] # set counterfactual mediator
  pY1.a <- predict.glm(fitY,newdata=new.Data_agg.a,type="response")
  new.Data_agg.a[,"Y.po" := rbinom(length(pY1.a),1,pY1.a)]
  rm(pY1.a)
  
  return(new.Data_agg.a[,list(WGT,Y.po)])
}

all.combis <- expand.grid(a0=0:1,a1=0:1,a2=0:1,a3=0:1,rho=rho.grid)
nrow(all.combis)

# initialize for parallel MC jobs
args <- 13
if (!grepl("apple",sessionInfo()[[1]]$platform)) {
  args <- commandArgs(trailingOnly=TRUE) # for CMD BATCH '--args 1'
  all.combis <- all.combis[rep(1:nrow(all.combis),each=20),]
  nrow(all.combis)
}
(seed <- as.integer(args[1]))
set.seed(seed)

(all.combis[seed,])

ptm=proc.time()[3]
OneMCEst <- rbindlist(lapply(1:1000, function(x) 
  OneEst(a0=all.combis[seed,"a0"],
         a1=all.combis[seed,"a1"],
         a2=all.combis[seed,"a2"],
         a3=all.combis[seed,"a3"],
         rho=all.combis[seed,"rho"])
))
proc.time()[3]-ptm
# 4h per 1000 sims

all.combis[seed,"Y.po"] <- weighted.mean(x=OneMCEst$Y.po,w=OneMCEst$WGT)
res <- all.combis[seed,]
subfolder <- "MM_ABPL_PND-4sensitivity/"
save(res, file=paste0(subfolder,"MM_ABPL_PND-4sensitivity-",seed,".Rdata"))

q()

# results #####################################################################
rm(list=ls())
library("data.table")
load("MM_ABPL_PND-4sensitivity.Rdata")
res.all <- res.mc <- NULL
subfolder <- "MM_ABPL_PND-4sensitivity/"
files_to_load <- list.files(subfolder)[grepl(".Rdata",list.files(subfolder))]
for (idx in 1:length(files_to_load)) {
  # load results
  load(file=paste0(subfolder,files_to_load[idx]))
  res.all[[idx]] <- res
  rm(res)
}
res.all <- rbindlist(res.all, fill=TRUE)
setkey(res.all)

res <- res.all[,mean(Y.po),by=eval(key(res.all)[-ncol(res.all)])]; rm(res.all)
setnames(res,"V1","Y.po")
setkey(res)

OneDecomposition <- function(de.a123,ie1.a023,ie2.a013,ie12.a012,rho.fixed,res) {
  # hypothetical exposure levels
  ## fixed rho only relevant for IE1 only and IE12
  res.rho <- res[rho==rho.fixed] 
  # risk ratios of marginal potential outcomes
  res.est <- NULL
  res.est[["de"]] <- res[a1==de.a123[1] & a2==de.a123[2] & a3==de.a123[3]]
  res.est[["ie1"]] <- res.rho[a0==ie1.a023[1] & a2==ie1.a023[2] & a3==ie1.a023[3]]
  res.est[["ie2"]] <- res[a0==ie2.a013[1] & a1==ie2.a013[2] & a3==ie2.a013[3]]
  res.est[["ie12"]] <- res.rho[a0==ie12.a012[1] & a1==ie12.a012[2] & a2==ie12.a012[3]]
  if (all(ie1.a023[1:2]==ie12.a012[c(1,3)])) {
    # combined IE1 and IE12 if a0 and a1 are the same
    ## can be guaranteed by only decomposing combined IE1 as below
    res.est[["ie1.c"]] <- res[a0==ie1.a023[1] & a2==ie1.a023[2]]
  }
  res.est <- lapply(res.est, function(one.est) {
    if (nrow(one.est)>2) {
      # average over rho for effects that do not require rho
      one.est <- one.est[,mean(Y.po),by=list(a0,a1,a2,a3)]
      setnames(one.est,"V1","Y.po")
      setkey(one.est)
    }
    a.sums <- rowSums(one.est[,list(a0,a1,a2,a3)])
    one.est[which.max(a.sums),Y.po]/one.est[which.min(a.sums),Y.po]  
  })
  return( unlist(res.est) )
}

# 6 possible decompositions when using NE models: DE, IE1, IE2
all.decomps.nemodel <- NULL
all.decomps.nemodel[[1]] <- list(c(0,0),c(1,1),c(1,0))
all.decomps.nemodel[[2]] <- list(c(0,0),c(1,0),c(1,1))
all.decomps.nemodel[[3]] <- list(c(1,0),c(0,0),c(1,1))
all.decomps.nemodel[[4]] <- list(c(1,1),c(0,0),c(0,1))
all.decomps.nemodel[[5]] <- list(c(1,1),c(0,1),c(0,0))
all.decomps.nemodel[[6]] <- list(c(0,1),c(1,1),c(0,0))

# for each decomposition: split IE1 under NE model into IE1.only and IE12
# and augment additional hypothetical exposure level a3==a1
all.decomps.4paths <- lapply(all.decomps.nemodel, function(one.decomp) {
  one.decomp.ie12 <- NULL
  for (a3 in 0:1) {
    one.decomp.a3 <- one.decomp
    # DE
    one.decomp.a3[[1]] <- c(one.decomp[[1]],one.decomp[[1]][1])
    # IE1 only
    one.decomp.a3[[2]] <- c(one.decomp[[2]],a3)
    # IE2
    one.decomp.a3[[3]] <- c(one.decomp[[3]],one.decomp[[3]][2])
    # IE12
    one.decomp.a3[[4]] <- c(one.decomp[[2]][1],1-a3,one.decomp[[2]][2])
    one.decomp.ie12[[a3+1]] <- one.decomp.a3
  }
  return(one.decomp.ie12)
})
rm(all.decomps.nemodel)

res.decomp <- NULL
rr.est <- matrix(NA,nrow=4,ncol=6)
for (ai in 1:length(all.decomps.4paths)) {
  res.decomp[[ai]] <- list()
  for (aj in 1:2) {
    one.decomp <- all.decomps.4paths[[ai]][[aj]]
    res.decomp[[ai]][[aj]] <- do.call(rbind,lapply(rho.grid,function(rho) {
      OneDecomposition(de.a123=one.decomp[[1]],
                       ie1.a023=one.decomp[[2]],
                       ie2.a013=one.decomp[[3]],
                       ie12.a012=one.decomp[[4]],
                       rho.fixed=rho,
                       res=res)
    }))
    rm(one.decomp)
  }
  rr.est[2:4,ai] <- unique(do.call(rbind,lapply(res.decomp[[ai]],function(x) 
    x[,c("de","ie1.c","ie2")])))
  rr.est[1,ai] <- prod(rr.est[2:4,ai])
}
row.names(rr.est) <- c("total","direct","SGA","Preterm")
round(rr.est,2)

pdf("MM_ABPL_PND-4sensitivity.pdf",height=18,width=8)
par(mfrow=c(6,2))
my_ylim <- range(lapply(1:length(res.decomp), function(ai) 
  lapply(res.decomp[[ai]],function(x) x[,c("ie1","ie12","ie1.c")])))
for (ai in 1:length(res.decomp)) {
  for (aj in 1:2) {
    plot(rho.grid,res.decomp[[ai]][[aj]][,"ie1.c"],
         type="n",ylim=my_ylim,
         xlab=expression(paste("Cross-world corr. ",rho)),ylab="RR",
         main=paste("Decomposition",ai))
    points(rho.grid,res.decomp[[ai]][[aj]][,"ie1"],type="b",pch=1)
    points(rho.grid,res.decomp[[ai]][[aj]][,"ie12"],type="b",pch=2)
    points(rho.grid,rep(rr.est["SGA",ai],length(rho.grid)),type="b",pch=20)
    abline(h=1,lty=3)
    if (ai==1 || ai==6) {
      legend("topright",bty="n",
             legend=c("only M1","only M1 -> M2","combined via M1"),
             pch=c(1:2,20))
    }
  }
}
dev.off()
    
