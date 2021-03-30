rm(list=ls())
library("data.table")
library("nnet")

# for parallel jobs ===========================================================
args <- 1
if (!grepl("apple",sessionInfo()[[1]]$platform)) {
  args <- commandArgs(trailingOnly=TRUE) # for CMD BATCH '--args 1'
}
(seed <- as.integer(args[1]))

nseeds <- 2500
(seed <- (seed %% nseeds))

load("MM_ABPL_PND-1prep-M2ordinal.Rdata")

OneNEmodelEstimator <- function(data) {
  # logistic regression model for binary M1
  fitM1 <- glm(as.formula(paste0("M1~",paste(paste0("A*",l_names),
                                             collapse="+"))),
               family = binomial("logit"), data = data, weights= WGT)
  # multinomial logistic regression model for ordinal M2, given M1
  fitM2 <- nnet::multinom(as.formula(paste0("M2~",paste(paste0("A*M1*",l_names),
                                                        collapse="+"))),
                          data = data, weights= WGT, 
                          maxit = 100000, MaxNWts = 100000)
  # logistic regression model for binary Y, given A, M1, M2
  fitY <- glm(as.formula(paste0("Y~",paste(paste0("A*M1*M2*",l_names),
                                           collapse="+"))),
              family = binomial("logit"), data = data, weights= WGT)
  # model for the probability of exposure
  fitA <- glm(as.formula(paste0("A~",paste(l_names,collapse="+"))),
              family = binomial("logit"), data = data, weights= WGT)
  
  # expanded dataset: fitting model for M1
  n <- nrow(data)
  a_aux <- lapply(1:n, function(i) {
    A <- data$A[i]
    cbind(data[rep(i,4),],
          "replicate"=1:4,
          "a0"=rep(c(A,1-A),times=2),
          "a1"=rep(c(A,1-A),each=2),
          "a2"=rep(A,4)) # observed exposure A
  })
  dat1 <- data.frame(do.call(rbind,a_aux))
  rm(a_aux)
  
  # expanded dataset: fitting model for M2
  dat2 <- dat1
  colnames(dat2)[colnames(dat2)=="a2"] <- "a1.new"
  colnames(dat2)[colnames(dat2)=="a1"] <- "a2"
  colnames(dat2)[colnames(dat2)=="a1.new"] <- "a1"
  
  # estimated weights using model for M1
  pM1a1.hat <- predict.glm(fitM1,
                           newdata=data.frame("A"=dat1$a1,dat1[,l_names]),
                           type="response")
  W1.numer <- dbinom(dat1$M1, size=1, prob=pM1a1.hat)
  rm(pM1a1.hat)
  
  pM1a2.hat <- predict.glm(fitM1,
                           newdata=data.frame("A"=dat1$a2,dat1[,l_names]),
                           type="response")
  W1.denom <- dbinom(dat1$M1, size=1, prob=pM1a2.hat)
  rm(pM1a2.hat)
  
  W1 <- W1.numer/W1.denom
  
  # estimated weights using model for M2
  pM2a2.hat <- predict(fitM2,
                       newdata=data.frame("A"=dat2$a2,"M1"=dat2$M1,
                                          dat2[,l_names]),
                       type="probs")
  W2.numer <- sapply(1:length(dat2$M2), function(i) pM2a2.hat[i,dat2$M2[i]])
  rm(pM2a2.hat)
  
  pM2a1.hat <- predict(fitM2,
                       newdata=data.frame("A"=dat2$a1,"M1"=dat2$M1,
                                          dat2[,l_names]),
                       type="probs")
  W2.denom <- sapply(1:length(dat2$M2), function(i) pM2a1.hat[i,dat2$M2[i]])
  rm(pM2a1.hat)
  
  W2 <- W2.numer/W2.denom
  
  # predicted outcomes
  dat1[, "Yhat"] <- predict.glm(fitY,newdata=data.frame(
    "A"=dat1$a0,dat1[,l_names],"M1"=dat1$M1,"M2"=dat1$M2),
    type="response")
  dat2[, "Yhat"] <- predict.glm(fitY,newdata=data.frame(
    "A"=dat2$a0,dat2[,l_names],"M1"=dat2$M1,"M2"=dat2$M2),
    type="response")
  
  ## inverse probability of treatment weights for population-average effects
  pA_hat <- predict.glm(fitA,newdata=dat1,type="response")
  W1 <- W1/dbinom(dat1$A, size=1, prob=pA_hat)
  rm(pA_hat)
  pA_hat <- predict.glm(fitA,newdata=dat2,type="response")
  W2 <- W2/dbinom(dat2$A, size=1, prob=pA_hat)
  rm(pA_hat)
  ## multiply by aggregated weights
  W1 <- W1*dat1$WGT
  W2 <- W2*dat2$WGT
  
  # expanded data with weights
  dat1 <- cbind(dat1,W1)
  dat2 <- cbind(dat2,W2)
  rm(W1,W2)
  
  # fit NE model
  fit_mod <- list()
  fit_mod[[1]] <- glm(Yhat ~ a0 * a1 * a2,
                      family = binomial("logit"), data = dat1, weights = W1)
  fit_mod[[2]] <- glm(Yhat ~ a0 * a1 * a2,
                      family = binomial("logit"), data = dat2, weights = W2)
  names(fit_mod) <- paste0("fit",c("M1","M2"))
  
  return( unlist(lapply(fit_mod,coef)) )
}

set.seed(9000*seed)
nboots <- 8
res_boots <- NULL
for (bb in 1:nboots) {
  if (seed==1 && bb==1) {
    Data_boot <- Data_agg
  } else {
    idx <- sample(1:nrow(Data_agg),size=N,replace=TRUE,prob=Data_agg$WGT)
    idx_boot <- table(idx); rm(idx)
    Data_boot <- Data_agg[as.integer(names(idx_boot)),]
    Data_boot$WGT <- as.integer(idx_boot)
  }
  ptm=proc.time()[3]
  res <- OneNEmodelEstimator(data=Data_boot)
  print(res)
  cat(proc.time()[3]-ptm,"\n")
  # 274.768
  res_boots[[bb]] <- res
  rm(res)
}
if (seed==1) {
  names(res_boots)[1] <- "obs"
}
res <- res_boots; rm(res_boots)

subfolder <- "MM_ABPL_PND-2analysis-M2ordinal/"
myfile <- "NEmodels-boot-"
save(res,file=paste0(subfolder,myfile,seed,".Rdata"))
q()

# results =====================================================================
rm(list=ls())
subfolder <- "MM_ABPL_PND-2analysis-M2ordinal/"
myfiles <- list.files(subfolder)
myfiles <- myfiles[grep(pattern=".Rdata",myfiles)]
boot.out <- NULL
res_obs_all <- NULL
for (ll in myfiles) {
  load(paste0(subfolder,ll))
  if ("obs" %in% names(res)) {
    res_obs_all <- c(res_obs_all,list(res$obs)) # observed estimate
    res <- res[-1]
    names(res) <- NULL
  } 
  boot.out <- c(boot.out, res)
  rm(res)
}
res_boot_all <- do.call(rbind,boot.out); rm(boot.out)
res_obs_all <- do.call(rbind,res_obs_all)

res_boot <- res_boot_all
res_obs <- res_obs_all
# unique bootstrap samples
nrow(unique(res_boot))
# check for NAs
which(rowSums(is.na(res_boot))>0)

# calculate marginal risk ratios for each (bootstrap) estimate ================
expit <- function(x) exp(x)/(1+exp(x))
meths <- paste0("fitM",1:2)
res_all <- list()
for (meth in meths) {
  thetas <- res_obs[,grep(meth,colnames(res_obs),value=TRUE)]
  thetas_boot <- res_boot[,grep(meth,colnames(res_boot),value=TRUE)][,names(thetas)]
  thetas <- rbind("obs"=thetas,thetas_boot)
  rm(thetas_boot)
  colnames(thetas) <- unlist(lapply(
    strsplit(colnames(thetas),split=paste0(meth,".")),"[",2))
  thetas_boot <- t(apply(thetas, 1, function(theta) {
    OneDecomposition <- function(de.a1a2,ie1.a0a2,ie2.a0a1,theta) {
      # hypothetical exposure levels
      a_mtx <- rbind(c(0,de.a1a2),c(1,de.a1a2),
                     c(ie1.a0a2[1],0,ie1.a0a2[2]),c(ie1.a0a2[1],1,ie1.a0a2[2]),
                     c(ie2.a0a1,0),c(ie2.a0a1,1))
      a_mtx <- data.frame(a_mtx)
      colnames(a_mtx) <- paste0("a",0:2)
      a_mtx <- cbind(1,a_mtx,
                     a_mtx$a0*a_mtx$a1,
                     a_mtx$a0*a_mtx$a2,
                     a_mtx$a1*a_mtx$a2,
                     a_mtx$a0*a_mtx$a1*a_mtx$a2)
      colnames(a_mtx) <- names(theta)
      # marginal potential outcomes
      yhats <- expit((as.matrix(a_mtx) %*% theta)[,1])
      # risk ratios
      theta1 <- c((yhats[2]*yhats[4]*yhats[6])/(yhats[1]*yhats[3]*yhats[5]),
                  yhats[2]/yhats[1],
                  yhats[4]/yhats[3],
                  yhats[6]/yhats[5])
      names(theta1) <- paste0(c("te","de","ie1","ie2"),".rr")
      return( theta1 )
    }
    theta_all <- list(
      OneDecomposition(c(0,0),c(1,1),c(1,0),theta),
      OneDecomposition(c(0,0),c(1,0),c(1,1),theta),
      OneDecomposition(c(1,0),c(0,0),c(1,1),theta),
      OneDecomposition(c(1,1),c(0,0),c(0,1),theta),
      OneDecomposition(c(1,1),c(0,1),c(0,0),theta),
      OneDecomposition(c(0,1),c(1,1),c(0,0),theta)
    )
    names(theta_all) <- paste0("d",1:6)
    # proportion mediated depends on the exact decomposition
    for (dd in 1:6) {
      ide.dd <- theta_all[[dd]]
      if (dd==1) {
        pm.dd <- c((ide.dd["ie1.rr"]-1)*ide.dd["ie2.rr"]*ide.dd["de.rr"],
                   (ide.dd["ie2.rr"]-1)*ide.dd["de.rr"],
                   (ide.dd["de.rr"]-1))
      } else if (dd==2) {
        pm.dd <- c((ide.dd["ie1.rr"]-1)*ide.dd["de.rr"],
                   (ide.dd["ie2.rr"]-1)*ide.dd["ie1.rr"]*ide.dd["de.rr"],
                   (ide.dd["de.rr"]-1))
      } else if (dd==3) {
        pm.dd <- c((ide.dd["ie1.rr"]-1),
                   (ide.dd["ie2.rr"]-1)*ide.dd["ie1.rr"]*ide.dd["de.rr"],
                   (ide.dd["de.rr"]-1)*ide.dd["ie1.rr"])
      } else if (dd==4) {
        pm.dd <- c((ide.dd["ie1.rr"]-1),
                   (ide.dd["ie2.rr"]-1)*ide.dd["ie1.rr"],
                   (ide.dd["de.rr"]-1)*ide.dd["ie1.rr"]*ide.dd["ie2.rr"])
      } else if (dd==5) {
        pm.dd <- c((ide.dd["ie1.rr"]-1)*ide.dd["ie2.rr"],
                   (ide.dd["ie2.rr"]-1),
                   (ide.dd["de.rr"]-1)*ide.dd["ie1.rr"]*ide.dd["ie2.rr"])
      } else if (dd==6) {
        pm.dd <- c((ide.dd["ie1.rr"]-1)*ide.dd["ie2.rr"]*ide.dd["de.rr"],
                   (ide.dd["ie2.rr"]-1),
                   (ide.dd["de.rr"]-1)*ide.dd["ie2.rr"])
      }
      pm.dd <- pm.dd/(ide.dd["te.rr"]-1)
      names(pm.dd) <- paste0("propmed.",c("M1","M2","DE"))
      theta_all[[dd]] <- c(ide.dd,pm.dd*100)
    }
    unlist(theta_all)
  }))
  res.IE <- t(rbind("obs"=thetas_boot[1,], 
                    apply(thetas_boot[-1,], 2, function(x) 
                      c(#"se"=sd(x),
                        "CI"=quantile(x,probs=c(.025,.975),na.rm=TRUE)))))
  res_all <- c(res_all,list(res.IE))
}
names(res_all) <- meths

write.csv(res_all,file="MM_ABPL_PND-2analysis-M2ordinal.csv")
