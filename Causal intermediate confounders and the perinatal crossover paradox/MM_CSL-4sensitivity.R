rm(list=ls())
# load estimates using observed data ##########################################
library("data.table")
library("writexl")
library("xtable")
source("helper-funs-MC.R")
source("helper-funs-sensitivity.R")

# subfolder <- "MMEIC_CSL-bootstrap/" # complete case analysis
subfolder <- "MMEIC_CSL-imputed-bootstrap/" # single stochastic imputation

Aname <- "PE"
Mname <- "PTD_TYPE"
Yname <- "PND"

meths <- c("PTDonly","withEIC","MMs")

# specify which method to use
args <- 1 # initialize
args <- commandArgs(trailingOnly=TRUE)
meth <- meths[as.integer(args[1])]
rm(args)

# load bootstrap results 
myfiles <- list.files(subfolder)
myfiles <- myfiles[grep(pattern=".Rdata",myfiles)]

myfiles <- myfiles[grep(pattern=meth,myfiles)]
myfiles <- myfiles[grep(pattern=paste0(Aname,"-"),myfiles,fixed=TRUE)]
myfiles <- myfiles[grep(pattern=paste0(Mname,"-"),myfiles,fixed=TRUE)]
myfiles <- myfiles[grep(pattern=Yname,myfiles,fixed=TRUE)]

out.filename <- strsplit(myfiles[1],split="[.]Rdata")[[1]]
cat(out.filename,"\n")

boot.out <- NULL

for (ll in myfiles) {
  load(paste0(subfolder,ll))
  boot.out <- c(boot.out, list(EST))
  rm(EST)
}

# values of sensitivity parameters based on observed associations ############
source2 <- function(file, start, end, ...) {
  file.lines <- scan(file, what=character(), skip=start-1, nlines=end-start+1, sep='\n')
  file.lines.collapsed <- paste(file.lines, collapse='\n')
  source(textConnection(file.lines.collapsed), ...)
}

## load data
source2("MM_CSL-1prep-completecases.R",2,27)

## initialize
Cnames <- xcategorical # covariates
var.list <- NULL
var.list[[Yname]] <- Yname
var.list[[Mname]] <- Mname
Mnames <- c("ABPL","SGA","CHORIO",Mname)

## single imputation of missing data
rm(DATA)
DATA <- DATA.withNAs
rm(DATA.withNAs)
set.seed(07102+1) # same seed as for observed data
source("MM_CSL-1prep-missingdata.R")

## fit regression models used for estimation
if (meth=="PTDonly") {
  Mnames <- Mname
  source2("MM_CSL-2fit_estimate-PTDonly.R",72,109)
} else if (meth=="withEIC") {
  Lidx <- 1:3
  source2("MM_CSL-2fit_estimate-PTD_EIC.R",74,134)
} else if (meth=="MMs") {
  Lidx <- NULL
  source2("MM_CSL-2fit_estimate-PTD_MMs.R",74,115)
}
# print coefficient estimates of fitted models
lapply(fit_MY, function(fit.logit) {
  if (class(fit.logit)[1]=="glm") {
    print(xtable(fit.logit,digits=2))
  } else if (class(fit.logit)[1]=="multinom") {
    # create own coefficient table
    fit.logit.EST <- coef(fit.logit)
    fit.logit.SE <- summary(fit.logit)$standard.errors
    fit.logit.Z <- fit.logit.EST/fit.logit.SE
    fit.logit.Pv <- 2*pnorm(abs(fit.logit.Z),lower.tail=FALSE)
    fit.logit.i <- data.frame(t(rbind(
      fit.logit.EST,fit.logit.SE,fit.logit.Z,fit.logit.Pv)))
    colnames(fit.logit.i) <- rep(
      c("Estimate","Std. Error","z value","Pr(>|z|)"),each=2)
    print(xtable(fit.logit.i,digits=2))
  }
})

# coefficients of outcome model
Y.coef <- coef(fit_MY[[Yname]])
Y.coef <- Y.coef[setdiff(names(Y.coef),c("(Intercept)"))]
Y.coef
range(Y.coef)

# coefficients of mediator model
M.coef <- coef(fit_MY[[grep("PTD",Mnames,value=TRUE)]])
if ("PTD_TYPE" %in% Mnames) {
  M.coef <- M.coef[,setdiff(colnames(M.coef),c("(Intercept)"))]
} else {
  M.coef <- M.coef[setdiff(names(M.coef),c("(Intercept)"))]  
}
M.coef
range(M.coef)

# prevalance of binary covariates given exposure and mediator
Y.Dmtx <- data.table(model.matrix(fit_MY[[Yname]]))
setkey(Y.Dmtx)
if ("PTD_TYPE" %in% Mnames) {
  U.A0M0.vals <- colMeans(Y.Dmtx[get(Aname)==0 & 
                                   get(paste0(Mname,"1"))==0 & 
                                   get(paste0(Mname,"2"))==0])
  U.A0M1.vals <- colMeans(Y.Dmtx[get(Aname)==0 & 
                                   get(paste0(Mname,"1"))==1 & 
                                   get(paste0(Mname,"2"))==0])
  U.A0M2.vals <- colMeans(Y.Dmtx[get(Aname)==0 & 
                                   get(paste0(Mname,"1"))==0 & 
                                   get(paste0(Mname,"2"))==1])
  U.A1M0.vals <- colMeans(Y.Dmtx[get(Aname)==1 & 
                                   get(paste0(Mname,"1"))==0 & 
                                   get(paste0(Mname,"2"))==0])
  U.A1M1.vals <- colMeans(Y.Dmtx[get(Aname)==1 & 
                                   get(paste0(Mname,"1"))==1 & 
                                   get(paste0(Mname,"2"))==0])
  U.A1M2.vals <- colMeans(Y.Dmtx[get(Aname)==1 & 
                                   get(paste0(Mname,"1"))==0 & 
                                   get(paste0(Mname,"2"))==1])
  U.A.vals <- rbind(cbind(U.A0M0.vals,U.A1M0.vals),
                    cbind(U.A0M1.vals,U.A1M1.vals),
                    cbind(U.A0M2.vals,U.A1M2.vals))
} else {
  U.A0M0.vals <- colMeans(Y.Dmtx[get(Aname)==0 & get(Mname)==0])
  U.A0M1.vals <- colMeans(Y.Dmtx[get(Aname)==0 & get(Mname)==1])
  U.A1M0.vals <- colMeans(Y.Dmtx[get(Aname)==1 & get(Mname)==0])
  U.A1M1.vals <- colMeans(Y.Dmtx[get(Aname)==1 & get(Mname)==1])
  U.A.vals <- rbind(cbind(U.A0M0.vals,U.A1M0.vals),
                    cbind(U.A0M1.vals,U.A1M1.vals))
}
## keep baseline covariates only
Cidx <- !grepl("(Intercept)", row.names(U.A.vals)) &
  !grepl(Aname, row.names(U.A.vals)) &
  apply(sapply(Mnames, function(m) !grepl(m, row.names(U.A.vals))), 1, all)
U.A.vals[Cidx,]
summary(U.A.vals[Cidx,])

# sensitivity parameters for unmeasured confounding ###########################
SENS_PARVALS <- expand.grid(
  "YUeff"=exp(seq(from=0, to=log(6), length.out=21)),
  "pU.A0"=quantile(U.A.vals[Cidx,1],probs=seq(from=0,to=1,length.out=16)),
  "pU.A1"=quantile(U.A.vals[Cidx,2],probs=seq(from=0,to=1,length.out=16))
)
nrow(SENS_PARVALS)

# apply bias formula to average potential outcomes and effects estimates
boot.BF.out <- NULL
ptm <- proc.time()[3]
for (bfi in 1:nrow(SENS_PARVALS)) {
  boot.bfi <- lapply(boot.out, function(EST) {
    if (meth=="PTDonly") {
      EST_BFI <- BIAS_FACTOR_RR_Monly(
        res=EST,
        YUeff=SENS_PARVALS[bfi,"YUeff"],
        pU.A0=SENS_PARVALS[bfi,"pU.A0"],
        pU.A1=SENS_PARVALS[bfi,"pU.A1"])
    } else {
      if (meth=="MMs") {
        # keep indirect effect via PTD only
        if (nrow(EST)==2^(3+1)) {
          # 3 mediators
          EST <- EST[a0==a1 & a0==a2]
          EST[, c("a1","a2") := NULL]
          setnames(EST,"a3","a1")
        } else if (nrow(EST)==2^(4+1)) {
          # 4 mediators
          EST <- EST[a0==a1 & a0==a2 & a0==a3]
          EST[, c("a1","a2","a3") := NULL]
          setnames(EST,"a4","a1")
        }
      }
      EST_BFI <- BIAS_FACTOR_RR_UEIC(
        res=EST,
        YUeff=SENS_PARVALS[bfi,"YUeff"],
        pU.A0=SENS_PARVALS[bfi,"pU.A0"],
        pU.A1=SENS_PARVALS[bfi,"pU.A1"])
    }
    EST_BFI <- do.call(rbind,PROP_MEDIATED_logRR_singleM(EST_BFI))
    return(cbind("decomp"=1:2,EST_BFI))
  })
  
  # same procedure to combine bootstraps as without sensitivity analysis ######
  OBS <- data.table(boot.bfi[[1]]) # observed estimate
  setkey(OBS)
  eff.names <- names(OBS)[-1]
  OBS[, "stat" := "0.obs"]
  
  RES <- NULL
  RES[[1]] <- OBS
  
  boot.BF <- data.table(do.call(rbind,boot.bfi[-1]))
  setkey(boot.BF)
  RES[[2]] <- boot.BF[, lapply(.SD, sd), by=decomp][, "stat" := "1.se"]
  RES[[3]] <- boot.BF[, lapply(.SD, quantile, probs=.025), by=decomp][, "stat" := "ci.l"]
  RES[[4]] <- boot.BF[, lapply(.SD, quantile, probs=.975), by=decomp][, "stat" := "ci.u"]
  
  RES <- rbindlist(RES)
  setkey(RES)
  
  RES.TABLES <- lapply(eff.names, function(en) 
    dcast(RES, decomp~stat, value.var=en))
  names(RES.TABLES) <- eff.names
  #############################################################################
  rm(OBS,eff.names,RES,boot.BF)
  boot.BF.out[[bfi]] <- cbind(SENS_PARVALS[rep(bfi,2),],
                              do.call(cbind,RES.TABLES))
  rm(boot.bfi,RES.TABLES)
  cat(bfi,"out of",nrow(SENS_PARVALS),"done;",
      round((proc.time()[3]-ptm)/60), "mins \n")
}

# select values of sensitivity parameters (other than gamma)
SENS_RES.ALL <- rbindlist(boot.BF.out)
setkey(SENS_RES.ALL)

# check that certain values return the same observed estimates
SENS_RES.ALL[YUeff==1 & de.rr.decomp==1, unique(de.rr.0.obs)]
SENS_RES.ALL[YUeff==1 & de.rr.decomp==2, unique(de.rr.0.obs)]

# find proportions closest to observed associations among DE>=1
pi.cands <- SENS_RES.ALL[de.rr.0.obs>=1 & de.rr.decomp==2, 
                         unique(list(pU.A0,pU.A1))]
pi.sqmin <- apply(pi.cands,1,function(pu) {
  min(sqrt(colMeans((unlist(pu)-t(U.A.vals[Cidx,]))^2)))
})
pi.cands[which.min(pi.sqmin)]

pU.A0.fixed <- unlist(pi.cands[which.min(pi.sqmin)][,1])
pU.A1.fixed <- unlist(pi.cands[which.min(pi.sqmin)][,2])
SENS_RES <- SENS_RES.ALL[abs(pU.A0-pU.A0.fixed)<1e-6 & 
                           abs(pU.A1-pU.A1.fixed)<1e-6 & 
                           de.rr.decomp==2]
setkey(SENS_RES)
rm(boot.BF.out)

# make plots ##################################################################
# load(file=paste0("MM_CSL-4sensitivity-",out.filename,".Rdata"))
pdf(paste0("plot-sensitivity-",out.filename,"-prevalences.pdf"),
    width=5,height=5)
plot(U.A.vals[Cidx,],
     pch=4,xlim=c(0,1),ylim=c(0,1), col="gray",
     xlab="p(C|A=0,M)",ylab="p(C|A=1,M)",
     main="Prevalences of covariates C")
abline(a=0,b=1,lty=3)
points(pU.A0.fixed,pU.A1.fixed,pch=20)
dev.off()

dd <- 2 # selected decomposition
SENS_RES.p <- SENS_RES[de.rr.decomp==dd, list(
  YUeff,
  de.rr.0.obs,de.rr.ci.l,de.rr.ci.u,
  ie1.rr.0.obs,ie1.rr.ci.l,ie1.rr.ci.u)]
pdf(paste0("plot-sensitivity-",out.filename,"-effect_ests-decomp-",dd,".pdf"),
    width=5,height=4)
plot(SENS_RES.p[,YUeff],SENS_RES.p[,de.rr.0.obs],
     pch=18, cex=1.1, 
     ylim=range(SENS_RES.p[,-1]),
     log="xy",yaxt="n",
     xlab=expression(gamma), ylab="Effect Estimate", 
     main="Sensitivity Analysis")
axis(side=2,las=2,at=(c(1:4,6,8,10)*0.25))
points(SENS_RES.p[,YUeff],SENS_RES.p[,ie1.rr.0.obs],
       pch=20, cex=1.1,
       col=2)
abline(h=1,v=1, lty=3, lwd=0.5)
# 95% CIs
for (j in 1:nrow(SENS_RES.p)) {
  lines(rep(SENS_RES.p[j,YUeff],2),
        unlist(SENS_RES.p[j,list(de.rr.ci.l,de.rr.ci.u)]),
        col=1)
  lines(rep(SENS_RES.p[j,YUeff],2),
        unlist(SENS_RES.p[j,list(ie1.rr.ci.l,ie1.rr.ci.u)]),
        col=2)
}
legend("bottomright",
       legend= c("Indirect", "Direct"),
       col=c(2:1),pch=c(20,18),bty="n",cex=1.1)
dev.off()

print(SENS_RES.p[de.rr.0.obs>=1, min(YUeff,na.rm=TRUE)])

save.image(file=paste0("MM_CSL-4sensitivity-",out.filename,".Rdata"))
