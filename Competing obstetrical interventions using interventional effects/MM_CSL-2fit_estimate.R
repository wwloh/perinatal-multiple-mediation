source("MM_CSL-1prep.R")
source("helper-funs-MC.R")

Lnames <- xcategorical # covariates
Aname <- "ABPL" # exposure
Mnames <- c("PTD_IND","PTD_SPT") # mediators
# Yname <- "NND" # outcome
# Yname <- "PND"
Yname <- "SB"

if (grepl("PTD34",DATA_FILENAME)) {
  Mnames <- c("PTD_IND34","PTD_SPT34") # mediators
}

if (Yname=="NND") {
  # remove missing outcomes
  DATA <- DATA[!is.na(NND)]
}

# bootstrap ###################################################################
# initialize for parallel jobs on cluster
args <- 1
if (!grepl("apple",sessionInfo()[[1]]$platform)) {
  args <- commandArgs(trailingOnly=TRUE) # for CMD BATCH '--args 1'  
}
(seed <- as.integer(args[1]))

set.seed(30322+seed)
if (seed>1) {
  # resample with replacement for bootstrap
  N <- nrow(DATA)
  boot.i <- sort(sample(x=N,size=N,replace=TRUE))
  DATA <- DATA[boot.i]
  setkey(DATA)
  rm(N,boot.i)
}

myfilename <- "MM_CSL-bootstrap/MM_CSL-bootstrap-"
if (grepl("PTD34",DATA_FILENAME)) {
  myfilename <- "MM_CSL-bootstrap/MM_CSL-PTD34-bootstrap-"
}
(myfilename <- paste0(myfilename,Aname,"-",Yname,"-seed_",seed,".Rdata"))

# fit regression models to observed data #####################################
# variables for each model 
var.list <- NULL
var.list[[Yname]] <- c(Yname,Aname,Mnames,Lnames)
var.list[[Mnames[1]]] <- c(Mnames[1],Aname,Lnames)
var.list[[Mnames[2]]] <- c(Mnames[2],Aname,Lnames)

# formula for logistic regression models
fit_MY <- lapply(var.list, function(vars) {
  DV <- vars[1]
  if (DV==Yname) {
    # outcome model given exposure, mediators, and covariates
    form <- as.formula(paste0(Yname,"~",paste(
      c(paste0(Aname,"*",Mnames),Lnames),collapse = "+")))
  } else {
    # mediator model given exposure and covariates
    form <- as.formula(paste0(DV,"~",paste(
      c(Aname,Lnames),collapse = "+")))
  }
  glm(formula = form, data = DATA, family = binomial("logit"))
})
names(fit_MY)
# check coefficient estimates of fitted models
lapply(fit_MY, function(fit.glm) round(summary(fit.glm)$coef,2))

# point-estimates and 95% CIs for E-values
FITTED.MODELS <- lapply(fit_MY, function(fit.glm) 
  round(cbind(
    "Est"=summary(fit.glm)$coef[,"Estimate"],
    "95% CI (lower)"=summary(fit.glm)$coef[,"Estimate"]-
      qnorm(.975)*summary(fit.glm)$coef[,"Std. Error"],
    "95% CI (upper)"=summary(fit.glm)$coef[,"Estimate"]+
      qnorm(.975)*summary(fit.glm)$coef[,"Std. Error"]
  ),2)
)
write.csv(FITTED.MODELS[[Yname]], 
          file=paste0("mm mediation models - fitted_",Yname,".csv"))
write.csv(FITTED.MODELS[[Mnames[1]]], 
          file=paste0("mm mediation models - fitted_",Mnames[1],".csv"))
write.csv(FITTED.MODELS[[Mnames[2]]], 
          file=paste0("mm mediation models - fitted_",Mnames[2],".csv"))


# prepare data and fitted models for estimating (in)direct effects ############
Data=data.frame(DATA)
Data=cbind(Data,"id"=1:nrow(Data))

ptm=proc.time()[3]
EST <- OneMCestimator(Data,mc_draws=200,fit_MY,Aname,Mnames,Yname,Lnames)
round((proc.time()[3]-ptm)/60)
# 160 minutes

save(EST,file=myfilename)
EST
round(do.call(rbind,PROP_MEDIATED_logRR(EST)),1)
q()
