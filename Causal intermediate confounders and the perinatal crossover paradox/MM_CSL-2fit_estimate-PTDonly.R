source("MM_CSL-1prep-completecases.R")
source("helper-funs-MC.R")

# define variables for analysis ###############################################
Cnames <- xcategorical # covariates

# exposure
Aname <- "PE" # all PE

# mediator
# Mnames <- c("PTD") # single mediator; PTD: 1 = either SPT or IND
Mnames <- c("PTD_TYPE") # PTD_TYPE: 1 = SPT, 2 = IND

# outcome
Yname <- "PND"

# complete case analysis or single imputation #################################
COMPLETE.CASE <- FALSE
if (COMPLETE.CASE==FALSE) {
  rm(DATA)
  DATA <- DATA.withNAs
  rm(DATA.withNAs)
}

# bootstrap ###################################################################
# initialize for parallel jobs on cluster
args <- 1
if (!grepl("apple",sessionInfo()[[1]]$platform)) {
  args <- commandArgs(trailingOnly=TRUE) # for CMD BATCH '--args 1'  
}
(seed <- as.integer(args[1]))

set.seed(07102+seed)
if (seed>1) {
  # resample with replacement for bootstrap
  N <- nrow(DATA)
  boot.i <- sort(sample(x=N,size=N,replace=TRUE))
  DATA <- DATA[boot.i]
  setkey(DATA)
  rm(N,boot.i)
}


# dependent variable for each model ###########################################
var.list <- NULL
var.list[[Yname]] <- Yname
for (mm in 1:length(Mnames)) {
  var.list[[Mnames[mm]]] <- Mnames[mm]
}
var.list

# for saving results ##########################################################
if (COMPLETE.CASE==TRUE) {
  myfilename <- "MMEIC_CSL-bootstrap/MMEIC-bootstrap"  
} else {
  myfilename <- "MMEIC_CSL-imputed-bootstrap/MMEIC-bootstrap"
  # single imputation of missing data
  source("MM_CSL-1prep-missingdata.R")
}
myfilename <- paste(c(
  myfilename,
  "PTDonly",
  Aname,
  grep("PTD",Mnames,value=TRUE),
  Yname,
  paste0("seed_",seed,".Rdata")
),collapse="-")
myfilename

# fit regression models to observed data #####################################
# formula for (multinomial) logistic regression models
fit_MY <- lapply(var.list, function(vars) {
  DV <- vars[1]
  if (DV==Yname) {
    # outcome model given exposure, mediators, and covariates
    form <- as.formula(paste0(Yname,"~",paste(c(
      combn(c(Aname,Mnames),m=2,FUN=paste,collapse="*"),Cnames),
      collapse = "+")))
    glm(formula = form, data = DATA, family = binomial("logit"))
  } else {
    # mediator model given exposure and covariates
    form <- as.formula(paste0(DV,"~",paste(
      c(Aname,Cnames),collapse = "+")))
    if (DV=="PTD_TYPE") {
      # multinomial logistic regression
      nnet::multinom(formula=form, data=DATA, maxit=1e4, trace=FALSE)  
    } else {
      # logistic regression
      glm(formula = form, data = DATA, family = binomial("logit"))  
    }
  }
})
names(fit_MY)
# check coefficient estimates of fitted models
lapply(fit_MY, function(fit.logit) {
  if (class(fit.logit)[1]=="glm") {
    round(summary(fit.logit)$coef[,"Estimate",drop=FALSE],2)  
  } else if (class(fit.logit)[1]=="multinom") {
    round(t(coef(fit.logit)),2)
  }
})

if ("PTD_TYPE" %in% Mnames) {
  # check convergence of multinomial logistic model
  fit_MY[[grep("PTD",Mnames,value=TRUE)]]$convergence
  if (fit_MY[[grep("PTD",Mnames,value=TRUE)]]$convergence) 
    q() # exit if failed to converge
}

# prepare data and fitted models for estimating (in)direct effects ############
Data=data.frame(DATA)
Data=cbind(Data,"id"=1:nrow(Data))

ptm=proc.time()[3]
EST <- OneMCestimator(Data,mc_draws=200,fit_MY,Aname,Mnames,Yname,Cnames)
round((proc.time()[3]-ptm)/60)
# 311

save(EST,file=myfilename)
EST

q()
