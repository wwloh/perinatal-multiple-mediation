source("MM_CSL-1prep-completecases.R")
source("helper-funs-MC.R")

# define variables for analysis ###############################################
Cnames <- xcategorical # covariates

# exposure
Aname <- "PE" # all PE

# intermediate variables
Mnames <- c("ABPL","SGA","CHORIO","PTD_TYPE") # PTD_TYPE: 1 = SPT, 2 = IND

# indices for exposure-induced confounders
Lidx <- 1:3

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
  "withEIC",
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
      # exposure-mediator interactions
      paste0(Aname,"*",Mnames),
      # (exposure-)mediator-mediator interactions
      paste0(Aname,"*",combn(Mnames,m=2,FUN=paste,collapse="*")),
      # covariates
      Cnames),collapse = "+")))
    glm(formula = form, data = DATA, family = binomial("logit"))
  } else if (grepl(pattern="PTD",x=DV)) {
    # mediator model given exposure-induced confounders, exposure and covariates
    form <- as.formula(paste0(DV,"~",paste(c(
      # exposure-EIC interactions
      paste0(Aname,"*",Mnames[Lidx]), 
      # (exposure-)EIC-EIC interactions
      paste0(Aname,"*",combn(Mnames[Lidx],m=2,FUN=paste,collapse="*")),
      # covariates
      Cnames),collapse = "+")))
    if (DV=="PTD_TYPE") {
      # multinomial logistic regression
      nnet::multinom(formula=form, data=DATA, maxit=1e4, trace=FALSE)  
    } else {
      # logistic regression
      glm(formula = form, data = DATA, family = binomial("logit"))  
    }
  } else {
    Li <- which(DV==Mnames[Lidx])
    if (Li>1L) {
      # exposure-induced confounders given other EIC, exposure, and covariates
      form <- as.formula(paste0(DV,"~",paste(
        c(paste0(Aname,"*",Mnames[1:(Li-1)]), # exposure-EIC interactions
          Cnames),collapse = "+")))  
    } else {
      # exposure-induced confounders given exposure and covariates
      form <- as.formula(paste0(DV,"~",paste(
        c(Aname,Cnames),collapse = "+")))
    }
    # logistic regression
    glm(formula = form, data = DATA, family = binomial("logit"))
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
EST <- OneMCestimator(Data,mc_draws=200,fit_MY,Aname,Mnames,Yname,Cnames,Lidx)
round((proc.time()[3]-ptm)/60)


save(EST,file=myfilename)
EST

q()
