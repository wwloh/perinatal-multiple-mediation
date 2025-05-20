# impute missing data #########################################################
library("mice")
## covariates with missing data
x.missing <- xcategorical[apply(is.na(DATA[, ..xcategorical]),2,any)]
## single imputation of covariates
MICE.OUT <- mice(data=data.frame(DATA[, ..xcategorical]), 
                 m=1,
                 method="rf",
                 blocks=as.list(x.missing),
                 maxit=10)
for (xx in x.missing) {
  # check number of missing observations matches
  if(nrow(MICE.OUT$imp[[xx]])==nrow(DATA[is.na(get(xx)), ..xx])) {
    # replace missing values with imputations
    DATA[is.na(get(xx)), eval(xx) := MICE.OUT$imp[[xx]]]
    DATA[,xx] <- factor(unlist(DATA[,..xx]))
    cat(nrow(MICE.OUT$imp[[xx]]), "missing values for", xx, "imputed \n")
  }
}
setkey(DATA)
summary(DATA[, ..xcategorical])
rm(x.missing,MICE.OUT)
