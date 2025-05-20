rm(list=ls())
library("data.table")
library("writexl")
source("helper-funs-MC.R")

# subfolder <- "MMEIC_CSL-bootstrap/" # complete case analysis
subfolder <- "MMEIC_CSL-imputed-bootstrap/" # single stochastic imputation

Aname <- "-PE"
Mnames <- c("PTD-","PTD_TYPE")
Yname <- "PND"

meths <- c("PTDonly","withEIC","MMs")

Mname <- Mnames[2]
for (meth in meths) {

myfiles <- list.files(subfolder)
myfiles <- myfiles[grep(pattern=".Rdata",myfiles)]

myfiles <- myfiles[grep(pattern=meth,myfiles)]
myfiles <- myfiles[grep(pattern=Aname,myfiles,fixed=TRUE)]
myfiles <- myfiles[grep(pattern=Mname,myfiles,fixed=TRUE)]
myfiles <- myfiles[grep(pattern=Yname,myfiles,fixed=TRUE)]
myfiles

boot.out <- NULL
for (ll in myfiles) {
  load(paste0(subfolder,ll))
  if (nrow(EST)==4) {
    EST <- do.call(rbind,PROP_MEDIATED_logRR_singleM(EST))
  } else {
    EST <- do.call(rbind,PROP_MEDIATED_logRR_4M(EST))
  }
  boot.out <- c(boot.out, list(cbind("decomp"=1:2,EST)))
  rm(EST)
}
OBS <- data.table(boot.out[[1]]) # observed estimate
setkey(OBS)
eff.names <- names(OBS)[-1]
OBS[, "stat" := "0.obs"]

RES <- NULL
RES[[1]] <- OBS

BOOT <- data.table(do.call(rbind,boot.out[-1]))
setkey(BOOT)
RES[[2]] <- BOOT[, lapply(.SD, sd), by=decomp][, "stat" := "1.se"]
RES[[3]] <- BOOT[, lapply(.SD, quantile, probs=.025), by=decomp][, "stat" := "ci.l"]
RES[[4]] <- BOOT[, lapply(.SD, quantile, probs=.975), by=decomp][, "stat" := "ci.u"]

RES <- rbindlist(RES)
setkey(RES)

RES.TABLES <- lapply(eff.names, function(en) 
  dcast(RES, decomp~stat, value.var=en))
names(RES.TABLES) <- eff.names

print(BOOT[,.N,by=decomp]) # number of bootstrap samples

print(lapply(RES.TABLES,round,1))

out.filename <- strsplit(myfiles[1],split="[.]Rdata")[[1]]
write_xlsx(round(do.call(cbind,RES.TABLES),2),paste0(out.filename,".xlsx"))
}
