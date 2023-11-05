rm(list=ls())
library("data.table")
source("helper-funs-MC.R")

# subfolder <- "MM_CSL-bootstrap/"
subfolder <- "MM_CSL-bootstrap-PTD34/"
# subfolder <- "MM_CSL-bootstrap-completecases/"
# subfolder <- "MM_CSL-bootstrap-PTD34-completecases/"
myfiles <- list.files(subfolder)
myfiles <- myfiles[grep(pattern=".Rdata",myfiles)]
Yname <- "PND"
# Yname <- "SB"
# Yname <- "NND"
myfiles <- myfiles[grep(pattern=Yname,myfiles)]
boot.out <- NULL
for (ll in myfiles) {
  load(paste0(subfolder,ll))
  boot.out <- c(boot.out, list(
    cbind("decomp"=1:6,do.call(rbind,PROP_MEDIATED_logRR(EST)))))
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

lapply(RES.TABLES,round,1)

library("writexl")
out.filename <- paste0("MM_CSL-ABPL_",Yname,"-bootstraps.xlsx")
if (grepl("PTD34",subfolder)) {
  out.filename <- paste0("MM_CSL-PTD34-ABPL_",Yname,"-bootstraps.xlsx")
}
if (grepl("completecases",subfolder)) {
  out.filename <- paste0("MM_CSL-completecases-ABPL_",Yname,"-bootstraps.xlsx")
  if (grepl("PTD34",subfolder)) {
    out.filename <- paste0("MM_CSL-completecases-PTD34-ABPL_",Yname,"-bootstraps.xlsx")
  }
}
write_xlsx(round(do.call(cbind,RES.TABLES),2),out.filename)
           