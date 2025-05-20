rm(list=ls())
#define packages to install and load
packages <- c("readxl","data.table","nnet")
#install all packages that are not already installed
install.packages(setdiff(packages, rownames(installed.packages())))
# load packages
lapply(packages, library, character.only = TRUE)

# Read the dataset
DATA_FILENAME <- "MM_CSL 1.xlsx"
DATA <- data.table(readxl::read_xlsx(DATA_FILENAME))
setkey(DATA)

# aggregate thin strata
DATA[INS > 2, INS := 3]
setkey(DATA)

# define PTD_TYPE as categorical variable
DATA[, PTD_TYPE := factor(PTD_TYPE)]

# keep only complete cases ####################################################
xcategorical <- c("Sitenum","AGEG","PARA_SB","SINGLE","RACE","INS","SMK",
                  "ALCH","BMIG","HTG","DIABET","GDM","YEARG")

summary(DATA[, ..xcategorical])

DATA.withNAs <- data.table(DATA) # save a copy with missing data
DATA[, DROP := FALSE]
for (xx in xcategorical) {
  DATA[is.na(get(xx)), DROP := TRUE]
  DATA[,xx] <- factor(unlist(DATA[,..xx]))
}
DATA[, table(DROP)]
DATA <- DATA[DROP==FALSE]
DATA[, DROP := NULL]
setkey(DATA)
summary(DATA[, ..xcategorical])

# risk factors in relation to preeclampsia ####################################
TABLE1 <- NULL
for (xx in c(xcategorical,"CHORIO","ABPL","SGA","PTD","INDUCT","PTD_TYPE")) {
  # raw counts by PE
  xx.tab <- table(unlist(DATA[, ..xx]),"PE"=DATA[, PE])
  # include totals
  xx.tab <- cbind("total"=rowSums(xx.tab[,1:2]),xx.tab)
  # proportions by column
  xx.tab <- data.frame(cbind(xx.tab, apply(xx.tab,2,function(x) x/sum(x)*100)))
  # reorder columns and relabel
  xx.tab <- xx.tab[,c(1,4,3,6,2,5)]
  colnames(xx.tab) <- c("tot_n","tot_pc","pe1_n","pe1_pc","pe0_n","pe0_pc")
  xx.tab <- cbind("Variable"=xx,"levels"=rownames(xx.tab),xx.tab)
  TABLE1 <- rbind(TABLE1,xx.tab)
  rm(xx.tab)
}
# library("writexl")
# write_xlsx(TABLE1,"MM_CSL-Table1-complete_cases.xlsx")

# check missingness ###########################################################
DATA.MI <- data.table(DATA.withNAs) # create missing indicator
for (xx in xcategorical) {
  DATA.MI[is.na(get(xx)), eval(xx) := 99]
  DATA.MI[,xx] <- factor(unlist(DATA.MI[,..xx]))
}
setkey(DATA.MI)
summary(DATA.MI[, ..xcategorical])

# risk factors in relation to preeclampsia ####################################
TABLE1 <- NULL
for (xx in c(xcategorical,"CHORIO","ABPL","SGA","PTD","INDUCT","PTD_TYPE")) {
  # raw counts by PE
  xx.tab <- table(unlist(DATA.MI[, ..xx]),"PE"=DATA.MI[, PE])
  # include totals
  xx.tab <- cbind("total"=rowSums(xx.tab[,1:2]),xx.tab)
  # proportions by column
  xx.tab <- data.frame(cbind(xx.tab, apply(xx.tab,2,function(x) x/sum(x)*100)))
  # reorder columns and relabel
  xx.tab <- xx.tab[,c(1,4,3,6,2,5)]
  colnames(xx.tab) <- c("tot_n","tot_pc","pe1_n","pe1_pc","pe0_n","pe0_pc")
  xx.tab <- cbind("Variable"=xx,"levels"=rownames(xx.tab),xx.tab)
  TABLE1 <- rbind(TABLE1,xx.tab)
  rm(xx.tab)
}
# library("writexl")
# write_xlsx(TABLE1,"MM_CSL-Table1-missing_indicator.xlsx")
# library("xtable")
# print(xtable(TABLE1,digits=c(0,0,0,0,1,0,1,0,1)),include.rownames=FALSE)

###############################################################################
# # Gestational age-specific risk of perinatal mortality
# # stratified by preeclampsia and normotensive persons
# # with risk ratio of mortality between preeclampsia and normotensive persons
# DATA.MI[,GA := round(GA)]
# plot.dt <- DATA.MI[, list(
#   "PE1"=.SD[PE==1,sum(PND)]/.SD[,sum(PE==1)]*1000L,
#   "PE0"=.SD[PE==0,sum(PND)]/.SD[,sum(PE==0)]*1000L,
#   "N"=.SD[,sum(PE==1)]+.SD[,sum(PE==0)]),
#   by=GA]
# plot.dt[,sum(N)]==nrow(DATA.withNAs) # check all observations included
# plot.dt[PE1==0, PE1 := NA]
# plot.dt[, RR := PE1/PE0]
# plot.dt <- plot.dt[GA<=42]
# pdf("plot-GAstratified-PE-PND.pdf",width=8,height=6)
# plot(plot.dt[,GA],plot.dt[,PE1],
#      ylim=range(plot.dt[,list(PE1,PE0,RR)],na.rm=TRUE),
#      pch=1,
#      xaxt="n",yaxt="n",log="y",
#      xlab="Gestational Age (weeks)", ylab="Perinatal Mortality per 1000 births")
# points(plot.dt[,GA],plot.dt[,PE0],pch=18)
# axis(side=1,at=plot.dt[,unique(GA)])
# axis(side=2,at=10^((-1):3),las=2,
#      labels=c("0.1","1","10","100","1000"))
# # add minor ticks to y-axis
# axis(side=2,at=10^(((-4):12)/4),las=2,lwd.ticks=0.5, labels=FALSE)
# abline(h=1,lty=3)
# lines(plot.dt[,GA],plot.dt[,PE1])
# lines(plot.dt[,GA],plot.dt[,PE0])
# points(plot.dt[,GA],plot.dt[,RR],pch=20,col=2)
# lines(plot.dt[,GA],plot.dt[,RR],col=2)
# legend("topright",legend=c("Normotensive","Preeclampsia","Risk ratio"),
#        pch=c(18,1,20),lty=1,col=c(1,1,2))
# dev.off()
# 
# rm(DATA.MI,TABLE1)