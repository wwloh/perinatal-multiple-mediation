rm(list=ls())
#define packages to install and load
packages <- c("readxl","data.table","nnet")
#install all packages that are not already installed
install.packages(setdiff(packages, rownames(installed.packages())))
# load packages
lapply(packages, library, character.only = TRUE)

# Read the dataset
DATA_FILENAME <- "MM_CSL 1.xlsx"
# DATA_FILENAME <- "MM_CSL_PTD34.xlsx"
DATA <- data.table(readxl::read_xlsx(DATA_FILENAME))
setkey(DATA)

# aggregate thin strata
DATA[INS > 2, INS := 3]
setkey(DATA)

# define mediators
DATA[, PTD_TYPE := as.factor(PTD_TYPE)]
DATA[, table(PTD_TYPE)]
DATA[, PTD_SPT := (PTD_TYPE==1)*1L]
DATA[, PTD_IND := (PTD_TYPE==2)*1L]
setkey(DATA)

# missing data as separate level for categorical predictors ###################
xcategorical <- c("Sitenum","INS","SMK","DIABET","GDM",
                  "RACE","CHORIO","PARA_SB","YEARG","HTG",
                  "AGEG","SINGLE","BMIG","PE")
summary(DATA[, ..xcategorical])
# certain predictors as factors
for (xx in xcategorical) {
  DATA[is.na(get(xx)), eval(xx) := 99]
  DATA[,xx] <- factor(unlist(DATA[,..xx]))
}
summary(DATA[, ..xcategorical])

# exposure, mediators and outcome #############################################

# define neonatal death
DATA[SB==0, NND.bin := PND]
table(DATA[, list(NND,NND.bin)], useNA="ifany")
DATA[, NND.bin := NULL]
table(DATA[, list(NND,PND,SB)], useNA="ifany")

table(DATA[, ABPL])
table(DATA[, PTD_TYPE])
table(DATA[, PND])
table(DATA[, SB])
table(DATA[, NND], useNA="ifany")

if (grepl("PTD34",DATA_FILENAME)) {
  table(DATA[, PTD_TYPE34])
}
