BIAS_FACTOR_RR_Monly <- function(res, YUeff, pU.A0, pU.A1) {
  # res = estimated average potential outcomes
  # YUeff = gamma
  # pU.A0 = pi_0; pU.A1 = pi_1
  
  res <- data.table(res)
  setkey(res)
  # calculate bias factors
  for (i in 1:nrow(res)) {
    if (res[i,a0]==res[i,a1]) {
      BF <- 1L
    } else {
      if (res[i,a0]==0) {
        BF <- (1+(YUeff-1)*pU.A0)/(1+(YUeff-1)*pU.A1)
      }
      if (res[i,a0]==1) {
        BF <- (1+(YUeff-1)*pU.A1)/(1+(YUeff-1)*pU.A0)
      }
    }
    res[i, Y.a := Y.a/BF]
  }
  return(res)
}

BIAS_FACTOR_RR_UEIC <- function(res, YUeff, LUcor=NULL, pU.A0, pU.A1) {
  # res = estimated average potential outcomes
  # YUeff = gamma
  # LUcor = zeta
  # pU.A0 = pi_0; pU.A1 = pi_1
  
  LUcor <- ifelse(is.null(LUcor), yes=YUeff, no=LUcor)
  res <- data.table(res)
  setkey(res)
  # calculate bias factors
  for (i in 1:nrow(res)) {
    if (res[i,a0]==0) {
      BF <- (1+(YUeff*LUcor-1)*pU.A0)/(1+(LUcor-1)*pU.A0)
    } else {
      BF <- (1+(YUeff*LUcor-1)*pU.A1)/(1+(LUcor-1)*pU.A1)
    }
    if (res[i,a1]==0) {
      BF <- BF/(1+(YUeff-1)*pU.A0)
    } else {
      BF <- BF/(1+(YUeff-1)*pU.A1)
    }
    res[i, Y.a := Y.a/BF]
  }
  return(res)
}
