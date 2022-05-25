CrossWorldCorrelation <- function(pZ1_0=c("a1"=0.1,"a3"=0.8),rho=0) {
  P1z1 <- c(0,pZ1_0["a1"],1) # CDF of Z1(a1)
  P12z1 <- c(0,pZ1_0["a3"],1) # CDF of Z1(a3)
  C12jk <- lapply(1:500, function(r) {
    C12 <- matrix(0,nrow=2,ncol=2)
    for (j in 2:3) {
      for (k in 2:3) {
        u1 <- runif(1,min=P1z1[k-1],max=P1z1[k])
        U1 <- qnorm(u1)
        U12 <- rnorm(1)
        U12C <- rho*U1 + sqrt(1-rho^2)*U12
        u12C <- pnorm(U12C)
        C12[j-1,k-1] <- (u12C > P12z1[j-1]) && (u12C <= P12z1[j])
      }
    }
    return(C12)
  })
  Pr_Z1a3_Z1a1 <- Reduce(f="+",x=C12jk)/length(C12jk)
  row.names(Pr_Z1a3_Z1a1) <- paste0("Z1a3_eq",0:1)
  colnames(Pr_Z1a3_Z1a1) <- paste0("Z1a1_eq",0:1)
  # return Pr( Z1(a3)=1|Z1(a1)=z)
  return(unlist(as.data.frame(Pr_Z1a3_Z1a1)[2,]))
}
