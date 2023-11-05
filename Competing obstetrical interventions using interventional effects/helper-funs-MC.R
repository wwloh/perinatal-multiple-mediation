# helper function to create duplicated data for one individual
Dupdata <- function(t) {
  # all possible combinations of hypothetical treatments
  out <- expand.grid(lapply(1:(t+1L), function(a.t) 0:1))
  colnames(out) <- paste0("a",0:t)
  out <- cbind("a.i"=1:nrow(out),out)
  return(out)
}

OneMCestimator <- function(Data,mc_draws=100,fit_MY,
                           Aname,Mnames,Yname,Lnames) {
  # number of distinct mediators
  t <- length(Mnames)
  # duplicated data for each individual
  dat <- data.table(Data)
  setkey(dat)
  alevels <- dat[,as.data.table(Dupdata(t=t)),by=id]
  setkey(alevels)
  dat <- merge(dat,alevels,all.x=TRUE)
  setkey(dat)
  rm(alevels)
  
  SampleMs <- function(mydt,av_mc=TRUE) {
    # mydt is a data.table for each id containing the relevant columns
    # e.g., mydt <- dat[id==1]
    mydt_mc <- mydt[rep(1:nrow(mydt),each=mc_draws)]
    mydt_mc <- cbind("mc"=rep(1:mc_draws,times=nrow(mydt)),mydt_mc)
    setkey(mydt_mc)
    
    # newdata for imputation: same observed L for all duplicated rows
    newdataMs <- data.frame(mydt_mc[1, ..Lnames])
    
    # initialize counterfactual mediator draws
    for (s in 1:t) {
      # exposure group-specific fitted models
      for (aa in 0:1) {
        fitMs_Aa <- fit_MY[[Mnames[s]]]
        # relevant rows in duplicated data
        s_aa <- mydt_mc[, get(paste0("a",s))]==aa
        sampMs_Aa <- sum(s_aa)
        # predicted counterfactual mediator
        newdataMs.aa <- cbind(aa,newdataMs)
        colnames(newdataMs.aa)[1] <- Aname
        meanMs_Aa <- predict.glm(
          object=fitMs_Aa,newdata=newdataMs.aa,type="response")
        if (fitMs_Aa$family$family=="binomial") {
          # binomial distribution
          drawnMs_Aa <- rbinom(n=sampMs_Aa,size=1,prob=meanMs_Aa)
        }
        mydt_mc[s_aa, paste0(Mnames[s],".a") := drawnMs_Aa]
        rm(drawnMs_Aa,meanMs_Aa,newdataMs.aa,fitMs_Aa,s_aa,sampMs_Aa)
      }
    }
    setkey(mydt_mc)
    
    # helper function to predict Y for different outcome models
    PredictY <- function(onedat) {
      setnames(onedat,"a0",Aname)
      Ya <- rep(NA_real_, nrow(onedat))
      # exclude observations with both mediators being 1
      indexY_am <- !apply(onedat[,..Mnames]==1,1,all)
      if(any(indexY_am)) {
        check.am <- tryCatch(system.time(
          pred.Ya <- predict.glm(
            object=fit_MY[[Yname]],
            newdata=data.frame(onedat[indexY_am]),type="response")
        )[3], error=function(cond) return(NA))
        if (!is.na(check.am)) {
          Ya[indexY_am] <- pred.Ya  
        }
      }
      return(Ya)
    }
    
    ## predict potential outcomes using randomly sampled mediator values
    mydt_mc[, eval(Mnames) := NULL]
    setnames(mydt_mc,paste0(Mnames,".a"),Mnames)
    setkey(mydt_mc)
    mydt_mc[, "Y.a" := PredictY(
      onedat=mydt_mc[, c("a0",Lnames,Mnames), with=FALSE])]
    setkey(mydt_mc)
    # set to observed outcome for observed treatments
    mydt_mc[apply(mydt_mc[,paste0("a",0:t),with=FALSE]==
                    unique(unlist(mydt_mc[,..Aname])),1,all),
            Y.a := get(Yname)]
    
    if (av_mc==TRUE) {
      ## average over MC draws
      mydt_mc_means <- mydt_mc[,lapply(.SD, mean, na.rm=TRUE),
                               by=c("a.i",paste0("a",0:t)),.SDcols="Y.a"]
      setkey(mydt_mc_means)
      return(mydt_mc_means[,Y.a])
    } else {
      return(mydt_mc)  
    }
  }
  
  # average potential outcomes for each individual only
  dat[, "Y.a" := SampleMs(.SD,av_mc=TRUE), by=id]
  setkey(dat)
  
  # average potential outcomes for each combination of hypothetical treatments
  res <- dat[, lapply(.SD, mean, na.rm=TRUE),
             by=c("a.i",paste0("a",0:t)),.SDcols="Y.a"]
  setkey(res)
  
  return( res )
}

PROP_MEDIATED <- function(res) {
  res <- data.table(res)
  setkey(res)
  OneDecomposition <- function(de.a1a2,ie1.a0a2,ie2.a0a1,res) {
    c("de.rr"=
        res[a0==1 & a1==de.a1a2[1] & a2==de.a1a2[2], Y.a]/
        res[a0==0 & a1==de.a1a2[1] & a2==de.a1a2[2], Y.a],
      "ie1.rr"=
        res[a0==ie1.a0a2[1] & a1==1 & a2==ie1.a0a2[2], Y.a]/
        res[a0==ie1.a0a2[1] & a1==0 & a2==ie1.a0a2[2], Y.a],
      "ie2.rr"=
        res[a0==ie2.a0a1[1] & a1==ie2.a0a1[2] & a2==1, Y.a]/
        res[a0==ie2.a0a1[1] & a1==ie2.a0a1[2] & a2==0, Y.a],
      "te.rr"=
        res[a0==1 & a1==1 & a2==1, Y.a]/res[a0==0 & a1==0 & a2==0, Y.a])
  }
  DECOMP <- list(
    OneDecomposition(c(0,0),c(1,1),c(1,0),res),
    OneDecomposition(c(0,0),c(1,0),c(1,1),res),
    OneDecomposition(c(1,0),c(0,0),c(1,1),res),
    OneDecomposition(c(1,1),c(0,0),c(0,1),res),
    OneDecomposition(c(1,1),c(0,1),c(0,0),res),
    OneDecomposition(c(0,1),c(1,1),c(0,0),res)
  )
  # proportion mediated depends on the exact decomposition
  for (dd in 1:6) {
    ide.dd <- DECOMP[[dd]]
    if (dd==1) {
      pm.dd <- c((ide.dd["ie1.rr"]-1)*ide.dd["ie2.rr"]*ide.dd["de.rr"],
                 (ide.dd["ie2.rr"]-1)*ide.dd["de.rr"],
                 (ide.dd["de.rr"]-1))
    } else if (dd==2) {
      pm.dd <- c((ide.dd["ie1.rr"]-1)*ide.dd["de.rr"],
                 (ide.dd["ie2.rr"]-1)*ide.dd["ie1.rr"]*ide.dd["de.rr"],
                 (ide.dd["de.rr"]-1))
    } else if (dd==3) {
      pm.dd <- c((ide.dd["ie1.rr"]-1),
                 (ide.dd["ie2.rr"]-1)*ide.dd["ie1.rr"]*ide.dd["de.rr"],
                 (ide.dd["de.rr"]-1)*ide.dd["ie1.rr"])
    } else if (dd==4) {
      pm.dd <- c((ide.dd["ie1.rr"]-1),
                 (ide.dd["ie2.rr"]-1)*ide.dd["ie1.rr"],
                 (ide.dd["de.rr"]-1)*ide.dd["ie1.rr"]*ide.dd["ie2.rr"])
    } else if (dd==5) {
      pm.dd <- c((ide.dd["ie1.rr"]-1)*ide.dd["ie2.rr"],
                 (ide.dd["ie2.rr"]-1),
                 (ide.dd["de.rr"]-1)*ide.dd["ie1.rr"]*ide.dd["ie2.rr"])
    } else if (dd==6) {
      pm.dd <- c((ide.dd["ie1.rr"]-1)*ide.dd["ie2.rr"]*ide.dd["de.rr"],
                 (ide.dd["ie2.rr"]-1),
                 (ide.dd["de.rr"]-1)*ide.dd["ie2.rr"])
    }
    pm.dd <- pm.dd/(ide.dd["te.rr"]-1)
    names(pm.dd) <- paste0("propmed.",c("M1","M2","DE"))
    DECOMP[[dd]] <- c(ide.dd,pm.dd*100)
  }
  return(DECOMP)
}

PROP_MEDIATED_logRR <- function(res) {
  res <- data.table(res)
  setkey(res)
  OneDecomposition <- function(de.a1a2,ie1.a0a2,ie2.a0a1,res) {
    c("de.rr"=
        res[a0==1 & a1==de.a1a2[1] & a2==de.a1a2[2], Y.a]/
        res[a0==0 & a1==de.a1a2[1] & a2==de.a1a2[2], Y.a],
      "ie1.rr"=
        res[a0==ie1.a0a2[1] & a1==1 & a2==ie1.a0a2[2], Y.a]/
        res[a0==ie1.a0a2[1] & a1==0 & a2==ie1.a0a2[2], Y.a],
      "ie2.rr"=
        res[a0==ie2.a0a1[1] & a1==ie2.a0a1[2] & a2==1, Y.a]/
        res[a0==ie2.a0a1[1] & a1==ie2.a0a1[2] & a2==0, Y.a],
      "te.rr"=
        res[a0==1 & a1==1 & a2==1, Y.a]/res[a0==0 & a1==0 & a2==0, Y.a])
  }
  DECOMP <- list(
    OneDecomposition(c(0,0),c(1,1),c(1,0),res),
    OneDecomposition(c(0,0),c(1,0),c(1,1),res),
    OneDecomposition(c(1,0),c(0,0),c(1,1),res),
    OneDecomposition(c(1,1),c(0,0),c(0,1),res),
    OneDecomposition(c(1,1),c(0,1),c(0,0),res),
    OneDecomposition(c(0,1),c(1,1),c(0,0),res)
  )
  # proportion mediated based on log of RR
  for (dd in 1:6) {
    ide.dd <- DECOMP[[dd]]
    pm.dd <- c(log(ide.dd["ie1.rr"]),
               log(ide.dd["ie2.rr"]),
               log(ide.dd["de.rr"]))
    pm.dd <- pm.dd/log(ide.dd["te.rr"])
    names(pm.dd) <- paste0("propmed.",c("M1","M2","DE"))
    DECOMP[[dd]] <- c(ide.dd,pm.dd*100)
  }
  return(DECOMP)
}
