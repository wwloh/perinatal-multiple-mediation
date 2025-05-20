# helper function to create duplicated data for one individual
Dupdata <- function(t) {
  # all possible combinations of hypothetical treatments
  out <- expand.grid(lapply(1:(t+1L), function(a.t) 0:1))
  colnames(out) <- paste0("a",0:t)
  out <- cbind("a.i"=1:nrow(out),out)
  return(out)
}

OneMCestimator <- function(Data,mc_draws=100,fit_MY,
                           Aname,Mnames,Yname,Cnames,Lidx=NULL) {
  # number of distinct mediators
  t <- length(Mnames)-length(Lidx)
  # indices for mediator(s)
  Midx <- seq_len(length(Mnames))
  if (!is.null(Lidx)) {
    Midx <- Midx[-Lidx]
  }
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
    
    # newdata for imputation: same observed C for all duplicated rows
    newdataMs <- data.frame(mydt_mc[, ..Cnames])
    
    # initialize counterfactual time-dependent confounder or mediator draws
    for (s in c(Lidx,Midx)) {
      for (aa in 0:1) {
        fitMs_Aa <- fit_MY[[Mnames[s]]]
        # relevant rows in duplicated data
        if (s %in% Lidx) {
          ## hypothetical exposure levels are same as direct effect's 
          s_aa <- mydt_mc[, get("a0")]==aa
        } else if (s %in% Midx) {
          s_aa <- mydt_mc[, get(paste0("a",which(s==Midx)))]==aa
        }
        # values of predictors for counterfactual outcomes
        newdataMs.aa <- cbind(aa,newdataMs[s_aa,])
        colnames(newdataMs.aa)[1] <- Aname
        if (s %in% Lidx && s>Lidx[1]) {
          # conditional on draws of previous time-dependent confounders
          otherLname <- Mnames[Lidx[1]:(which(s==Lidx)-1L)]
          newdataMs.aa <- cbind(
            newdataMs.aa,
            mydt_mc[s_aa, paste0(otherLname,".a"), with=FALSE])
          # change column name to same as original variable name
          colnames(newdataMs.aa)[
            match(paste0(otherLname,".a"),colnames(newdataMs.aa))] <- otherLname
          rm(otherLname)
        }
        if (!is.null(Lidx) && s %in% Midx) {
          # conditional on draws of all time-dependent confounders
          ## based on hypothetical exposure level of indirect effect
          newdataMs.aa <- cbind(
            newdataMs.aa,
            mydt_mc[mydt_mc[, get("a0")]==aa, 
                    paste0(Mnames[Lidx],".a"), with=FALSE])
          # change column names to same as original variable names
          colnames(newdataMs.aa)[
            match(paste0(Mnames[Lidx],".a"),colnames(newdataMs.aa))] <- Mnames[Lidx]
        }
        # random draw of counterfactual TDC or mediator
        if (class(fitMs_Aa)[1]=="glm") {
          # binomial distribution
          meanMs_Aa <- predict.glm(
            object=fitMs_Aa,newdata=newdataMs.aa,type="response")
          drawnMs_Aa <- rbinom(n=length(meanMs_Aa),size=1,prob=meanMs_Aa)
        } else if (class(fitMs_Aa)[1]=="multinom") {
          # multinomial distribution
          meanMs_Aa <- predict(object=fitMs_Aa,newdata=newdataMs.aa,type="probs")
          drawnMs_Aa <- apply(meanMs_Aa,1,stats::rmultinom,n=1,size=1)
          drawnMs_Aa <- as.integer(apply(drawnMs_Aa,2,function(dM)
            colnames(meanMs_Aa)[which(dM==1L)]))
        }
        mydt_mc[s_aa, paste0(Mnames[s],".a") := drawnMs_Aa]
        rm(drawnMs_Aa,meanMs_Aa,newdataMs.aa,s_aa)
      }
    }
    setkey(mydt_mc)
    
    # helper function to predict Y for different outcome models
    PredictY <- function(onedat) {
      setnames(onedat,"a0",Aname)
      Ya <- rep(NA_real_, nrow(onedat))
      onedat.df <- data.frame(onedat)
      # categorical mediators: set factors with same levels as in outcome model
      for (s in c(Lidx,Midx)) {
        Ms.levels <- levels(fit_MY[[Yname]]$model[,Mnames[s]])
        if (!is.null(Ms.levels)) {
          onedat.df[, Mnames[s]] <- factor(
            onedat.df[, Mnames[s]],levels=Ms.levels)
        }
        rm(Ms.levels)
      }
      check.am <- tryCatch(system.time(
        pred.Ya <- predict.glm(
          object=fit_MY[[Yname]],
          newdata=onedat.df,
          type="response")
      )[3], error=function(cond) return(NA))
      if (!is.na(check.am)) {
        Ya <- pred.Ya  
      }
      return(Ya)
    }
    
    ## predict potential outcomes using randomly sampled mediator values
    mydt_mc[, eval(Mnames) := NULL]
    setnames(mydt_mc,paste0(Mnames,".a"),Mnames)
    setkey(mydt_mc)
    mydt_mc[, "Y.a" := PredictY(
      onedat=mydt_mc[, c("a0",Cnames,Mnames), with=FALSE])]
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

IEestimator_TVC <- function(res, Mnames, Lidx) {
  # Lidx = indices in Mnames for time-varying confounders
  res <- data.table(res)
  setkey(res)
  # number of distinct mediators
  t <- length(Mnames)
  Midx <- 1:t
  Midx <- Midx[!(Midx %in% Lidx)]
  ## retain only rows with the same treatment levels as direct effect
  res.noL <- res[,lapply(.SD,"==",a0),.SDcols=paste0("a",Lidx),
                 by=c("a.i",paste0("a",0:t))]
  setnames(res.noL,
           old=ncol(res.noL)-((length(Lidx)-1):0),
           new=paste0("a",Lidx,"_eq_a0"))
  setkey(res.noL)
  res <- res[apply(res.noL[,paste0("a",Lidx,"_eq_a0"),with=FALSE],1,all)]
  setkey(res)
  
  res[, paste0("a",Lidx) := NULL]
  setnames(res, paste0("a",Midx), paste0("a",1:length(Midx)))
  
  return( res )
}

PROP_MEDIATED_logRR_4M <- function(res) {
  res <- data.table(res)
  setkey(res)
  OneDecomposition <- function(de.a1a2a3a4,ie1.a0a2a3a4,ie2.a0a1a3a4,
                               ie3.a0a1a2a4,ie4.a0a1a2a3,res){
    c("de.rr"=
        res[a0==1 & a1==de.a1a2a3a4[1] & a2==de.a1a2a3a4[2] &
              a3==de.a1a2a3a4[3] & a4==de.a1a2a3a4[4], Y.a]/
        res[a0==0 & a1==de.a1a2a3a4[1] & a2==de.a1a2a3a4[2] &
              a3==de.a1a2a3a4[3] & a4==de.a1a2a3a4[4], Y.a],
      "ie1.rr"=
        res[a0==ie1.a0a2a3a4[1] & a1==1 & a2==ie1.a0a2a3a4[2] &
              a3==ie1.a0a2a3a4[3] & a4==ie1.a0a2a3a4[4], Y.a]/
        res[a0==ie1.a0a2a3a4[1] & a1==0 & a2==ie1.a0a2a3a4[2] &
              a3==ie1.a0a2a3a4[3] & a4==ie1.a0a2a3a4[4], Y.a],
      "ie2.rr"=
        res[a0==ie2.a0a1a3a4[1] & a1==ie2.a0a1a3a4[2] & a2==1 & 
              a3==ie2.a0a1a3a4[3] & a4==ie2.a0a1a3a4[4], Y.a]/
        res[a0==ie2.a0a1a3a4[1] & a1==ie2.a0a1a3a4[2] & a2==0 & 
              a3==ie2.a0a1a3a4[3] & a4==ie2.a0a1a3a4[4], Y.a],
      "ie3.rr"=
        res[a0==ie3.a0a1a2a4[1] & a1==ie3.a0a1a2a4[2] & a2==ie3.a0a1a2a4[3] & 
              a3==1 & a4==ie3.a0a1a2a4[4], Y.a]/
        res[a0==ie3.a0a1a2a4[1] & a1==ie3.a0a1a2a4[2] & a2==ie3.a0a1a2a4[3] & 
              a3==0 & a4==ie3.a0a1a2a4[4], Y.a],
      "ie4.rr"=
        res[a0==ie4.a0a1a2a3[1] & a1==ie4.a0a1a2a3[2] & a2==ie4.a0a1a2a3[3] & 
              a3==ie4.a0a1a2a3[4] & a4==1, Y.a]/
        res[a0==ie4.a0a1a2a3[1] & a1==ie4.a0a1a2a3[2] & a2==ie4.a0a1a2a3[3] & 
              a3==ie4.a0a1a2a3[4] & a4==0, Y.a],
      "te.rr"=
        res[a0==1 & a1==1 & a2==1 & a3==1 & a4==1, Y.a]/
        res[a0==0 & a1==0 & a2==0 & a3==0 & a4==0, Y.a])
  }
  DECOMP <- list(
    OneDecomposition(c(0,0,0,0),c(1,0,0,0),c(1,1,0,0),c(1,1,1,0),c(1,1,1,1),
                     res),
    OneDecomposition(c(1,1,1,1),c(0,1,1,1),c(0,0,1,1),c(0,0,0,1),c(0,0,0,0),
                     res)
  )
  # proportion mediated based on log of RR
  for (dd in 1:length(DECOMP)) {
    ide.dd <- DECOMP[[dd]]
    pm.dd <- c(log(ide.dd["ie1.rr"]),
               log(ide.dd["ie2.rr"]),
               log(ide.dd["ie3.rr"]),
               log(ide.dd["ie4.rr"]),
               log(ide.dd["de.rr"]))
    pm.dd <- pm.dd/log(ide.dd["te.rr"])
    names(pm.dd) <- paste0("propmed.",c("M1","M2","M3","M4","DE"))
    DECOMP[[dd]] <- c(ide.dd,pm.dd*100)
  }
  return(DECOMP)
}

PROP_MEDIATED_logRR_singleM <- function(res) {
  res <- data.table(res)
  setkey(res)
  OneDecomposition <- function(de.a1,ie1.a0,res){
    c("de.rr"=
        res[a0==1         & a1==de.a1[1], Y.a]/
        res[a0==0         & a1==de.a1[1], Y.a],
      "ie1.rr"=
        res[a0==ie1.a0[1] & a1==1, Y.a]/
        res[a0==ie1.a0[1] & a1==0, Y.a],
      "te.rr"=
        res[a0==1         & a1==1, Y.a]/
        res[a0==0         & a1==0, Y.a])
  }
  DECOMP <- list(
    OneDecomposition(c(0),c(1),res),
    OneDecomposition(c(1),c(0),res)
  )
  # proportion mediated based on log of RR
  for (dd in 1:length(DECOMP)) {
    ide.dd <- DECOMP[[dd]]
    pm.dd <- c(log(ide.dd["ie1.rr"]),
               log(ide.dd["de.rr"]))
    pm.dd <- pm.dd/log(ide.dd["te.rr"])
    names(pm.dd) <- paste0("propmed.",c("M1","DE"))
    DECOMP[[dd]] <- c(ide.dd,pm.dd*100)
  }
  return(DECOMP)
}
