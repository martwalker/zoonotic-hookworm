###########################################
## functions
###########################################
funcs <- list(ddred = function(n, b){
  (n)^(b-1) 
},

densdep = function(W, k, b){
  
  W<-max(W,pars["tol"])
  
  max <- max(qnbinom(0.999, size=k, mu=W),1)
  
  probs <- dnbinom(seq(1,max), size = k, mu= W)
  
  sumsol <- sum( funcs$ddred(1:max, b=b)*probs ) / (1-(1+W/k)^(-k))
  
  return(as.numeric(sumsol))
},


# mating probability function
matingf = function(W, k, rho) {
  W<-max(W,pars["tol"])
  1-(1+rho*W/k)^(-(k+1))
},

##  egg output is here parameterized in terms of adult worms
prev = function(thresh, W, k, a, b) {
  W<-max(W,pars["tol"])
  1-pnbinom((thresh/a)^(1/b), mu = W, size=k)
},

R=function(Rmat) {
  eigen(Rmat)$values[1]
},

partrans = function(pars) {
  with(as.list(c(pars)), {
  
  ## calculate intra- and inter-host contribution to transmission
  
  #omega22 <- 1/omega11 #omega11 = R011/R022; omega22 = R022/R011
  #R011 <- R0*omega11/(1+omega11)
  #R022 <- R0*omega22/(1+omega22)
  
  #omega12 <- 1/omega21 #omega21 = R021/R012; omega12 = R012/R021
  #R021 <- sqrt(R011*R022*omega21)
  #R012 <- sqrt(R011*R022*omega12)
    
  omega22 <- (1-omega11)
  omega12 = (1-omega21)
  R011 <- omega11*R0
  R022 <- omega22*R0
  R012 <- (omega11*omega22*R0^2)^omega12
  R021 <- (omega11*omega22*R0^2)^omega21
  
  Rmat <- matrix(c(R011, R012, R021, R022), nrow=2, byrow = T)
  
  if (round(det(Rmat), log10(1/tol))!=0) {
    stop ("transmission matrix not invertible")
  } else {
    return(Rmat)
  }
  }) 
  },
# event function for tx (dogs)
eventfun1 = function(t, y, pars) {
  with(as.list(y), {
    
    ymat <- matrix(y, ncol=2, nrow=1)
    
    #ymat[,1] <- max(ymat[,1], pars["tol"])
    #ymat[,2] <- max(ymat[,2], pars["tol"])
    
    if (pars["c1"]>0) {
    ymat[,1] <- ymat[,1]*(1-pars["epsilon1"]*pars["c1"])
    }
    
    ymat[,1] <- max(ymat[,1], pars["tol"])
    ymat[,2] <- max(ymat[,2], pars["tol"])
    
    return(c(ymat))
  })
  
},

# event function for tx (humans)
eventfun2 = function(t, y, pars) {
  with(as.list(y), {
    
    ymat <- matrix(y, ncol=2, nrow=1)
    
   
    
    if(pars["c2"]>0) {
    ymat[,2] <- ymat[,2]*(1-pars["epsilon2"]*pars["c2"])
    }
    
    ymat[,1] <- max(ymat[,1], pars["tol"])
    ymat[,2] <- max(ymat[,2], pars["tol"])
    
    return(c(ymat))
  })
  
},

# event function for tx (humans + dogs)
eventfun12 = function(t, y, pars) {
  with(as.list(y), {
    
    ymat <- matrix(y, ncol=2, nrow=1)
    
   
    
    if(pars["c1"]>0) {
    ymat[,1] <- ymat[,1]*(1-pars["epsilon1"]*pars["c1"])
    }
    if(pars["c2"]>0) {
    ymat[,2] <- ymat[,2]*(1-pars["epsilon2"]*pars["c2"])
    }
    
    ymat[,1] <- max(ymat[,1], pars["tol"])
    ymat[,2] <- max(ymat[,2], pars["tol"])
    
    return(c(ymat))
  })
  
},

# join the event functions together
eventfun = function(t, y, pars) {
  with(as.list(c(y,pars)), {
    t1 <- seq(start.tx, start.tx+(n.tx1-1)*freq.tx1, freq.tx1)
    t2 <- seq(start.tx, start.tx+(n.tx2-1)*freq.tx2, freq.tx2)
    
    ## both events 
    if (c1>0 & c2>0 & any(abs(t-t2) < dt) & any(abs(t-t1) < dt)) {funcs$eventfun12(t,y,pars)} 
    else if (c1>0 & c2==0 & any(abs(t-t1) < dt)) {funcs$eventfun1(t,y,pars)} 
    else if (c2>0 & c1==0 & any(abs(t-t2) < dt)) {funcs$eventfun2(t,y,pars)} 
    else y
  })
},

rootfunbase = function(t, y, pars) {
  with(as.list(y,pars), {
    ymat<- matrix(y, ncol=2)
    return(c(ymat[,1] - pars["tol"], ymat[,2]-pars["tol"]))
    
  })
},

rootfunint = function(t, y, pars) {
  with(as.list(y, pars),{
   
    ## calculate intra- and inter-host contribution to transmission
    Rmat <- funcs$partrans(pars=pars)
    R011 <- Rmat[1,1]
    R022 <- Rmat[2,2]
    R012 <- Rmat[1,2]
    R021 <- Rmat[2,1]
    
    ymat<- matrix(y, ncol=2)
    
    ymat[,1] <- max(ymat[,1], pars["tol"])
    ymat[,2] <- max(ymat[,2], pars["tol"])
    
    W1 <- ymat[,1]; W2<- ymat[,2] 
    
    # dynamic overdispersion parameter for the distribution of worms (NBD)
    kdyn1 <- as.vector(e$kdyn1)
    kdyn2 <- as.vector(e$kdyn2)
    kpost1 <- as.vector(e$kpost1)
    kpost2 <- as.vector(e$kpost2)
    kinf1 <- as.vector(e$kinf1)
    kinf2 <- as.vector(e$kinf2)
    Wpost1 <- as.vector(e$Wpost1)
    Wpost2 <- as.vector(e$Wpost2)
    Winf1 <- as.vector(e$Winf1)
    Winf2 <- as.vector(e$Winf2)
    
    # mating probability
    mp1 <- funcs$matingf(W=W1,k=kdyn1,rho=pars["rho"])
    mp2 <- funcs$matingf(W=W2,k=kdyn2,rho=pars["rho"])
    
    # density dependent fecundity 
    dd1 <- funcs$densdep(W=W1, k=kdyn1, b=pars["b"])
    dd2 <- funcs$densdep(W=W2, k=kdyn2, b=pars["b"])
    
    # effective reproduction numbers recipient i donor j
    RE11 <- R011*dd1*mp1
    RE22 <- R022*dd2*mp2
    
    RE12 <- R012*dd2*mp2
    RE21 <- R021*dd1*mp1
    
    REmat <- matrix(c(RE11, RE12, RE21, RE22), nrow=2, byrow = T)
    
    RE <- funcs$R(REmat)
    
    #print(RE)
    
    return(rnorm(1,sign(RE-0.99)+1, sd=pars["tol"]))
  })
},

# model equations
mod = function(t, y, pars) {
  with(as.list(c(y,pars)), {
    
    ## calculate intra- and inter-host contribution to transmission
    Rmat <- funcs$partrans(pars=pars)
    R011 <- Rmat[1,1]
    R022 <- Rmat[2,2]
    R012 <- Rmat[1,2]
    R021 <- Rmat[2,1]

    ymat<-dymat<-matrix(y, ncol=2, nrow=1) 
    
    ## prevents numerical errors when worm numbers very low
    ymat[,1]<-max(ymat[,1], pars["tol"])
    ymat[,2] <- max(ymat[,2], pars["tol"])
    
    # dog and humans mean female worm burden
    W1 <- ymat[,1]; W2<- ymat[,2] 
    
    # dynamic overdispersion parameter for the distribution of worms (NBD)
    kdyn1 <- as.vector(e$kdyn1)
    kdyn2 <- as.vector(e$kdyn2)
    kpost1 <- as.vector(e$kpost1)
    kpost2 <- as.vector(e$kpost2)
    kinf1 <- as.vector(e$kinf1)
    kinf2 <- as.vector(e$kinf2)
    Wpost1 <- as.vector(e$Wpost1)
    Wpost2 <- as.vector(e$Wpost2)
    Winf1 <- as.vector(e$Winf1)
    Winf2 <- as.vector(e$Winf2)
    
    # mating probability
    mp1 <- funcs$matingf(W=W1,k=kdyn1,rho=rho)
    mp2 <- funcs$matingf(W=W2,k=kdyn2,rho=rho)
    
    # density dependent fecundity 
    dd1 <- funcs$densdep(W=W1, k=kdyn1, b=b)
    dd2 <- funcs$densdep(W=W2, k=kdyn2, b=b)
    
    # prevalence of infection (in faecal sample, i.e. epg or gentic material)
    Wp1 <- funcs$prev(thresh=0, W=W1, k=kdyn1,a=a, b=b)
    Wp2 <- funcs$prev(thresh=0, W=W2, k=kdyn2, a=a, b=b)
    
    # prevalence of moderate to heavy infection (in faecal sample, i.e. epg or gentic material)
    Wph1 <- funcs$prev(thresh=z, W=W1, k=kdyn1, a=a, b=b)
    Wph2 <- funcs$prev(thresh=z, W=W2, k=kdyn2, a=a, b=b)
    
    # epg
    epg1 <- W1*a*dd1*mp1
    epg2 <- W2*a*dd2*mp2
    
    # effective reproduction numbers recipient i donor j
    RE11 <- R011*dd1*mp1
    RE22 <- R022*dd2*mp2
    
    RE12 <- R012*dd2*mp2
    RE21 <- R021*dd1*mp1
    
    REmat <- matrix(c(RE11, RE12, RE21, RE22), nrow=2, byrow=T)
    
    RE <- funcs$R(REmat)
    
    # ODEs
    #dymat[,1] <- muW*(RE11+RE12) - (muW+mu1)*W1
    #dymat[,2] <- muW*(RE22+RE21) - (muW+mu2)*W2
    
    dymat[,1] <- (muW+mu1)*RE11*W1+(muW+mu2)*RE12*W2 - (muW+mu1)*W1
    dymat[,2] <- (muW+mu1)*RE21*W1+(muW+mu2)*RE22*W2 - (muW+mu2)*W2
    

    if (pars["equib"]!=1) {
      ## treatment time?
      tmp <- funcs$eventfun(t,c(W1,W2),pars=pars)
      if (tmp[1]!=W1 & pars["c1"]>0) {
        kdyn1 <- kinf1*W1/( (1+kinf1)*Winf1 - kinf1*W1)
        Wpost1 <- W1
        kpost1 <- kdyn1
      }   else if (t>pars["start.tx"] & pars["c1"]>0) {
        kdyn1 <- W1^2*(Winf1-Wpost1)^2/
          ((Winf1^2/kinf1)*(W1-Wpost1)^2 + (Wpost1^2/kpost1)*(W1-Winf1)^2)
      }
      if (tmp[2]!=W2 & pars["c2"]>0) {
        kdyn2 <- kinf2*W2/( (1+kinf2)*Winf2 - kinf2*W2)
        Wpost2 <- W2
        kpost2 <- kdyn2
      } else if (t>pars["start.tx"] & pars["c2"]>0) {
        kdyn2 <- W2^2*(Winf2-Wpost2)^2/
          ((Winf2^2/kinf2)*(W2-Wpost2)^2 + (Wpost2^2/kpost2)*(W2-Winf2)^2)
      } else {
        kdyn1<-kdyn1
        kdyn2<-kdyn2
      }
      
    }
    
    ## store updates in environment
    e$kdyn1 <<- kdyn1
    e$kdyn2 <<- kdyn2
    e$kinf1 <<- kinf1
    e$kinf2 <<- kinf2
    e$kpost1 <<- kpost1
    e$kpost2 <<- kpost2
    e$Wpost1 <<- Wpost1
    e$Wpost2 <<- Wpost2
    e$Winf1 <<- Winf1 
    e$Winf2 <<- Winf2 
    
    #print(c(W1, W2, RE))
    
    return(list(y=(rbind(dymat)),
                W1=W1, W2=W2, Wp1=Wp1, Wp2=Wp2, Wph1=Wph1, Wph2=Wph2,
                epg1=epg1, epg2=epg2, kdyn1=kdyn1, kdyn2=kdyn2,
                RE11=RE11, RE12=RE12, RE21=RE21, RE22=RE22,
                RE=RE, R011=R011, R012=R012, R021=R021, R022=R022))
    
  })
},

runmod = function(pars, inits=c(5,5)) {
  

  ## initialize dynamic variable environment
  e <<- new.env()
  e$kdyn1 <- as.vector(pars["k1"])
  e$kdyn2 <- as.vector(pars["k2"])
  e$kinf1 <- as.vector(pars["k1"])
  e$kinf2 <- as.vector(pars["k2"])
  e$kpost1 <- as.vector(pars["k1"])
  e$kpost2 <- as.vector(pars["k2"])
  e$Winf1 <- inits[1]
  e$Winf2 <- inits[2]
  e$Wpost1 <- inits[1]
  e$Wpost2 <- inits[2]
  
  with(as.list(c(pars)), {
    
    ## calculate intra- and inter-host contribution to transmission
    Rmat <- funcs$partrans(pars=pars)
    R011 <- Rmat[1,1]
    R022 <- Rmat[2,2]
    R012 <- Rmat[1,2]
    R021 <- Rmat[2,1]

    y0 <- inits

    if (equib==1) {
      t1 <- seq(0,start.tx,by=dt)
      
      out<-lsoda(y=y0, times=t1, func=funcs$mod, parms=pars, rootfun=funcs$rootfunbase)
      
    } else {
      
      stop.tx <- max(c(start.tx+(n.tx1-1)*freq.tx1, start.tx+(n.tx2-1)*freq.tx2  ))  
      min.freq <- min(c(freq.tx1, freq.tx2))
      
      t1 <- seq(0,stop.t,by=dt)
      
      out<-lsoda(y=y0, times=t1, func=funcs$mod, parms=pars,
                 events = list(func=funcs$eventfun, 
                               time=seq(start.tx,stop.tx,min.freq), parms=pars, root=F))
    }
    
    out <- as.data.frame(out)
    out <- out[,-c(2,3)]
    
    every <- floor(dtout/dt)
    
    out[seq(1,nrow(out), every),]
    
  })
},
findstable2 = function(par) {
  maxW <- max(as.numeric((parameters["R0"]/parameters["rho"]*parameters["mu2"]-parameters["mu2"])/
                           (parameters["R0hs"]*parameters["mu1"]*parameters["N1N2"])),1)
  {
    W <- maxW-cntrlpar$accendem/2
    while(I(W - cntrlpar$accendem)>0) {
      W0 <- W
      W1 <- W0-cntrlpar$accendem
      dW0 <- funcs$derivs(W0, par=parameters)
      dW1 <- funcs$derivs(W1, par=parameters)
      i <- diff(sign(c(dW0, dW1)))
      if(i!=0) {
        return(W0-(W1-W0)/2)
        break
      }
      W <- W - cntrlpar$accendem
      #print(W)
      
    }
    return(NA)
  }
},

findbreak2 = function(par, Wstar)
{
  if(is.na(Wstar)) {
    return(NA)
  } else {
    W<-1.0E-12
    while(W<Wstar)
    {
      W0 <- W
      W1 <- W0+cntrlpar$accbreak
      dW0 <- funcs$derivs(W0, par=parameters)
      dW1 <- funcs$derivs(W1, par=parameters)
      i <- diff(sign(c(dW0, dW1)))
      if(i!=0) {
        return(W0+(W1-W0)/2)
        break
      }
      W <- W + cntrlpar$accbreak
    }
  }
  return(NA)
}, 

plotit = function(pars, inits=c(5,5)) {
  tmp <- funcs$runmod(pars,inits=inits)
  par(mfrow=c(1,3))
  plot(tmp[,"Wp1"]~tmp[,"time"], ylim=c(0,1), type="l")
  lines(tmp[,"Wp2"]~tmp[,"time"], lty="dashed")
  plot(tmp[,"W1"]~tmp[,"time"], type="l", ylim=c(0, max(c(tmp[,"W1"], tmp[,"W2"]))))
  lines(tmp[,"W2"]~tmp[,"time"], lty="dashed")
  plot(tmp[,"RE"]~tmp[,"time"], type="l")
},

PRCC=function(outcome,covariates){
  
  rank_outcome<-rank(outcome)
  rank_covariates<-as.data.frame(apply(covariates, 2, rank))
  
  PRCC_out<-c()
  for (par in 1:ncol(covariates)){
    xx<-rank_covariates[,par]		  
    xy<-rank_covariates[,-par]	
    
    xj<-xx-predict(lm(xx~.,data=xy))
    yy<-rank_outcome-predict(lm(rank_outcome~.,data=xy))
    
    PRCC_out[par]<-cov(xj,yy)/(sqrt(var(xj)*var(yy)))	
  }
  names(PRCC_out)<-colnames(covariates)
  PRCC_out
}


)
