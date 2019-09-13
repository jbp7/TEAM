###################################################################################
#
# Authors: J. Pura, C. Chan, PhD, J. Xie, PhD
#
# PURPOSE:
# Contains functions to run the simulations in the main text of
# Pura, Chan, Xie (2019). Multiple Hypothesis Testing on an Aggregation Tree Method 
#
###################################################################################

############# LIBRARIES ################
library(ks)
library(MRS)
source("~/TEAM.R")

run_simulation_layer <- function(nreps,params,L,K,alpha.vec){
  
  #Initialize output: nreps repetitions, L layers, number of FDR values, 4 performance metrics
  res_layer <- array(NA,c(nreps,L,length(alpha.vec),4))
  
  theta0 = params$n1/(params$n1+params$n0)
  
  for(b in 1:nreps){
    
    #Run TEAM by layer
    set.seed(b)
    X1 = rnorm.mixt(n=params$n1,mus=params$mu1,sigmas=params$sigma1,props=params$prop1)
    X0 = rnorm.mixt(n=params$n0,mus=params$mu0,sigmas=params$sigma0,props=params$prop0)
    
    m = 2^K
    
    #Performance
    Omega = seq(m)
    truth = get.truth(dat=c(X0,X1),m=2^K,cutoff=theta0,
                      func=eff.size,params=params)
    
    nonnull.bins = truth$nonnull
    null.bins = setdiff(Omega,nonnull.bins)
    
    for (i in seq_along(alpha.vec)){
      
      res = TEAM(x1=X1,x0=X0,theta0=theta0,K=K,L=L,alpha=alpha.vec[i])
      
      #Performance by layer
      for (l in 1:L){
        res_layer[1,l,i,1] = sum(res$S.list[[l]]%in%nonnull.bins) #TP
        res_layer[1,l,i,2] = sum(res$S.list[[l]]%in%null.bins) #FP
        res_layer[1,l,i,3] = sum(setdiff(Omega,res$S.list[[l]])%in%null.bins) #TN
        res_layer[1,l,i,4] = sum(setdiff(Omega,res$S.list[[l]])%in%nonnull.bins) #FN
      }
    }
    
  }
  
  return(res_layer)
}

run_simulation_strata <- function(nreps,stratbrks,params,L,K,alpha.vec){
  
  #Initialize output: nreps repetitions, L layers, number of strata, number of FDR values, 4 performance metrics
  
  res_strata <- array(NA,c(nreps,L,length(stratbrks)-1,length(alpha.vec),4)) 
  
  theta0 = params$n1/(params$n1+params$n0)
  
  for(b in 1:nreps){
    
    set.seed(b)
    X1 = rnorm.mixt(n=params$n1,mus=params$mu1,sigmas=params$sigma1,props=params$prop1)
    X0 = rnorm.mixt(n=params$n0,mus=params$mu0,sigmas=params$sigma0,props=params$prop0)
    
    m = 2^K
    theta0 = params$n1/(params$n1+params$n0)
    
    #Performance
    Omega = seq(m)
    truth = get.truth(dat=c(X0,X1),m=2^K,cutoff=theta0,
                      func=eff.size,params=params)
    
    nonnull.bins = truth$nonnull
    null.bins = setdiff(Omega,nonnull.bins)
    
    #Stratify true nonnulls 
    nonnull.strat = list(NULL)
    for(s in seq(length(stratbrks)-1)){
      nonnull.strat[[s]] =  nonnull.bins[which(stratbrks[s] < truth$eff.sizes[nonnull.bins] & truth$eff.sizes[nonnull.bins] <= stratbrks[s+1])]
    }
    
    for (i in seq(length(alpha.vec))){
      
      res = TEAM(x1=X1,x0=X0,theta0=theta0,K=K,L=L,alpha=alpha.vec[i])
      
      #Performance stratified by effect size
      for(l in 1:L){
        for(s in 1:length(nonnull.strat)){
          res_strata[b,l,s,i,1] = sum(res$S.list[[l]]%in%nonnull.strat[[s]]) #TP
          res_strata[b,l,s,i,2] = sum(res$S.list[[l]]%in%null.bins) #FP
          res_strata[b,l,s,i,3] = sum(setdiff(Omega,res$S.list[[l]])%in%null.bins) #TN
          res_strata[b,l,s,i,4] = sum(setdiff(Omega,res$S.list[[l]])%in%nonnull.strat[[s]]) #FN
        }
      }
    }
  
  }
  
  return(res_strata)
}


get.truth <- function(dat,m,cutoff,two.sided=FALSE,func,...){
  brk.x = quantile(dat,seq(0,1,length.out = m+1))
  endpts = matrix(sort(c(brk.x[-length(brk.x)],brk.x[-1])),ncol = 2,byrow = TRUE)
  eff.sizes = apply(endpts,1,function(x) func(x[1],x[2],...)) 
  
  if(two.sided){
    return(list("eff.sizes"=eff.sizes,"region.endpts" = endpts, "nonnull" = which(abs(eff.sizes) > cutoff)))
  }
  else{
    return(list("eff.sizes"=eff.sizes,"region.endpts" = endpts, "nonnull" = which(eff.sizes > cutoff)))
  }
}

evaluate.sim <- function(sim.mat){
  num.rej = apply(sim.mat[,,1]+sim.mat[,,2],2,mean,na.rm=TRUE)
  TD = apply(sim.mat[,,1],2,mean,na.rm=TRUE) #true discoveries
  FD = apply(sim.mat[,,2],2,mean,na.rm=TRUE)
  FDP = apply(sim.mat[,,2]/max(sim.mat[,,1]+sim.mat[,,2],1),2,mean,na.rm=TRUE)
  FN = apply(sim.mat[,,4],2,mean,na.rm=TRUE)
  FNP = apply(sim.mat[,,4]/max(sim.mat[,,3]+sim.mat[,,4],1),2,mean,na.rm=TRUE)
  Sensitivity = apply(sim.mat[,,1]/max
                (sim.mat[,,1]+sim.mat[,,4],1),2,mean,na.rm=TRUE)
  nonnulls = mean(sim.mat[,,1]+sim.mat[,,4])
  FPR = apply(sim.mat[,,2]/max(sim.mat[,,2]+sim.mat[,,3],1),2,mean,na.rm=TRUE)
  return(list("num.rej"=num.rej,"TD"=TD,"FN"=FN,"FD"=FD,"FNP"=FNP,"FDP"=FDP,
              "Sensitivity"=Sensitivity,"FPR"=FPR,"nonnulls"=nonnulls))
}

pmnorm <- function(x, mu, sigma, pmix) {
  (1-pmix[1])*pnorm(x,mu[1],sigma[1]) + pmix[1]*pnorm(x,mu[2],sigma[2])
}

eff.size <- function(left,right,params){
  num = pmnorm(right,mu=params$mu1,sigma=params$sigma1,pmix=params$prop1[2]) - 
    pmnorm(left,mu=params$mu1,sigma=params$sigma1,pmix=params$prop1[2])
  denom = pmnorm(right,mu=params$mu0,sigma=params$sigma0,pmix=params$prop0[2]) - 
    pmnorm(left,mu=params$mu0,sigma=params$sigma0,pmix=params$prop0[2])
  
  return(params$n1*num/(params$n1*num+params$n0*denom))
}



