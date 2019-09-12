############# LIBRARIES ################
library(ks)
library(plyr)
library(ggplot2)
library(MRS)
source("~/TEAM.R")


##### SIMULATION FUNCTIONS #####

#Params input determines the setting (local shift, dispersion, or shift+dispersion)
#See corresponding paper

run_sim <- function(n1,n0,theta0,params,K,alpha.vec,L,cutoff.val,seed){
  
  TEAM.Sim = array(NA,c(1,L,length(alpha.vec),4))
  
  set.seed(seed)
  X1 = rnorm.mixt(n=n1,mus=params$mu1,sigmas=params$sigma1,props=params$prop1)
  X0 = rnorm.mixt(n=n0,mus=params$mu0,sigmas=params$sigma0,props=params$prop0)
  
  m = 2^K
  
  #Performance
  Omega = seq(m)
  truth = get.truth(dat=c(X0,X1),m=2^K,cutoff=cutoff.val,
                    func=eff.size,params=params)
  
  
  nonnull.bins = truth$nonnull
  null.bins = setdiff(Omega,nonnull.bins)
  
  #Stratify true nonnulls by log OR
  
  for (i in seq(length(alpha.vec))){
    
    if(missing(theta0)){
      res = TEAM(x1=X1,x0=X0,K=K,alpha=alpha.vec[i],L=L)
    }
    else{
      res = TEAM(x1=X1,x0=X0,theta0=theta0,K=K,alpha=alpha.vec[i],L=L)
    }
  
      #Performance of TEAM
      for (l in 1:L){
        TEAM.Sim[1,l,i,1] = sum(res$S.list[[l]]%in%nonnull.bins) #TP
        TEAM.Sim[1,l,i,2] = sum(res$S.list[[l]]%in%null.bins) #FP
        TEAM.Sim[1,l,i,3] = sum(setdiff(Omega,res$S.list[[l]])%in%null.bins) #TN
        TEAM.Sim[1,l,i,4] = sum(setdiff(Omega,res$S.list[[l]])%in%nonnull.bins) #FN
      }
  }
  
  return(TEAM.Sim)
  
}

run_sim_strat <- function(params,theta0,K,alpha.vec,L,stratbrks,seed){
  TEAM.Sim.strat = array(NA,c(1,L,numstrats,length(alpha.vec),4))

  set.seed(seed)
  X1 = rnorm.mixt(n=params$n1,mus=params$mu1,sigmas=params$sigma1,props=params$prop1)
  X0 = rnorm.mixt(n=params$n0,mus=params$mu0,sigmas=params$sigma0,props=params$prop0)
  
  m = 2^K
  
  #Performance
  Omega = seq(m)
  truth = get.truth(dat=c(X0,X1),m=2^K,cutoff=theta0,
                    func=eff.size,params=params)
  
  nonnull.bins = truth$nonnull
  null.bins = setdiff(Omega,nonnull.bins)
  
  #Stratify true nonnulls 
  nonnull.strat = list(NULL)
  for(s in seq(numstrats)){
    nonnull.strat[[s]] =  nonnull.bins[which(stratbrks[s] < truth$eff.sizes[nonnull.bins] & truth$eff.sizes[nonnull.bins] <= stratbrks[s+1])]
  }
  
  for (i in seq(length(alpha.vec))){
    
    res = TEAM(x1=X1,x0=X0,theta0=theta0,K=K,alpha=alpha.vec[i],L=L)
    
    #Performance stratified by effect size
    
      for(l in 1:L){
        for(s in 1:length(nonnull.strat)){
          TEAM.Sim.strat[1,l,s,i,1] = sum(res$S.list[[l]]%in%nonnull.strat[[s]]) #TP
          TEAM.Sim.strat[1,l,s,i,2] = sum(res$S.list[[l]]%in%null.bins) #FP
          TEAM.Sim.strat[1,l,s,i,3] = sum(setdiff(Omega,res$S.list[[l]])%in%null.bins) #TN
          TEAM.Sim.strat[1,l,s,i,4] = sum(setdiff(Omega,res$S.list[[l]])%in%nonnull.strat[[s]]) #FN
        }
      }
  }

  return(TEAM.Sim.strat)
  
}


run_sim_compare <- function(n1,n0,theta0,params,K,L,mrs.fdr.vec,cutoff.val,seed){
  
  #compare 3 methods: #SLM, TEAM, MRS
  Comp.Perf.Sim = array(NA,c(1,3,length(mrs.fdr.vec),4)) 
  
  
  set.seed(seed)
  X1 = rnorm.mixt(n=n1,mus=params$mu1,sigmas=params$sigma1,props=params$prop1)
  X0 = rnorm.mixt(n=n0,mus=params$mu0,sigmas=params$sigma0,props=params$prop0)
  
  m = 2^K
  
  #Performance
  Omega = seq(m)
  truth = get.truth(dat=c(X0,X1),m=2^K,cutoff=cutoff.val,
                    func=eff.size,params=params,two.sided=TRUE)

  nonnull.bins = truth$nonnull
  null.bins = setdiff(Omega,nonnull.bins)
  #emp.FDR.MRS = NULL
  
  for (i in seq(length(mrs.fdr.vec))){
    
    #Run MRS on two-sided alternative and obtain empirical frequentist FDR (on last L layers)
    #MRS
    X = c(X1,X0)
    G = rep(c(2,1),times=c(n1,n0))
    data = list(X,G)
    names(data) = c("X","G")
    
    ans = mrs(data$X,data$G,K=K)
    
    #MRS Performance
    t = FDRThreshold(ans=ans,fdr = mrs.fdr.vec[i],by = 0.1)
    
    subset.ind = intersect(which(ans$RepresentativeTree$Levels%in%seq(K-L+1,K)),
                           which(ans$RepresentativeTree$EffectSizes[,1]>t))
  
    
    MRS.endpts = ans$RepresentativeTree$Regions[subset.ind,]
    
    res.MRS = NULL
    if(length(subset.ind)>0){
      if(is.matrix(MRS.endpts)){
        res.MRS = which(apply(apply(truth$region.endpts,1,
                                    function(x) (x[1]<=MRS.endpts[,2])&(x[2]>=MRS.endpts[,1])),2,any))
      }
      else{
        res.MRS = which(apply(truth$region.endpts,1,
                              function(x) (x[1]<=MRS.endpts[2])&(x[2]>=MRS.endpts[1])),2,any)
      }
    }
    
    #MRS
    Comp.Perf.Sim[1,3,i,1] = sum(res.MRS%in%nonnull.bins) #TP
    Comp.Perf.Sim[1,3,i,2] = sum(res.MRS%in%null.bins) #FP
    Comp.Perf.Sim[1,3,i,3] = sum(setdiff(Omega,res.MRS)%in%null.bins) #TN
    Comp.Perf.Sim[1,3,i,4] = sum(setdiff(Omega,res.MRS)%in%nonnull.bins) #FN
    
    emp.FDR.MRS = Comp.Perf.Sim[1,3,i,2]/max(length(res.MRS),1)
    
    #Two-sided alternatives
    res.gt = TEAM(x1=X1,x0=X0,theta0=theta0,K=K,alpha=emp.FDR.MRS,L=L)
    res.lt = TEAM(x1=X0,x0=X1,theta0=theta0,K=K,alpha=emp.FDR.MRS,L=L)
    
    res = list("n"=res.gt$n,
               "m"=res.gt$m,
               "m.excl"=union(res.gt$m.excl,res.lt$m.excl),
               "S.list"=do.call(Map,c(union,list(res.gt$S.list,res.lt$S.list)))
    )
    

    #SLM
    Comp.Perf.Sim[1,1,i,1] = sum(res$S.list[[1]]%in%nonnull.bins) #TP
    Comp.Perf.Sim[1,1,i,2] = sum(res$S.list[[1]]%in%null.bins) #FP
    Comp.Perf.Sim[1,1,i,3] = sum(setdiff(Omega,res$S.list[[1]])%in%null.bins) #TN
    Comp.Perf.Sim[1,1,i,4] = sum(setdiff(Omega,res$S.list[[1]])%in%nonnull.bins) #FN
    
    #TEAM
    Comp.Perf.Sim[1,2,i,1] = sum(res$S.list[[L]]%in%nonnull.bins) #TP
    Comp.Perf.Sim[1,2,i,2] = sum(res$S.list[[L]]%in%null.bins) #FP
    Comp.Perf.Sim[1,2,i,3] = sum(setdiff(Omega,res$S.list[[L]])%in%null.bins) #TN
    Comp.Perf.Sim[1,2,i,4] = sum(setdiff(Omega,res$S.list[[L]])%in%nonnull.bins) #FN
  }
  
  return(Comp.Perf.Sim)
}

############# AUXILIARY SIMULATION FUNCTIONS ###############
FDRThreshold <- function(ans, fdr = 0.01, by = 0.1) {
  eff = ans$RepresentativeTree$EffectSizes[, 1]
  q = rev(seq(min(eff,na.rm = TRUE),
              max(eff,na.rm = TRUE),
              by = by))
  it = 0
  emp.fdr = 0
  while (emp.fdr < fdr) {
    it = it + 1
    D = which(eff > q[it])
    emp.fdr = ifelse(length(D) > 0,
                     sum(1 - ans$RepresentativeTree$AltProbs[D])/length(D),
                     0)
  }
  return(q[it])
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
  Power = apply(sim.mat[,,1]/max
                (sim.mat[,,1]+sim.mat[,,4],1),2,mean,na.rm=TRUE)
  nonnulls = mean(sim.mat[,,1]+sim.mat[,,4])
  FPR = apply(sim.mat[,,2]/max(sim.mat[,,2]+sim.mat[,,3],1),2,mean,na.rm=TRUE)
  return(list("num.rej"=num.rej,"TD"=TD,"FN"=FN,"FD"=FD,"FNP"=FNP,"FDP"=FDP,
              "Power"=Power,"FPR"=FPR,"nonnulls"=nonnulls))
}

pmnorm <- function(x, mu, sigma, pmix) {
  (1-pmix[1])*pnorm(x,mu[1],sigma[1]) + pmix[1]*pnorm(x,mu[2],sigma[2])
}

eff.size <- function(left,right,params){
  num = pmnorm(right,mu=params$mu1,sigma=params$sigma1,pmix=params$prop1[2]) - 
    pmnorm(left,mu=params$mu1,sigma=params$sigma1,pmix=params$prop1[2])
  denom = pmnorm(right,mu=params$mu0,sigma=params$sigma0,pmix=params$prop0[2]) - 
    pmnorm(left,mu=params$mu0,sigma=params$sigma0,pmix=params$prop0[2])
  params$n1*num/(params$n1*num+n0*denom)
}

get_effsize_strats <- function(n1,n0,params,K,cutoff,seed){
  
  eff.size.arr = array(NA,c(1,2^K))
  
  set.seed(seed)
  X1 = rnorm.mixt(n=n1,mus=params$mu1,sigmas=params$sigma1,props=params$prop1)
  X0 = rnorm.mixt(n=n0,mus=params$mu0,sigmas=params$sigma0,props=params$prop0)
  
  #Performance
  truth = get.truth(dat=c(X0,X1),m=2^K,cutoff=cutoff,
                    func=eff.size,params=params)
  
  nonnull.bins = truth$nonnull
  
  eff.size.arr[1,] = truth$eff.sizes
  return(eff.size.arr)
}


