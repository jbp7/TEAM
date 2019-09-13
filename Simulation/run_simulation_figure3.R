###################################################################################
#
# Authors: J. Pura, C. Chan, PhD, J. Xie, PhD
#
# PURPOSE:
# Generates Figure 3 in the main text of
# Pura, Chan, Xie (2019). Multiple Hypothesis Testing on an Aggregation Tree Method 
#
# DETAILS:
# This code generates the data, runs TEAM, and plots the avg Sensitivity
# across layers and FDR for each of the three local settings: shift, dispersion, 
# and shift+dispersion. The 25 quantiles seen in Figure 3 are fixed beforehand.
# For each setting, we generate a single random draw from that setting's 
# corresponding normal mixture and compute the quantiles. Even though the quantiles
# are random, they do not change much across repetitions.
#
###################################################################################

library(ggplot2)
library(cowplot)
library(zoo)
library(grid)

source("~/simulation_functions.R")

# Number of repetitions
nreps=1000

# Parameters for TEAM
L = 5
K = 14
alpha.vec = 0.01*c(5,10,15,20)

# Dataset 
n = 180
n1 = n0 = 0.5*n*2^K

theta0 = n1/(n1+n0)

# Parameters for each of the three settings
params.shift = list("n0"=n0, "n1"=n1,
                    "mu1"=c(0.2,0.89),"mu0"=c(0.2,0.88),
                    "sigma1"=c(0.04,0.01),"sigma0"=c(0.04,0.01),
                    "prop1"=c(0.97,0.03),"prop0"=c(0.97,0.03))

params.disp = list("n0"=n0, "n1"=n1,
                   "mu1"=c(0.4,0.8),"mu0"=c(0.4,0.8),
                   "sigma1"=c(0.04,0.03),"sigma0"=c(0.04,0.02),
                   "prop1"=c(0.97,0.03),"prop0"=c(0.97,0.03))

params.combo = list("n0"=n0, "n1"=n1,
                    "mu1"=c(0.4,0.82),"mu0"=c(0.4,0.8),
                    "sigma1"=c(0.04,0.05),"sigma0"=c(0.04,0.04),
                    "prop1"=c(0.97,0.03),"prop0"=c(0.98,0.02))

############# Run the simulations #############

#### Local Shift

#Generate strata based on a random draw from distribution and 
#fix the strata for nreps repetitions
start_time <- Sys.time()
set.seed(1)
X1 = rnorm.mixt(n=params.shift$n1,
                mus=params.shift$mu1,
                sigmas=params.shift$sigma1,
                props=params.shift$prop1)

X0 = rnorm.mixt(n=params.shift$n0,
                mus=params.shift$mu0,
                sigmas=params.shift$sigma0,
                props=params.shift$prop0)


truth.shift = get.truth(dat=c(X0,X1),m=2^K,cutoff=theta0,
                  func=eff.size,params=params.shift)

#specify 25 quantiles
stratbrks.shift = quantile(truth.shift$eff.sizes[truth.shift$eff.sizes>theta0],
                           seq(0,1,length.out = 26))


shift.strata.sim = run_simulation_strata(nreps=nreps,
                                         stratbrks=stratbrks.shift,
                                         params=params.shift,
                                         L=L,
                                         K=K,
                                         alpha.vec = alpha.vec)


#### Local Dispersion

#Generate strata based on a random draw from distribution and 
#fix the strata for nreps repetitions
set.seed(1)
X1 = rnorm.mixt(n=params.disp$n1,
                mus=params.disp$mu1,
                sigmas=params.disp$sigma1,
                props=params.disp$prop1)

X0 = rnorm.mixt(n=params.disp$n0,
                mus=params.disp$mu0,
                sigmas=params.disp$sigma0,
                props=params.disp$prop0)

truth.disp = get.truth(dat=c(X0,X1),m=2^K,cutoff=theta0,
                        func=eff.size,params=params.disp)

#specify 25 quantiles
stratbrks.disp = quantile(truth.disp$eff.sizes[truth.disp$eff.sizes>theta0],
                           seq(0,1,length.out = 26))

disp.strata.sim = run_simulation_strata(nreps=nreps,
                                      stratbrks=stratbrks.disp,
                                      params=params.disp,
                                      L=L,
                                      K=K,
                                      alpha.vec = alpha.vec)

#### Local Shift+Dispersion

#Generate strata based on a random draw from distribution and 
#fix the strata for nreps repetitions
set.seed(1)
X1 = rnorm.mixt(n=params.combo$n1,
                mus=params.combo$mu1,
                sigmas=params.combo$sigma1,
                props=params.combo$prop1)

X0 = rnorm.mixt(n=params.combo$n0,
                mus=params.combo$mu0,
                sigmas=params.combo$sigma0,
                props=params.combo$prop0)

truth.combo = get.truth(dat=c(X0,X1),m=2^K,cutoff=theta0,
                        func=eff.size,params=params.combo)

#specify 25 quantiles
stratbrks.combo = quantile(truth.combo$eff.sizes[truth.combo$eff.sizes>theta0],
                           seq(0,1,length.out = 26))

combo.strata.sim = run_simulation_strata(nreps=nreps,
                                       stratbrks=stratbrks.combo,
                                        params=params.combo,
                                        L=L,
                                        K=K,
                                        alpha.vec = alpha.vec)


############# Plot the results #############

#Figure settings
point.size=1.5
line.size=0.5
legend.text.size=16
figure.text.size=12

alphabysim = vector("list",length(alpha.vec))
layer.labels = sapply(seq(L),function(x) paste("L=",x,sep = ""))
for(a in seq_along(alpha.vec)){
  if(a==1){ #Plot first row of Figure 3, corresponding to first FDR value in alpha.vec
   
    #Local Shift
    numstrats.shift = length(stratbrks.shift)-1
    TEAM.shift.strat = vector("list",numstrats.shift)
    for(s in seq(numstrats.shift)){
      TEAM.shift.strat[[s]] = evaluate.sim(shift.strata.sim[,,s,a,])
    }
    TEAM.shift.strat.dat = do.call("rbind",lapply(TEAM.shift.strat,function(x) as.data.frame(x)))
    strat.labels = round(zoo::rollmean(stratbrks.shift,2),3)
    TEAM.shift.strat.dat$strat = rep(strat.labels,each=length(layer.labels))
    TEAM.shift.strat.dat$layer = rep(layer.labels,length(strat.labels))
    
    p1 = ggplot(TEAM.shift.strat.dat, aes(x=strat, y=Sensitivity, group=layer,colour = layer,shape=layer))
    tmp=p1 + geom_point(size=point.size) + 
      geom_line(size=line.size) + 
      theme_bw() +
      theme(text = element_text(size=figure.text.size),
            panel.grid=element_blank(),
            legend.position = "bottom",
            plot.title = element_text(vjust=1.5,hjust = 0.5,size=12,face="plain"),
            legend.title=element_text(size=12),
            legend.text = element_text(size = 12)) +
      guides(colour=guide_legend(nrow=1,byrow=TRUE)) +
      ggtitle("Local Shift\n") +
      xlab(expression(theta[i])) + 
      ylab("Avg Sensitivity") +
      scale_y_continuous(breaks=seq(0, 1, 0.2), limits=c(0, 1)) +
      scale_colour_discrete(name="Max Layer")+
      scale_shape_discrete(name="Max Layer")
    
    Shift.strat.Sensitivity = tmp + theme(legend.position = "none")
    
    layer.legend = get_legend(tmp + theme(plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),legend.direction = "horizontal",legend.justification="center"))
    
    #Local Dispersion
    numstrats.disp = length(stratbrks.disp)-1
    TEAM.disp.strat = vector("list",numstrats.disp)
    for(s in seq(numstrats.disp)){
      TEAM.disp.strat[[s]] = evaluate.sim(disp.strata.sim[,,s,a,])
    }
    TEAM.disp.strat.dat = do.call("rbind",lapply(TEAM.disp.strat,function(x) as.data.frame(x)))
    strat.labels = round(zoo::rollmean(stratbrks.disp,2),3)
    TEAM.disp.strat.dat$strat = rep(strat.labels,each=length(layer.labels))
    TEAM.disp.strat.dat$layer = rep(layer.labels,length(strat.labels))
    
    p2 = ggplot(TEAM.disp.strat.dat, aes(x=strat, y=Sensitivity, group=layer,colour = layer,shape=layer))
    Disp.strat.Sensitivity=p2 + geom_point(size=point.size) + 
      geom_line(size=line.size) + 
      theme_bw() +
      theme(text = element_text(size=figure.text.size),
            axis.title.y = element_blank(),
            panel.grid=element_blank(),
            legend.position = "none",
            plot.title = element_text(vjust=1.5,hjust = 0.5,size=12,face="plain"),    
            legend.title=element_text(size=12),
            legend.text = element_text(size = 12)) +
      guides(colour=guide_legend(nrow=1,byrow=TRUE)) +
      ggtitle("Local Dispersion\n") +
      xlab(expression(theta[i])) +
      scale_y_continuous(breaks=seq(0, 1, 0.2), limits=c(0, 1))
    
    #Local Shift+Dispersion
    numstrats.combo = length(stratbrks.combo)-1
    TEAM.combo.strat = vector("list",numstrats.combo)
    for(s in seq(numstrats.combo)){
      TEAM.combo.strat[[s]] = evaluate.sim(combo.strata.sim[,,s,a,])
    }
    TEAM.combo.strat.dat = do.call("rbind",lapply(TEAM.combo.strat,function(x) as.data.frame(x)))
    strat.labels = round(zoo::rollmean(stratbrks.combo,2),3)
    TEAM.combo.strat.dat$strat = rep(strat.labels,each=length(layer.labels))
    TEAM.combo.strat.dat$layer = rep(layer.labels,length(strat.labels))
    
    p3 = ggplot(TEAM.combo.strat.dat, aes(x=strat, y=Sensitivity, group=layer,colour = layer,shape=layer))
    tmp2=p3 + geom_point(size=point.size) + 
      geom_line(size=line.size) + 
      theme_bw() +
      theme(text = element_text(size=figure.text.size),
            panel.grid=element_blank(),
            axis.title.y = element_blank(),
            legend.position = "none",
            plot.title = element_text(vjust=1.5,hjust = 0.5,size=12,face="plain"),
            legend.title=element_text(size=12),
            legend.text = element_text(size = 12)) +
      guides(colour=guide_legend(nrow=1,byrow=TRUE)) +
      ggtitle("Local Shift + Dispersion\n") +
      xlab(expression(theta[i])) +
      scale_y_continuous(breaks=seq(0, 1, 0.2), limits=c(0, 1))
    
      Combo.strat.Sensitivity = grid.arrange(tmp2, right = textGrob(bquote(alpha==.(alpha.vec[a])), rot = -90, vjust = 1))
  } else{ #Plot subsequent rows of Figure 3, corresponding to subsequent FDR values in alpha.vec
    
    #Local Shift
    TEAM.shift.strat = vector("list",numstrats.shift)
    for(s in seq(numstrats.shift)){
      TEAM.shift.strat[[s]] = evaluate.sim(shift.strata.sim[,,s,a,])
    }
    TEAM.shift.strat.dat = do.call("rbind",lapply(TEAM.shift.strat,function(x) as.data.frame(x)))
    strat.labels = round(zoo::rollmean(stratbrks.shift,2),3)
    TEAM.shift.strat.dat$strat = rep(strat.labels,each=length(layer.labels))
    TEAM.shift.strat.dat$layer = rep(layer.labels,length(strat.labels))
    
    p1 = ggplot(TEAM.shift.strat.dat, aes(x=strat, y=Sensitivity, group=layer,colour = layer,shape=layer))
    Shift.strat.Sensitivity=p1 + geom_point(size=point.size) + 
      geom_line(size=line.size) + 
      theme_bw() +
      theme(text = element_text(size=figure.text.size),
            panel.grid=element_blank(),
            legend.position = "none") +
      guides(colour=guide_legend(nrow=1,byrow=TRUE)) +
      xlab(expression(theta[i])) + 
      ylab("Avg Sensitivity") +
      scale_y_continuous(breaks=seq(0, 1, 0.2), limits=c(0, 1))
    
    #Local Dispersion
    TEAM.disp.strat = vector("list",numstrats.disp)
    for(s in seq(numstrats.disp)){
      TEAM.disp.strat[[s]] = evaluate.sim(disp.strata.sim[,,s,a,])
    }
    TEAM.disp.strat.dat = do.call("rbind",lapply(TEAM.disp.strat,function(x) as.data.frame(x)))
    strat.labels = round(zoo::rollmean(stratbrks.disp,2),3)
    TEAM.disp.strat.dat$strat = rep(strat.labels,each=length(layer.labels))
    TEAM.disp.strat.dat$layer = rep(layer.labels,length(strat.labels))
    
    p2 = ggplot(TEAM.disp.strat.dat, aes(x=strat, y=Sensitivity, group=layer,colour = layer,shape=layer))
    Disp.strat.Sensitivity=p2 + geom_point(size=point.size) + 
      geom_line(size=line.size) + 
      theme_bw() +
      theme(text = element_text(size=figure.text.size),
            axis.title.y = element_blank(),
            panel.grid=element_blank(),
            legend.position = "none") +
      guides(colour=guide_legend(nrow=1,byrow=TRUE)) +
      xlab(expression(theta[i])) +
      scale_y_continuous(breaks=seq(0, 1, 0.2), limits=c(0, 1))
    
    #Local Shift+Dispersion
    TEAM.combo.strat = vector("list",numstrats.combo)
    for(s in seq(numstrats.combo)){
      TEAM.combo.strat[[s]] = evaluate.sim(combo.strata.sim[,,s,a,])
    }
    TEAM.combo.strat.dat = do.call("rbind",lapply(TEAM.combo.strat,function(x) as.data.frame(x)))
    strat.labels = round(zoo::rollmean(stratbrks.combo,2),3)
    TEAM.combo.strat.dat$strat = rep(strat.labels,each=length(layer.labels))
    TEAM.combo.strat.dat$layer = rep(layer.labels,length(strat.labels))
    
    p3 = ggplot(TEAM.combo.strat.dat, aes(x=strat, y=Sensitivity, group=layer,colour = layer,shape=layer))
    tmp2=p3 + geom_point(size=point.size) + 
      geom_line(size=line.size) + 
      theme_bw() +
      theme(text = element_text(size=figure.text.size),
            panel.grid=element_blank(),
            axis.title.y = element_blank(),
            legend.position = "none") +
      guides(colour=guide_legend(nrow=1,byrow=TRUE)) +
      xlab(expression(theta[i])) +
      scale_y_continuous(breaks=seq(0, 1, 0.2), limits=c(0, 1))
    Combo.strat.Sensitivity = grid.arrange(tmp2, right = textGrob(bquote(alpha==.(alpha.vec[a])), rot = -90, vjust = 1))
  }
  
  alphabysim[[a]] = plot_grid(Shift.strat.Sensitivity,
                              Disp.strat.Sensitivity,
                              Combo.strat.Sensitivity,
                              ncol=3,nrow=1)
  
}

strat.byalpha.legend = plot_grid(alphabysim[[1]],
                                  alphabysim[[2]],
                                  alphabysim[[3]],
                                  alphabysim[[4]],
                                  layer.legend,ncol=1,
                                  rel_heights = c(1.25,1,1,1,0.3))

ggsave("~/Figure3.pdf",
       strat.byalpha.legend,width=8,height=10,units="in")





