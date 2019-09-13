###################################################################################
#
# Authors: J. Pura, C. Chan, PhD, J. Xie, PhD
#
# PURPOSE:
# Generates Figure 2 in the main text of
# Pura, Chan, Xie (2019). Multiple Hypothesis Testing on an Aggregation Tree Method 
#
# DETAILS:
# This code generates generates the data, runs TEAM, and plots the avg FDP, avg FN,
# and avg total discoveries across layers and FDR for each of the three local settings: 
# local shift, dispersion, and shift+dispersion. 
#
###################################################################################

library(ggplot2)
library(cowplot)
library(ggpubr)

source("~/simulation_functions.R")

# Number of repetitions
nreps=1000

# Parameters for TEAM
L = 5
K = 14
alpha.vec = c(0.01,seq(2.5,20,by=2.5))

# Dataset 
n = 180
n1 = n0 = 0.5*n*2^K

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

shift.layer.sim = run_simulation_layer(nreps=nreps,
                                        params=params.shift,
                                        L=L,
                                        K=K,
                                        alpha.vec=alpha.vec)

disp.layer.sim = run_simulation_layer(nreps=nreps,
                                       params=params.disp,
                                       L=L,
                                       K=K,
                                       alpha.vec=alpha.vec)

combo.layer.sim = run_simulation_layer(nreps=nreps,
                                        params=params.combo,
                                        L=L,
                                        K=K,
                                        alpha.vec=alpha.vec)


############# Plot the results #############

##### Create density plots ######

#Figure settings
point.size=3
line.size=1
legend.text.size=16
figure.text.size=12
inset.text.size=5

#Local Shift
x1.pdf = dnorm.mixt(seq(0,1,length.out = 10000),mus=params.shift$mu1,sigmas=params.shift$sigma1,props=params.shift$prop1)
x0.pdf = dnorm.mixt(seq(0,1,length.out = 10000),mus=params.shift$mu0,sigmas=params.shift$sigma0,props=params.shift$prop0)

df.pdf = data.frame(X1=seq(0,1,length.out=10000),Y1=x1.pdf,X0=seq(0,1,length.out=10000),Y0=x0.pdf)

main.plot.s = ggplot(df.pdf, aes(x=X1,y=Y1,color='red')) + 
  geom_line() +
  geom_line(data=df.pdf,aes(x=X0,y=Y0),linetype="dashed",color='black') +
  theme(text=element_text(size=figure.text.size),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position="none",
        plot.title = element_text(vjust=3,size=14,face="plain")) +
  ggtitle("Local Shift\n") +
  xlab("Location") +
  ylab("Density") +
  geom_rect(mapping=aes(xmin=0.8,xmax=.97,ymin=-0.1,ymax=1.75),
            alpha=0,linetype=2,color='black',size=0.2)

sub.plot.s = ggplot(subset(df.pdf,between(X1,0.8,0.97))) + 
  geom_line(aes(x=X1,y=Y1,color='red'))+
  geom_line(data=subset(df.pdf,between(X0,0.8,97)),aes(x=X0,y=Y0),linetype="dashed",color='black') +
  theme(text = element_text(size=inset.text.size),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title.x=element_blank(), 
        axis.title.y=element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position="none") 

main.sub.s = ggdraw() +
  draw_plot(main.plot.s, 0, 0, 1, 1) +
  draw_plot(sub.plot.s, x=0.4, y=0.45, width=0.55,height=0.45)

#Local Dispersion
x1.pdf = dnorm.mixt(seq(0,1,length.out = 10000),mus=params.disp$mu1,sigmas=params.disp$sigma1,props=params.disp$prop1)
x0.pdf = dnorm.mixt(seq(0,1,length.out = 10000),mus=params.disp$mu0,sigmas=params.disp$sigma0,props=params.disp$prop0)

df.pdf = data.frame(X1=seq(0,1,length.out=10000),Y1=x1.pdf,X0=seq(0,1,length.out=10000),Y0=x0.pdf)
main.plot.d = ggplot(df.pdf, aes(x=X1,y=Y1,color='red')) + 
  geom_line() +
  geom_line(data=df.pdf,aes(x=X0,y=Y0),linetype="dashed",color='black') +
  theme(text=element_text(size=figure.text.size),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position="none",
        plot.title = element_text(vjust=3,size=14,face="plain")) +
  ggtitle("Local Dispersion\n") +
  xlab("Location") + 
  geom_rect(mapping=aes(xmin=0.65,xmax=.95,ymin=-0.1,ymax=0.8),alpha=0,linetype=2,color='black',size=0.2)

sub.plot.d = ggplot(subset(df.pdf,between(X1,0.6,1))) + 
  geom_line(aes(x=X1,y=Y1,color='red'))+
  geom_line(data=subset(df.pdf,between(X0,0.6,1)),aes(x=X0,y=Y0),linetype="dashed",color='black') +
  theme(text = element_text(size=inset.text.size),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title.x=element_blank(), 
        axis.title.y=element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position="none") 

main.sub.d = ggdraw() +
  draw_plot(main.plot.d, 0, 0, 1, 1) +
  draw_plot(sub.plot.d, x=0.53, y=0.45, width=0.42,height=0.45)

#Local Shift+Disp
x1.pdf = dnorm.mixt(seq(0,1,length.out = 10000),mus=params.combo$mu1,sigmas=params.combo$sigma1,props=params.combo$prop1)
x0.pdf = dnorm.mixt(seq(0,1,length.out = 10000),mus=params.combo$mu0,sigmas=params.combo$sigma0,props=params.combo$prop0)

df.pdf = data.frame(X1=seq(0,1,length.out=10000),Y1=x1.pdf,X0=seq(0,1,length.out=10000),Y0=x0.pdf)

main.plot.sd = ggplot(df.pdf, aes(x=X1,y=Y1,color='red')) + 
  geom_line() +
  geom_line(data=df.pdf,aes(x=X0,y=Y0),linetype="dashed",color='black') +
  theme(text=element_text(size=figure.text.size),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position="none",
        plot.title = element_text(vjust=3,size=14,face="plain")) +
  ggtitle("Local Shift + Dispersion\n") +
  xlab("Location") +
  geom_rect(mapping=aes(xmin=0.65,xmax=.95,ymin=-0.1,ymax=0.5),alpha=0,linetype=2,color='black',size=0.2)

sub.plot.sd = ggplot(subset(df.pdf,between(X1,0.6,1))) + 
  geom_line(aes(x=X1,y=Y1,color='red'))+
  geom_line(data=subset(df.pdf,between(X0,0.6,1)),aes(x=X0,y=Y0),linetype="dashed",color='black') +
  theme(text = element_text(size=inset.text.size),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title.x=element_blank(), 
        axis.title.y=element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position="none") 

main.sub.sd = ggdraw() +
  draw_plot(main.plot.sd, 0, 0, 1, 1) +
  draw_plot(sub.plot.sd, x=0.519, y=0.45, width=0.43,height=0.45)

##### Create FDR, FN, and Total Discovery plots ######

### Local Shift

# Performance by layers
TEAM.shift.layer = vector("list",L)
for(l in seq(L)){
  TEAM.shift.layer[[l]] = evaluate.sim(shift.layer.sim[,l,,])
}
TEAM.shift.layer.dat = do.call("rbind",lapply(TEAM.shift.layer,function(x) as.data.frame(x)))

layer.labels = sapply(seq(L),function(x) paste("L=",x,sep = ""))
TEAM.shift.layer.dat$layer = as.factor(rep(layer.labels,each=length(alpha.vec)))
TEAM.shift.layer.dat$alpha.val = rep(alpha.vec,times=length(layer.labels))

#Avg FDP
p1 = ggplot(TEAM.shift.layer.dat, aes(x=alpha.val, y=FDP, group=layer,colour = layer,shape=layer))
Shift.FDP=p1 + geom_point(size=point.size) +
  geom_line(aes(x=alpha.val,y=alpha.val),colour='black',size=line.size) + 
  theme_bw() +
  theme(text = element_text(size=figure.text.size),
        legend.position = "none",
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  xlab(expression(alpha)) +
  ylab("Avg FDP")

#Avg FN
p2 = ggplot(TEAM.shift.layer.dat, aes(x=alpha.val, y=FN, group=layer,colour = layer,shape=layer))
Shift.FN = p2 + geom_point(size=point.size) +
  geom_line(size=line.size) + 
  theme_bw() +
  theme(text = element_text(size=figure.text.size),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "none") +
  xlab(expression(alpha)) +
  ylab("Avg FN")

#Avg total discoveries

p3 = ggplot(TEAM.shift.layer.dat, aes(x=alpha.val, y=num.rej, group=layer,colour = layer,shape=layer))
Shift.numrej = p3 + geom_point(size=point.size) +
  geom_line(size=line.size) + 
  theme_bw() +
  theme(text = element_text(size=figure.text.size),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "none") +
  xlab(expression(alpha)) +
  ylab("Avg Tot. Disc.")


### Local Dispersion
# Performance by layers
TEAM.disp.layer = vector("list",L)
for(l in seq(L)){
  TEAM.disp.layer[[l]] = evaluate.sim(disp.layer.sim[,l,,])
}
TEAM.disp.layer.dat = do.call("rbind",lapply(TEAM.disp.layer,function(x) as.data.frame(x)))

layer.labels = sapply(seq(L),function(x) paste("L=",x,sep = ""))
TEAM.disp.layer.dat$layer = as.factor(rep(layer.labels,each=length(alpha.vec)))
TEAM.disp.layer.dat$alpha.val = rep(alpha.vec,times=length(layer.labels))

#Avg FDP
p1 = ggplot(TEAM.disp.layer.dat, aes(x=alpha.val, y=FDP, group=layer,colour = layer,shape=layer))
Disp.FDP=p1 + geom_point(size=point.size) +
  geom_line(aes(x=alpha.val,y=alpha.val),colour='black',size=line.size) + 
  theme_bw() +
  theme(text = element_text(size=figure.text.size),
        legend.position = "none",
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title.y = element_blank()) +
  xlab(expression(alpha))

#Avg FN
p2 = ggplot(TEAM.disp.layer.dat, aes(x=alpha.val, y=FN, group=layer,colour = layer,shape=layer))
Disp.FN = p2 + geom_point(size=point.size) +
  geom_line(size=line.size) + 
  theme_bw() +
  theme(text = element_text(size=figure.text.size),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        axis.title.y=element_blank()) +
  xlab(expression(alpha))

#Avg Total Discoveries

p3 = ggplot(TEAM.disp.layer.dat, aes(x=alpha.val, y=num.rej, group=layer,colour = layer,shape=layer))
Disp.numrej = p3 + geom_point(size=point.size) +
  geom_line(size=line.size) + 
  theme_bw() +
  theme(text = element_text(size=figure.text.size),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        axis.title.y=element_blank()) +
  xlab(expression(alpha))


### Local Shift+Disp

# Performance by layers
TEAM.combo.layer = vector("list",L)
for(l in seq(L)){
  TEAM.combo.layer[[l]] = evaluate.sim(combo.layer.sim[,l,,])
}
TEAM.combo.layer.dat = do.call("rbind",lapply(TEAM.combo.layer,function(x) as.data.frame(x)))

layer.labels = sapply(seq(L),function(x) paste("L=",x,sep = ""))
TEAM.combo.layer.dat$layer = as.factor(rep(layer.labels,each=length(alpha.vec)))
TEAM.combo.layer.dat$alpha.val = rep(alpha.vec,times=length(layer.labels))

#Avg FDP 
p1 = ggplot(TEAM.combo.layer.dat, aes(x=alpha.val, y=FDP, group=layer,colour = layer,shape=layer))
Combo.FDP=p1 + geom_point(size=point.size) +
  geom_line(aes(x=alpha.val,y=alpha.val),colour='black',size=line.size) + 
  theme_bw() +
  theme(text = element_text(size=figure.text.size),
        legend.position = "none",
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title.y = element_blank()) +
  xlab(expression(alpha))

#Avg FN
p2 = ggplot(TEAM.combo.layer.dat, aes(x=alpha.val, y=FN, group=layer,colour = layer,shape=layer))
tmp = p2 + geom_point(size=point.size) +
  geom_line(size=line.size) + 
  theme_bw() +
  theme(text = element_text(size=figure.text.size),
        axis.title.y = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position="bottom",
        legend.title=element_text(size=legend.text.size),
        legend.text = element_text(size = legend.text.size)) +
  xlab(expression(alpha)) +
  scale_colour_discrete(name="Max Layer")+
  scale_shape_discrete(name="Max Layer")

Combo.FN = tmp + theme(legend.position = "none")
legend = get_legend(tmp + theme(legend.direction = "horizontal",legend.justification="center" ,legend.box.just = "bottom"))

#Avg Total Discoveries
p3 = ggplot(TEAM.combo.layer.dat, aes(x=alpha.val, y=num.rej, group=layer,colour = layer,shape=layer))
Combo.numrej = p3 + geom_point(size=point.size) +
  geom_line(size=line.size) + 
  theme_bw() +
  theme(text = element_text(size=figure.text.size),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        axis.title.y=element_blank()) +
  xlab(expression(alpha))

#Combine pdf, FDR, FN, and legend
sim.all = ggarrange(main.sub.s,main.sub.d,main.sub.sd,
                     Shift.FDP,Disp.FDP,Combo.FDP,
                     Shift.FN,Disp.FN,Combo.FN,
                     Shift.numrej,Disp.numrej,Combo.numrej,
                     ncol = 3, align = "v")


sim.all.legend = plot_grid(sim.all[[1]],sim.all[[2]],
                            sim.all[[3]],sim.all[[4]],
                            legend,ncol=1,rel_heights = c(1,1,1,1,0.3))

ggsave("~/Figure2.pdf",sim.all.legend,width=12,height=12,units="in")












