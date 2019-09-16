#####################################################################################
#
# Authors: J. Pura, C. Chan, PhD, J. Xie, PhD
#
# PURPOSE:
# 1) Generates Figures 4 and 5 in the main text of
#    Pura, Chan, Xie (2019). Multiple Testing Embedded in an Aggregation Tree to Identify 
#    where Two Distributions Differ 
#
# 2) Generates the values mentioned in paragraphs 3 and 4 of Section 6 in the main text.
#    These are obtained from the resulting UpSet plot analysis
#
# 3) Generates values in Table S1 of the Supplementary Material
# 
# 4) Generates Figures S1-S9 of the Supplementary Material.
#
# DETAILS:
# This code loads the EQAPOL data, runs TEAM and MRS, and plots the kernel density 
# estimates and discovered differential regions for the CMV-stimulated cells and 
# negative co-stimulated (Costim) cells across the four functional markers. Each 
# individual contains two files: IDXX_CMV_pp65.csv (case) and IDXX_Costim.csv (control), 
# where XX is the ID. UpSet plots of the discoveries for all 16 functionals marker 
# combinations are generated. 
#
# Figures 4 and 5 are generated for ID 11 using the following data files: 
# ID11_CMV_pp65.csv and ID11_Costim.csv. Note that TEAM identifies differential regions 
# where the case density exceeds the control density, while MRS identifies differential 
# regions where either density exceeds the other. The data is centered, prior to running
# TEAM or MRS, by subtracting off the median of the pooled data.
# 
#
########################################################################################

library(grid)
library(gridExtra)
library(dplyr)
library(cowplot)
library(MRS)
library(ggplotify)
library(ggplot2)
library(UpSetR)
library(TEAM)
source("~/case_study_functions.R")

####################  READ IN DATA #######################

#Set working directory
setwd("~/Dropbox/EQAPOL_normal/gated_export_xform_csv/EQAPOL_csv/")

#read in data using wildcards (glob2rx)
num_subj = 11

filelist = list()
for(i in 1:num_subj){
  filelist[[i]] = list.files(".", pattern = glob2rx(paste0("ID",i,"_","*",".csv")))
} 

eqapol_datalist = list()
for(j in 1:num_subj){
  eqapol_datalist[[j]] = lapply(filelist[[j]], FUN=read.table, 
                                header=FALSE,skip=1,sep=",")
}

#assign channel names
channels = c("FSC-A",
              "FSC-H",
              "FSC-W",
              "SSC-A",
              "SSC-H",
              "SSC-W",
              "CD57 FITC",
              "CD4 PerCP Cy55",
              "CD14 CD19 vAmine",
              "CD3 AmCyan",
              "CD27 APC",
              "TNFa A700",
              "CD8 APC Cy7",
              "IL 2 PE",
              "CD45RO ECD",
              "CD107 PE Cy5",
              "IFNg PE Cy7",
              "Time")


#Assign group and channel names
for (k in 1:num_subj){
  names(eqapol_datalist[[k]])= c("CMV","Costim")
  eqapol_datalist[[k]] = lapply(eqapol_datalist[[k]],setNames,channels)
}

################# Run TEAM and MRS ################
#For each subject:
#1) Pool case and control data and subtract off the pooled median from each value 
#2) Run TEAM and MRS on last four channels (i.e. functional channels)

L = 3
alpha = 0.05
center = TRUE
func_channel_names = c("TNFa A700","IL 2 PE","IFNg PE Cy7","CD107 PE Cy5")

TEAM_disc_by_subj = vector("list",num_subj)
MRS_disc_by_subj = vector("list",num_subj)

for (i in seq_along(subj)){
  #print(subj[i])
  cell_list_TEAM = vector("list",length(func_channel_names))
  cell_list_MRS = vector("list",length(func_channel_names))
  TEAMres = vector("list",length(func_channel_names))
  MRSres = vector("list",length(func_channel_names))
  
  for (j in seq_along(func_channel_names)){
    #print(func_channel_names[j])
    x1 = eqapol_datalist[[i]]$CMV[[func_channel_names[j]]]
    x0 = eqapol_datalist[[i]]$Costim[[func_channel_names[j]]]
    
    n1 = length(x1)
    n0 = length(x0)
    
    #Set number of bins to be a power of 2
    #and each bin has at least 100 pooled observations: #floor(log2((n1+n0)/100))
    #K=11 seems to work well for most of the data
    K = 11 
    
    ####Center data
    if(center){
      #center by median
      median.comb = median(c(x1,x0))
      data.cent = c(x1,x0)-median.comb
      x1 = head(data.cent,n1)
      x0 = tail(data.cent,n0)
    }
    
    theta0=n1/(n1+n0)
    
    #TEAM
    TEAMres[[j]] = TEAM(x1,x0,theta0=theta0,K=K,alpha=alpha,L=L)

    #CASE cell ID's discovered in each functional channel by TEAM 
    cell_list_TEAM[[j]] = which(head(TEAMres$dat$index,n1)%in%TEAMres$S.list[[L]])
    
    #MRS
    X = c(x1,x0)
    G = rep(c(2,1),times=c(n1,n0))
    data = list(X,G)
    names(data) = c("X","G")
    
    ans = mrs(data$X,data$G,K=K)
    
    t = FDRThreshold(ans=ans,fdr = 0.05,by = 0.1)
    
    selected_regions = which(ans$RepresentativeTree$EffectSizes[, 1] > t)
    
    #CASE cell ID's discovered in each functional channel by MRS  
    cell_list_MRS[[j]] = unique(unlist(ans$RepresentativeTree$DataPoints[selected_regions]))
    MRSres[[j]] = list("ans"=ans,
                       "regions"=selected_regions)
  }
  TEAM_disc_by_subj[[i]] = list("points"=cell_list_TEAM,"res"=TEAMres)
  MRS_disc_by_subj[[i]] = list("points"=cell_list_MRS,"res"=MRSres)
  
  names(TEAM_disc_by_subj[[i]]$points) = func_channel_names
  names(TEAM_disc_by_subj[[i]]$res) = func_channel_names
  names(MRS_disc_by_subj[[i]]$points) = func_channel_names
  names(MRS_disc_by_subj[[i]]$res) = func_channel_names
}

######### Plot Results #########
channels.names = c("TNF-alpha","IL-2","IFN-gamma","CD107")

###### Density Plots ########
Dens_plots_by_subj = vector("list",num_subj)

for(i in 1:num_subj){
  plot.channel.1D = vector("list",length(func_channel_names))
  for (j in seq_along(func_channel_names)){
    plot.channel.1D[[j]] = plot1Ddens(TEAM_disc_by_subj[[i]]$res[[j]],
                                      L,channels.names[j])
  }
  #Save grobs
  Dens_plots_by_subj[[i]] = plot_grid(plotlist=plot.channel.1D,ncol=2)
  ggsave(paste0("ID",i,"_DensPlot.pdf"),
         Dens_plots_by_subj[[i]],width=8,height=5,units="in")
}

###### Generate UpSet Plots from TEAM Results ######

for(i in setdiff(seq_along(subj),6)){ #ID 6 did not have any discoveries
  
  listInput_TEAM = TEAM_disc_by_subj[[i]]$points
  names(listInput_TEAM) = channels.names

  upset(fromList(listInput_TEAM),order.by="degree",
        sets=channels.names, keep.order=TRUE,
        sets.x.label="Cells discovered/marker",
        mainbar.y.label = "Marker intersections (cells)",
        mb.ratio = c(0.65, 0.35), 
        text.scale = c(1,0.8,1,0.8,1,1),
        decreasing=FALSE, 
        empty.intersections = TRUE,
        scale.sets = "identity",
        scale.intersections="identity")
    
    grid.edit('arrange',name='arrange2')
    vpTEAM = grid.grab()
    
    gr = rectGrob(width=0.2,height=0.28, name="gr",gp=gpar(alpha=0.2,fill="red"))
    gr2 = editGrob(gr, vp=viewport(x=0.86, y=0.21), name="gr2")
    
    upsetTEAManno = gTree(children=gList(vpTEAM,gr2))
  

  ggsave(paste0("ID",i,"_UpSetPlotsTEAM.pdf"),
         upsetTEAManno,width=8,height=5,units="in")
}

###### Generate UpSet Plot from MRS Result for ID=11 as in Figure 5 ######

  listInput_MRS = MRS_disc_by_subj[[11]]$points
  names(listInput_MRS) = channels.names
  
  upset(fromList(listInput_MRS),order.by="degree",
        sets=channels.names, keep.order=TRUE,
        sets.x.label="Cells discovered/marker",
        mainbar.y.label = "Marker intersections (cells)",
        mb.ratio = c(0.65, 0.35), 
        text.scale = c(1,0.8,1,0.8,1,1),
        decreasing=FALSE, 
        empty.intersections = TRUE,
        scale.sets = "identity",
        scale.intersections="identity")
  
  grid.edit('arrange',name='arrange2')
  vpMRS = grid.grab()
  
  gr = rectGrob(width=0.2,height=0.28, name="gr",gp=gpar(alpha=0.2,fill="red"))
  gr2 = editGrob(gr, vp=viewport(x=0.86, y=0.21), name="gr2")
  
  upsetMRSanno = gTree(children=gList(vpMRS,gr2))
  
  ggsave(paste0("ID11_UpSetPlotsMRS.pdf"),
         upsetMRSanno,width=8,height=5,units="in")

  
  
