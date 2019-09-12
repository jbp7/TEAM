library(grid)
library(gridExtra)
library(dplyr)
library(cowplot)
library(MRS)
library(ggplotify)
library(ggplot2)
library(UpSetR)
source("~/TEAM.R")

####################  READ IN DATA #######################

#read in data using wildcards (glob2rx)

num_subj = 11

filelist = list()
for(i in 1:num_subj){
  filelist[[i]] <- list.files(".", pattern = glob2rx(paste0("*",ID,"*",".txt")))
} 

eqapol_datalist = list()
for(j in 1:num_subj){
  eqapol_datalist[[j]] = lapply(filelist[[j]], FUN=read.table, 
                                header=FALSE,skip=1,sep=",")
}

#assign channel names
channels <- c("FSC-A",
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
  names(eqapol_datalist[[k]])<- c("CMV","Costim")
  eqapol_datalist[[k]] <- lapply(eqapol_datalist[[k]],setNames,channels)
}

################# RUN TEAM and MRS ################
#Run on last four channels (i.e. functional channels)
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
    #and each bin has at least 100 events 
    K = 11 #floor(log2((n1+n0)/100))
    
    ####Center data
    if(center){
      #center by median
      median.comb = median(c(x1,x0))
      data.cent = c(x1,x0)-median.comb
      x1 = head(data.cent,n1)
      x0 = tail(data.cent,n0)
    }
    
    theta0=n1/(n1+n0)
    #TEAM (two-sample)
    TEAMres = TEAM(x1,x0,theta0=theta0,K=K,alpha=alpha,L=L)

    #Points falling into the bins
    cell_list_TEAM[[j]] = which(head(TEAMres$dat$index,n1)%in%TEAMres$S.list[[L]])
    
    #MRS
    X = c(x1,x0)
    G = rep(c(2,1),times=c(n1,n0))
    data = list(X,G)
    names(data) = c("X","G")
    
    ans = mrs(data$X,data$G,K=K)
    
    t = FDRThreshold(ans=ans,fdr = 0.05,by = 0.1)
    
    selected_regions = which(ans$RepresentativeTree$EffectSizes[, 1] > t)
    
    #Determine cells (events) in the discovered bins
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

###### DENSITY PLOTS ########
Dens_plots_by_subj = vector("list",num_subj)
for(i in 1:num_subj){
  plot.channel.1D = vector("list",length(func_channel_names))
  for (j in seq_along(func_channel_names)){
    plot.channel.1D[[j]] = plot1Ddens(TEAM_disc_by_subj[[i]]$res[[j]],L,channels.names[j])
  }
  #Save grobs
  Dens_plots_by_subj[[i]] = plot_grid(plotlist=plot.channel.1D,ncol=2)
  ggsave(paste0("ID",i,"_DensPlot.pdf"),
         upsetTEAManno,width=8,height=5,units="in")
}

###### UPSET PLOTS ######

#TEAM (one-sided) and MRS (two-sided)
channels.names = c("TNF-alpha","IL-2","IFN-gamma","CD107")

for(i in setdiff(seq_along(subj),c(6,11))){
  listInput_TEAM <- TEAM_disc_by_subj[[i]]$points#TEAM_disc_by_subj[[i]]$points#
  names(listInput_TEAM) = channels.names
  
  # listInput_MRS <- MRS_disc_by_subj[[i]]$points
  # names(listInput_MRS) = channels.names

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
    
    gr <- rectGrob(width=0.2,height=0.28, name="gr",gp=gpar(alpha=0.2,fill="red"))
    gr2 <- editGrob(gr, vp=viewport(x=0.86, y=0.21), name="gr2")
    
    upsetTEAManno <- gTree(children=gList(vpTEAM,gr2))
    
    # upset(fromList(listInput_MRS),order.by="degree",
    #       sets=channels.names, keep.order=TRUE,
    #       sets.x.label="Cells discovered/marker",
    #       mainbar.y.label = "Marker intersections (cells)",
    #       mb.ratio = c(0.65, 0.35),
    #       text.scale = c(1,0.8,1,0.8,1,1),
    #       decreasing=FALSE,
    #       empty.intersections = TRUE,
    #       scale.sets = "identity",
    #       scale.intersections="identity")
    # 
    # grid.edit('arrange',name='arrange2')
    # vpMRS = grid.grab()
    # 
    # gr <- rectGrob(width=0.2,height=0.28, name="gr",gp=gpar(alpha=0.2,fill="red"))
    # gr2 <- editGrob(gr, vp=viewport(x=0.86, y=0.21), name="gr2")
    # 
    # upsetMRSanno <- gTree(children=gList(vpMRS,gr2))
    # 
    # upsetplots <- grid.arrange(
    #   grobs=list(
    #     upsetTEAManno,
    #     upsetMRSanno
    #   ),
    #   nrow=1
    # )

  ggsave(paste0("ID",i,"_UpSetPlotsTEAM.pdf"),
         upsetTEAManno,width=8,height=5,units="in")
}

####### PLOTTING FUNCTIONS ###########

get.interval = function(x,df){
  tmp1 = unique(strsplit(gsub( "[][(]" , "", df$quant[which(df$index%in%x)]) , ","))
  tmp2 = lapply(tmp1,as.numeric)
  range(unlist(tmp2))
}

plot1Ddens <- function(res,L,channel){
  require(gridExtra)
  require(grid)
  
  ### Some pre-processing
  df = res$dat
  n1 = sum(df$lab)
  n0 = nrow(df)-n1
  x1 = head(df$X,n1)
  x0 = tail(df$X,n0)
  
  #Identify LDR intervals
  indx = res$S.list[[L]] 
  
  if(length(indx)>0){
    #Get contiguous bin indices
    indx.list = split(sort(indx), cumsum(c(1, diff(sort(indx)) != 1)))
    
    LDR.intervals=lapply(indx.list,get.interval,df=df) 
    start.intervals = unlist(lapply(LDR.intervals,`[`,1))
    end.intervals = unlist(lapply(LDR.intervals,`[`,2))
    
    df.intervals = data.frame(x1=start.intervals,y1=rep(-0.2,length(start.intervals)),
                              x2=end.intervals,y2=rep(-0.2,length(end.intervals)))
    
    ### Main plot
    p1 = ggplot(df, aes(x=X)) + 
      geom_line(aes(group=lab, colour=factor(lab),linetype = factor(3-lab)),stat="density",show.legend = FALSE) +
      scale_color_manual(values=c("black","red")) +
      ggtitle(channel) +
      theme(plot.title = element_text(hjust = 0.5),
            axis.ticks=element_blank(), panel.background=element_blank(),
            axis.text.x=element_blank(), axis.text.y=element_blank(),
            axis.title.x=element_blank(), axis.title.y=element_blank(),
            panel.border = element_blank(), plot.background = element_blank(),
            line = element_blank()) +
      geom_segment(data=df.intervals,aes(x = x1, y = y1, 
                                         xend = x2, yend = y2),
                   size = 1.03,lineend="butt")
    p1
  }
  #If no discoveries just plot densities
  else{
    p1 = ggplot(df, aes(x=X)) + 
      geom_line(aes(group=lab, colour=factor(lab),linetype = factor(3-lab)),stat="density",show.legend = FALSE) +
      scale_color_manual(values=c("black","red")) +
      ggtitle(channel) +
      theme(plot.title = element_text(hjust = 0.5),
            axis.ticks=element_blank(), panel.background=element_blank(),
            axis.text.x=element_blank(), axis.text.y=element_blank(),
            axis.title.x=element_blank(), axis.title.y=element_blank(),
            panel.border = element_blank(), plot.background = element_blank(),
            line = element_blank())
    p1
  }
}

