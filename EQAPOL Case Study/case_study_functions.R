###################################################################################
#
# Authors: J. Pura, C. Chan, PhD, J. Xie, PhD
#
#
# DETAILS:
# This code contains functions to plot the kernel density estimates of the CMV (case)
# and Costim (control) cells. Additionally, it contains a function to compute the
# threshold that controls the FDR defined in Soriano and Ma (2018).
# 
#
###################################################################################

####### Plotting Functions ###########

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


##### FDR threshold for MRS #####
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

