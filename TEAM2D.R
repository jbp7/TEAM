############# LIBRARIES ################
require(plyr)
require(dplyr)

############# MAIN FUNCTION ################
TEAM2D = function(df0,df1,m,alpha,L){
  #df0 = data.frame for controls containing two columns: X and Y values 
  #df1 = data.frame for cases containing two columns: X and Y values 
  #m = number of breaks in EACH dimension
  #alpha = FDR level
  #L = number of layers
  
  N0 = nrow(df0)
  N1 = nrow(df1)
  
  dat = bind_rows(df0,df1)
  dat$lab = rep(c(0,1),times=c(N0,N1))
  
  #number of bins per dim
  #m = length(breaks$x)-1 #2^(K/2)
  
  #number of pooled observations per bin
  n = round((N0+N1)/(m^2))
  
  breaks <- create.xy.breaks(dat=dat,m=m)

  #Set breaks for each dimension
  brk.x <- breaks$x
  brk.y <- breaks$y
  
  init.l1.mat.xs <- matrix(NA,nrow=m,ncol=m) 
  init.l1.mat.n <- matrix(NA,nrow=m,ncol=m)
  for(i in seq_len(m)){
    init.l1.mat.n[,i] = breaks$y[[i]]$n
    init.l1.mat.xs[,i] = breaks$y[[i]]$x
  }

  #Global/Initial variables
  c.hats = vector("integer",L) #set of critical values
  S = NULL #set of rejections regions
  #m.excl = NULL #set of excluded regions
  S.list=vector("list",L)
  #Each element of the list is the index of bins in a row
  init.l1.indx.list = split(matrix(seq(m^2),ncol=m),seq(m))
  #init.l1.indx.vec = unlist(init.l1.indx.list)
  
  theta0 = (N1+0.5)/(N1+N0)
  
  #Loop through layers
  for (l in seq(L)){#2^(seq(L)-1)){
    #print(paste("l = ",l))
    if (l==1){
      
      m.1 = m^2
      # #print(m.1)
      x.1 = c(init.l1.mat.xs)
      
      c.hat.l = est.c.hat(l=1,n=n,
                          theta0=theta0,x.l=x.1,
                          c.hats=NULL,alpha=alpha,m.l=m.1)    
      
      #identify indices of rejected regions - discoveries
      S = c(which(x.1 > c.hat.l),S)#c(init.l1.indx.vec[which(x.1 > c.hat.l)],S) 
      
      S.list[[l]] = unique(unlist(c(S.list,list(S))))
      ### #print(c("Rejected at l= 1 :",S))
      ### #print(c("Number Rejected at l=1 :",length(S)))
    }
    else{# l > 1
      
      m.l = length(curr.layer.bin.indx.list)
      
      #print(c("m.l : ",m.l))
      
      #Change indexing from by column to by row
      curr.l1.indx.list.row = lapply(curr.l1.indx.list,ind.1D.each.row,m)
      
      #Counts in each bin
      curr.layer.xs.list = mapply(function(x,y) x[y],split(init.l1.mat.xs,seq(m)),
                                   curr.l1.indx.list.row, SIMPLIFY = FALSE)
      x.l = unlist(lapply(create.layer.bins(curr.layer.xs.list,l-1),sum))
      
      # #print(summary(x.l))
      #Total observations in each bin
      curr.layer.n.list = mapply(function(x,y) x[y],split(init.l1.mat.n,seq(m)),
                                 curr.l1.indx.list.row, SIMPLIFY = FALSE)

      c.hat.l = est.c.hat(l=l,n=n,
                          theta0=theta0,x.l=x.l,
                          c.hats=c.hats,alpha=alpha,m.l=m.l)       
      
      #identify indices of rejected regions - need to map to indices of l = 1
      #index of rejected counts corresponding to layer l>1
      rej.x = which(x.l > c.hat.l) 
      
      if (length(rej.x)==0){
        #print(c("Rejected at l=",l,":","NULL"))
        S = S
      } else{
        #indices corresponding to layer 1
        #ind.rej = unique(c(unlist(indx.list[rej.x],rej.xtra)))
        ind.rej = unique(unlist(curr.layer.bin.indx.list[rej.x]))
        #print(c("Rejected at l=",l,":",ind.rej))
        #print(paste("Number Rejected at l=",l,":",length(ind.rej)))
        S = c(ind.rej,S)
      }
      S.list[[l]] = unique(unlist(c(S.list,list(S))))
    }
    
    #append c hat vector
    #c.hats = c(c.hats,c.hat.l)
    c.hats[l] = c.hat.l
    
    #print(paste("c.hat: ",c.hat.l))
    
    #update current regions
    curr.l1.indx.list = lapply(init.l1.indx.list,function(x) setdiff(x,S))
    #bin indices for the layer l+1
    curr.layer.bin.indx.list = create.layer.bins(curr.l1.indx.list,l)
    ## #print(head(curr.layer.bin.indx.list))
  }
  return(list("n"=n, "m"=m,"S.list"=S.list,
              "c.hats"=c.hats,"data.orig"=dat,
              "breaks"=breaks,
              "mat.xs"=init.l1.mat.xs))
}

############### AUXILIARY FUNCTIONS ################

ind.1D.each.row <- function(ind.1D,factor){
  if(all(ind.1D%%factor==0,na.rm=TRUE)){
    return(ind.1D/factor)
  }
  else{
    return(floor(ind.1D/factor)+1)
  }
}

create.xy.breaks <- function(dat,m){
  require(dplyr)
  brk.x <- quantile(dat$X,seq(0,1,length.out = m+1))
  brk.y <- vector("list",m)
  
  for(i in seq_len(m)){
    x.subset = which(between(dat$X,brk.x[i],brk.x[i+1]))
    brk.y[[i]]$breaks = quantile(dat$Y[x.subset],seq(0,1,length.out = m+1))
    brk.y[[i]]$binnedpts = cut(dat$Y[x.subset],brk.y[[i]]$breaks,include.lowest = TRUE)
    brk.y[[i]]$x = aggregate(lab[x.subset]~brk.y[[i]]$binnedpts,data=dat,FUN=sum)$lab
    brk.y[[i]]$n = table(brk.y[[i]]$binnedpts)
  }
  return(list(x=brk.x,y=brk.y))
}

create.layer.bins <- function(indx.list,l){
  #Pad the end with NA for grouping in higher layers
  #indx.list.padded = lapply(indx.list,function(x) c(x,rep(NA,m-length(x))))
  indx.mat = as.matrix(ldply(indx.list,rbind)[-1])#do.call(rbind,indx.list.padded)
  
  #Group layer 1 bins across columns then across rows
  dim1 = floor(l/2)#ceiling(l/2)-1
  dim2 = ceiling(l/2)#floor(l/2)
  # #print(c("square size: ",2^dim1," x ",2^dim2))
  bin.indx.list = matsplitter(indx.mat,2^dim1,2^dim2)
  
  #Remove bins with NA
  #vector of booleans indicating if layer l bin has NA
  bin.indx.NA = unlist(lapply(bin.indx.list,function(x) all(!is.na(x))))
  
  return(bin.indx.list[bin.indx.NA])
}

create.layer.bins2 <- function(indx.list,l){
  #Pad the end with NA for grouping in higher layers
  #indx.list.padded = lapply(indx.list,function(x) c(x,rep(NA,m-length(x))))
  indx.mat = as.matrix(ldply(indx.list,rbind)[-1])#do.call(rbind,indx.list.padded)
  
  #Group layer 1 bins across columns then across rows
    dim1 = floor(l/2)#ceiling(l/2)-1
    dim2 = ceiling(l/2)#floor(l/2)
  # #print(c("square size: ",2^dim1," x ",2^dim2))
  bin.indx.list = matsplitter(indx.mat,2^dim1,2^dim2)
  
  #Remove bins with NA
  #vector of booleans indicating if layer l bin has NA
  bin.indx.NA = unlist(lapply(bin.indx.list,function(x) all(!is.na(x))))
  
  return(bin.indx.list[bin.indx.NA])
}


matsplitter<-function(M, r, c) {
  rg <- (row(M)-1)%/%r+1 # %/% = divide and round up
  cg <- (col(M)-1)%/%c+1
  rci <- (rg-1)*max(cg) + cg
  N <- prod(dim(M))/(r*c)
  lapply(1:N, function(x) M[rci==x])
}

est.c.hat <- function(l,n,theta0,x.l,c.hats,alpha,m.l){
  
  #print(paste("a.l:",floor(n*2^(l-1)*theta0+sqrt(2*theta0*(1-theta0)*n*2^(l-1)*log(m.l)))))
  #print(paste("2c-1:",2*c.hats[l-1]-1))
  
  # a.l =  max(2*c.hats[l-1]-1,
  #            floor(n*2^(l-1)*theta0+sqrt(2*theta0*(1-theta0)*n*2^(l-1)*log(m.l))))
  
  a.l =  floor(n*2^(l-1)*theta0+sqrt(2*theta0*(1-theta0)*n*2^(l-1)*log(m.l)))
  
  
  #print(paste("max.x:",max(x.l)))
  
  filter1 = which(x.l < 2^(l-1)*n*theta0)
  filter2 = which(x.l > a.l)
  c.vec = sort(unique(x.l[-c(filter1,filter2)]),decreasing = TRUE)
  
  x.indx = 0
  emp.fdp = 0
  
  
  if(length(c.vec)>0){
    while(all(emp.fdp <= alpha, x.indx < length(c.vec))){ 
      x.indx = x.indx + 1
      emp.fdp = est.FDP.hat.l(min.x=c.vec[x.indx],
                              max.x=max(n,2*c.hats[l-1]),
                              c.prev=ifelse(length(c.hats[l-1])>0,c.hats[l-1],NULL),
                              n.l=n*2^(l-1),
                              x.l=x.l,
                              theta0=theta0,
                              l=l)
      #print(paste("x.indx",x.indx,"c[x.indx]",c.vec[x.indx]))
      ##print(emp.fdp)
    }
    ##print(paste("x.indx-1",x.indx-1))
    c.hat = ifelse(is.na(c.vec[x.indx-1]) || length(c.vec[x.indx-1])==0,
                   a.l,
                   c.vec[x.indx-1])
    
  }else{
    c.hat = a.l
  }
  return(c.hat)
}



est.FDP.hat.l <- function(min.x,max.x,c.prev,n.l,x.l,theta0,l){
  
  #print(paste("n.l:",n.l))
  
  num = 0 #for ### ###printing only
  denom = 0
  ##print(paste("min.x = ", min.x,"max.x = ", max.x))
  if(is.na(min.x)) return(1) #need to add condition to break out of function if last value of x.vec is reached
  
  if(min.x + 1 > max.x){ 
    FDP.est = 0
  }
  else{
    tmp = seq(from=min.x+1,to=max.x) #corresponds to x.l > 2^(l-1)*n*theta0, # > xvec[j]
    ## ###print(c("tmp:",tmp))
    if(l==1){
      num = sum(dbinom(tmp,n.l,theta0))
      denom = max(sum(x.l>min.x),1)
      FDP.est = length(x.l)*num/denom 
    }
    else{  
      f1 <- function(z){
        #sums over the dbinoms for vector of a_{1,r}'s...
        #need row products
        return(ifelse(!is.matrix(z),
                      sum(prod(dbinom(x=z,size=n.l/2,prob=theta0)),na.rm = TRUE),
                      sum(apply(dbinom(x=z,size=n.l/2,prob=theta0),1,prod),na.rm = TRUE)))
      }
      f2 <- function(w){
        #sum over the values of c_k's
        #print(head(valid.counts(w,c.prev=c.prev)))
        val = ifelse(is.list(valid.counts(w,c.prev=c.prev)),
                     sum(sapply(valid.counts(w,c.prev=c.prev),f1),na.rm=TRUE),
                     ifelse(is.null(valid.counts(w,c.prev=c.prev)),0,
                            f1(valid.counts(w,c.prev=c.prev))))
        return(val)
      }
      
      #numerator
      num = f2(tmp)/pbinom(q=c.prev,size=n.l/2,prob=theta0)^2
      
      #denominator
      denom = max(sum(x.l>min.x),1)
      FDP.est = length(x.l)*num/denom
    }
  }
  #print(paste("num:",num,"denom:",denom,"m.l:",length(x.l),"FDP:",FDP.est))
  return(FDP.est)
}

expand.mat = function(mat, vec) {
  #### #print(paste(nrow(mat),length(vec)))
  #### #print(paste("nrow=",nrow(mat)))
  out = matrix(0, nrow = as.numeric(nrow(mat)) * as.numeric(length(vec)),
               ncol = as.numeric(ncol(mat) + 1)) #deal with integer overflow
  for (i in 1:ncol(mat)) out[, i] = mat[, i]
  out[, ncol(mat) + 1] = rep(vec, each = nrow(mat))
  return(out)
}

valid.counts = function(x,c.prev){
  #Outputs a 2-column matrix of valid counts 
  vec = seq(from = 0, to = c.prev)
  mat = matrix(vec, ncol = 1)
  
  mat  = expand.mat(mat, vec)
  mat = mat[rowSums(mat) %in% x, ]
  
  if(!is.matrix(mat)){
    return(mat)
  }
  else if(is.matrix(mat) & nrow(mat)<1){
    return(NULL)
  }
  else{
    #mat = mat[rowSums(mat) %in% x, ]
    idx = split(seq(nrow(mat)), rowSums(mat))
    return(lapply(idx, function(i, x) x[i, ], mat))
  }
}
