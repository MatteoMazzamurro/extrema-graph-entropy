#This R script supports the finding of Cambie, Dong, Mazzamurro's paper
#"Extremal values of degree-based entropies of bipartite graphs"
#############################
#Structure:
#0. load packages
#1. essential functions
#2. analysis
#2.1 B(n,m,y) yields minimal entropy for some y (small n)
#2.2 properties of minimal entropy B(n,m,y)
#############################
#0. packages
{
  library(partitions)
}
#1. essential functions
{
  #we define functions for the following tasks:
  #compute generalised Shannon entropy h_c
  #find all difference graphs
  #find the degree sequence of B(n,m,y)
  #find the degree sequence and properties of interest of B(n,m,y) for all valid values y
  
  #generalised Shannon entropy h_c
  h_c <- function(x,c){
    #avoid issues with log(0)
    if(c==0){
      x[x==0] <- 1
    }
    ans<-sum((x+c)*log(x+c))
    return(ans)
  }
  
  #By Lemma 2 in (http://arxiv.org/abs/2108.13884v1), a bipartite graph (U,V,E) yielding minimal entropy is a difference graph
  #This function finds all difference graphs given order n and size m
  #each difference graph corresponds to a partition and its conjugate
  #(see N.V.R. Mahadev and U.N. Peled, Threshold graphs and related topics, Elsevier Science B.V., Amsterdam, 1995)
  #output of the function is a matrix with all valid degree sequences, and for each sequence,
  #the cardinalities u and v of the two sets of nodes U and V, and the size m (for convenience).
  difference_graphs<-function(n,m){
    #initialise empty dataframe of degree sequencex
    difference_graphs_ds<-NULL
    #split nodes into two parts (avoiding empty set on either side)
    parts_n<-restrictedparts(n,2)[,-1]
    #loop on each partition
    for (i in 1:ncol(parts_n)) {
      #take u>=v
      u<-parts_n[1,i]
      v<-parts_n[2,i]
      #find all partitions of m into u parts (it's best to work with smaller value u)
      u_parts<-restrictedparts(m,u)
      #see whether all parts are <=v
      u_parts_valid_col<-(u_parts[1,]<=v)
      #avoid issues when there are no graphic partitions for these values of u and v
      if(any(u_parts_valid_col)){
        graphic_u_parts<-as.matrix(u_parts[,u_parts_valid_col])
        graphic_v_parts<-as.matrix(conjugate(graphic_u_parts))
        #first three entries are always u,v, and m
        deg_sequences<-cbind(u,v,m,t(graphic_u_parts),t(graphic_v_parts))
      }else{
        deg_sequences<-NULL
      }
      #fill dataframe with valid degree sequences
      difference_graphs_ds<-rbind(difference_graphs_ds,deg_sequences)
    }
    return(difference_graphs_ds)
  }
  
  #function to compute degree sequence of B(n,m,y), a promising candidate for the minimum entropy graph 
  #(see http://arxiv.org/abs/2108.13884v1)
  #In the article "#"Extremal values of degree-based entropies of bipartite graphs"
  #we prove that this is indeed the case for some choice(s) of 1<=y<=\sqrt(m)
  B_n_m_y<-function(n,m,y){
    q<-floor(m/y)
    r<-m-y*q
    #make sure y and q are suitable values
    #we need y+q<=n when r=0 and y+q+1<=n when r>0
    if(y+q+(r>0)<=n){
      deg_sec<-c(rep(y,q),r,rep(q+1,r),rep(q,y-r))
      #If you want the degree sequence to be sorted, uncomment the following line:
      #deg_sec<-sort(deg_sec,decreasing=TRUE)
    }else{
      deg_sec<-NULL 
    }
    return(deg_sec)
  }
  
  #For given n and m, the function B(n,m) computes B(n,m,y) for 1<=y<=\sqrt(m) 
  #This allows to see for which y the entropy is minimal. 
  #When there exists y such that m=qy, with q+y<=n, then B(n,m,y)=K_{q,y}U\notK_{n-q-y} yields minimal entropy 
  #(as proved in http://arxiv.org/abs/2108.13884v1) 
  #This means that for some values of n and m, several values of y (and q) may yield minimum entropy.
  #e.g.: for n=20, m=60, one can take (y,q)\in{(4,15),(5,12),(6,10)} and for all such values B(n,m,y) yields minimum entropy
  #When no such y exists, then determining which values of y and q yield minimal entropy is non-trivial.
  #For each choice of y, we record the value x, where if we let r=m-qy, then x=q when r=0, and x=q+1 when r>0.  
  #This allows us to study the divergence from the "ideal" situation m=qy using both m-qy and xy-m.
  B_n_m<-function(n,m){
    #dataframe of result
    m_y_q_x_h<-data.frame(
      n=numeric(),
      m=numeric(),
      y=numeric(),
      q=numeric(),
      x=numeric(),
      m_qy=numeric(),
      xy_m=numeric(),
      h=numeric()
    )
    i<-1
    for(y in 1:sqrt(m)){
      #find degree sequence
      q<-floor(m/y)
      B_n_m_y_ds<-B_n_m_y(n,m,y)
      B_n_m_y_ds<-B_n_m_y_ds[B_n_m_y_ds>0]
      r<-m-y*q
      #avoid complications when y is not a valid choice
      if(!is.null(B_n_m_y_ds[1])){
        #which position the largest degree occupies in the second part depends on whether m=qy (r==0) 
        #x<-B_n_m_y_ds[ifelse(r==0,q+1,q+2)]
        x<-ifelse(r==0,q,q+1)
        B_n_m_y_hc<-h_c(B_n_m_y_ds,0)
        m_y_q_x_h[i,]<-c(n,m,y,q,x,m-q*y,x*y-m,B_n_m_y_hc)
        i<-i+1
        #If you want a description of the outcome
        #print(paste("y=",y,"; q=",q,"; x=",x,"; m-qy=",m-q*y,"; xy-m=",x*y-m,"; Deg sequence: ",paste(B_n_m_y_ds,collapse=" "),"; h=",B_n_m_y_hc,sep=""))
      }
    }
    return(m_y_q_x_h)
  }
}
#2. analysis
{
  #2.1 verify that some B(n,m,y) graph achieves minimal entropy 
  {
    #min_bip_n_m finds the bipartite graph with minimum entropy by brute force among all (n,m)-difference graphs
    #Here m<=max_m=floor(n/2)*ceiling(n/2)
    min_bip_n_m<-function(n,m){
      #find all (n,m)-difference graphs
      dg<-difference_graphs(n,m)
      #vector of entropies
      dg_h<-numeric(nrow(dg))
      for(i in 1:nrow(dg)){
        dg_h[i]<-h_c(dg[i,-c(1,2,3)],0)
      }
      dg<-cbind(dg,dg_h)
      #deg sequences of minimal entropy bipartite graphs and h (last entry)
      min_bip_ent<-dg[dg_h==max(dg_h),]
      return(min_bip_ent)
    }
    
    #Given n, verify that B(n,m,y) really yields minimum entropy for all m and for some value of y
    #use a brute force approach of finding minimal entropy graph among all difference graphs
    #set n
    n<-15
    #max value of m depends on n
    max_m<-floor(n/2)*ceiling(n/2)
    
    #first compute all difference graphs by brute force 
    #for each value of m, store deg seq and entropy of the one with minimum entropy
    #(first 3 entries from function difference_graphs(n,m) are always cardinalities u, v and size m, hence the plus 4)
    min_bip_n<-matrix(nrow=(max_m-n+1),ncol=n+4)
    for (m in n:max_m) {
      #find all (n,m)-difference graphs
      dg<-difference_graphs(n,m)
      #vector of entropies
      dg_h<-numeric(nrow(dg))
      for(i in 1:nrow(dg)){
        dg_h[i]<-h_c(dg[i,-c(1,2,3)],0)
      }
      dg<-cbind(dg,dg_h)
      #deg sequence of a minimal entropy bipartite graph and h (last entry)
      min_bip_n[m-n+1,]<-as.vector(dg[which.max(dg_h),])
    }
    
    #compute minimal entropy B(n,m,y) graphs for each m 
    min_b_n<-numeric(max_m-n+1)
    for(m in n:max_m){
      m_y_q_x_h<-na.omit(B_n_m(n,m))
      min_b_n[m-n+1]<-max(m_y_q_x_h$h)
    }
    
    #check if the entropies are the same 
    is_b_n_m_min_logical<-all(min_bip_n[,n+4]==min_b_n)
  }
  
  #2.2 find minimal entropy graphs among B(n,m,y)'s and study their properties
  {
    #Find for which values of y B(n,m,y) yields minimum entropy
    #initialize empty dataframe of significant values.
    #the last three entry denote whether the graph achieves minimal entropy, qy-m, and xy-m respectively
    explore_B_n_m<-data.frame(
      n=numeric(),
      m=numeric(),
      y=numeric(),
      q=numeric(),
      x=numeric(),
      m_qy=numeric(),
      xy_m=numeric(),
      h=numeric(),
      is_h_max=logical(),
      is_m_qy_min=logical(),
      is_xy_m_min=logical()
    )
    
    #to ensure correct comparison, introduce a small tolerance to account for machine precision
    tol<-1e-5
    n_min<-2
    n_max<-50
    for (n in n_min:n_max){
      print(n)
      m_max<-floor(n/2)*ceiling(n/2)
      for(m in n:m_max){
        m_y_q_x_h<-B_n_m(n,m)
        #find max for h and min for m_qy and xy_m
        h_max<-max(m_y_q_x_h$h)
        m_qy_min<-min(m_y_q_x_h$m_qy)
        xy_m_min<-min(m_y_q_x_h$xy_m)
        #add to dataframe
        m_y_q_x_h$is_h_max<-(abs(m_y_q_x_h$h-h_max)<tol)
        m_y_q_x_h$is_m_qy_min<-(m_y_q_x_h$m_qy==m_qy_min)
        m_y_q_x_h$is_xy_m_min<-(m_y_q_x_h$xy_m==xy_m_min)
        explore_B_n_m<-rbind(explore_B_n_m,m_y_q_x_h)
      }
    }
    
    #Remarks:
    #When n<=50, at least one of xy-m or m-qy is minimal when h is maximal
    explore_B_n_m[explore_B_n_m$is_h_max==TRUE&explore_B_n_m$is_m_qy_min==FALSE&explore_B_n_m$is_xy_m_min==FALSE,]
    
    #This is not always the case. See for example
    n<-17726
    m<-318728
    m_y_q_x_h<-B_n_m(n,m)
    #find max for h and min for m_qb and xy_m
    h_max<-max(m_y_q_x_h$h)
    m_qy_min<-min(m_y_q_x_h$m_qy)
    xy_m_min<-min(m_y_q_x_h$xy_m)
    #add to dataframe
    m_y_q_x_h$is_h_max<-(abs(m_y_q_x_h$h-h_max)<tol)
    m_y_q_x_h$is_m_qy_min<-(m_y_q_x_h$m_qy==m_qy_min)
    m_y_q_x_h$is_xy_m_min<-(m_y_q_x_h$xy_m==xy_m_min)
    example_B<-rbind(example_B,m_y_q_x_h)
    example_B[example_B$is_h_max==TRUE,]
  }
}