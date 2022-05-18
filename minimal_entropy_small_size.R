#load necessary packages
library(partitions)
library(igraph)

#define function h_c
h_c <- function(x,c){
  #avoid issues with log(0)
  if(c==0){
    x[x==0] <- 1
  }
  ans <- sum((x+c)*log(x+c))
  return(ans)
}

#allow tolerance in the comparison to account for machine precision
tol<-1e-6

#initialize empty list of degree sequences
h_c_max <- list()

for(m in 3:10){
  #find partitions of 2m 
  parts_m <- as.matrix(parts(2*m))
  #select partitions that are valid degree sequences
  is_graphical_m <- apply(parts_m,MARGIN=2,is_graphical)
  graphical_parts_m <- as.matrix(parts_m[,is_graphical_m])
  #find which degree sequence(s) maximize(s) h_c
  h_c_m <- apply(graphical_parts_m,MARGIN=2,h_c,c=1)
  h_c_m_max <- graphical_parts_m[,which((max(h_c_m)-h_c_m)<tol)]
  #add the degree sequence(s) maximizing h_c_m to the list
  h_c_max <- append(h_c_max,list(h_c_m_max))
}
