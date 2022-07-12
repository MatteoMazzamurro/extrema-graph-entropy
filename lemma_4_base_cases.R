#Verification of base cases for Lemma 4
#in Resolution of Yan’s Conjecture on Entropy of Graphs 
#(https://arxiv.org/abs/2205.03357)

#Load necessary library
library(rootSolve)

#Define the functions
#LL(x,c)=(x+1)*(x+c)*log(x+c)-(c)*log(c)-x*(x-1+c)*log(x-1+c) and
#RL(x,c)=((x+1)*x/2+c)*log((x+1)*x/2+c)-((x-1)*x/2+c)*log((x-1)*x/2+c)+x*((1+c)*log(1+c)-(c)*log(c)).
#We want to find for which values of x with x>=2, we have LL(x,c)<RL(x,c) when c=1,2,3

#Define RL-LL
RL_LL<-function(x,c){
  RL_LL_x_c<-((x+1)*x/2+c)*log((x+1)*x/2+c)-((x-1)*x/2+c)*log((x-1)*x/2+c)+x*((1+c)*log(1+c)-(c)*log(c))-((x+1)*(x+c)*log(x+c)-(c)*log(c)-x*(x-1+c)*log(x-1+c))
  return(RL_LL_x_c)
}

#Define specific functions for c=1,2,3 respectively 
RL_LL_1<-function(x) RL_LL(x,1)
RL_LL_2<-function(x) RL_LL(x,2)
RL_LL_3<-function(x) RL_LL(x,3)

#Find roots for x>=2
uniroot.all(RL_LL_1,c(2,100))
uniroot.all(RL_LL_2,c(2,100))
uniroot.all(RL_LL_3,c(2,100))
