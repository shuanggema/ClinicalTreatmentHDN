source("C:/Users/hm336/Box Sync/Work/Network/Aim2/Code/Ancillary functions.R")

#alpha_x,alpha_y,alpha_z,beta_x,beta_y,beta_z are p*p matrix, sigma,tao_y and tao_z are p*1 vector
#p is the number of the nodes
#link is a function
NW<-function(alpha_x,alpha_y,alpha_z,beta_x,beta_y,beta_z,tao_y,tao_z,p,link,n,thr,sigma){
  x<-rbinom(p,1,0.5)
  y<-rnorm(p)
  z<-rpois(p,1)
  
  i1=thr+n
  
  x_d<-matrix(1,nrow=i1,ncol=p)
  y_d<-matrix(1,nrow=i1,ncol=p)
  z_d<-matrix(1,nrow=i1,ncol=p)
  
  link<-link
  for (i in 1:i1){
    
   for (j in 1:p){
     s<-beta_x[j,-j]%*%x[-j]+beta_y[j,-j]%*%y[-j]+beta_z[j,-j]%*%z[-j]
     l<-alpha_x[j]+s
     prob<-link(l)
     x_d[i,j]<-rbinom(1,1,prob)
     y_d[i,j]<-ifelse(x_d[i,j]==0,0,rnorm(1,alpha_y[j]+tao_y[j]*s,sigma[j]))
     z_d[i,j]<-ifelse(x_d[i,j]==0,0,rpois(1,exp(alpha_z[j]+tao_z[j]*s)))
     x[j]<-x_d[i,j]
     y[j]<-y_d[i,j]
     z[j]<-z_d[i,j]
   }
    
  }
  return(list(x_d[(thr+1):i1,],y_d[(thr+1):i1,],z_d[(thr+1):i1,]))
}

p<-100
n<-100000

alpha_x<-rnorm(p,-3)#Change the mean can modify the proportion of zeros
alpha_y<-rnorm(p)
alpha_z<-rnorm(p)

# beta_x<-acf_mat(p,0.5)
# beta_y<-diag(p)
# beta_z<-diag(p)

####################### block setting
block_num<-5 #number of blocks
p_list<-rep(20,5)#number of nodes in each block
ac1_list<-rep(0.5,5)#ac parameter in each block
ac2_list<-rep(0.1,5)
p<-sum(p_list)#total num of nodes

beta_x<-block_construct(block_num=block_num,p_list=p_list,ac_list=ac1_list,form=tri_strape_mat)
beta_y<-block_construct(block_num=block_num,p_list=p_list,ac_list=ac2_list,form=tri_strape_mat)
beta_z<-block_construct(block_num=block_num,p_list=p_list,ac_list=ac2_list,form=tri_strape_mat)
#######################

tao_y<-rep(0.5,p)
tao_z<-rep(0.05,p)
sigma<-rep(1,p)

data=NW(alpha_x,alpha_y,alpha_z,beta_x,beta_y,beta_z,tao_y,tao_z,p,link=logistic1,n,thr=200000,sigma)
X<-data[[1]]
Y<-data[[2]]
Z<-data[[3]]

###########################################################################

zero.p<-NULL
for(i in 1:p){
  zero.p[i]<-1-sum(X[,i])/n
}
zero.p
mean(zero.p)

library(Matrix)
sum(Y)/nnzero(Y)
sum(Z)/nnzero(Z)

