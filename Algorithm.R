set<-1
rep<-1

setwd("/gpfs/ysm/home/hm336/aim2/")

source("Ancillary functions.R")

load(paste("xyz_single_setting",set,"_rep",rep,"_data.RData", sep=""))
n<-dim(X)[1]
p<-dim(X)[2]

library(doMPI)

cl = startMPIcluster()           
registerDoMPI(cl)
clusterSize(cl)

X.sample<-X[sample,]
Y.sample<-Y[sample,]
Z.sample<-Z[sample,]

n.sample<-dim(X.sample)[1]

##################################################### initial value
q1<-Sys.time()

myFun1<-function(X.sample, Y.sample, Z.sample, n.sample, p){
  tmp.list1<-foreach(j=1:p)%dopar%{ 
    
    t<-0.1
    thr<-0.01
    lambda<-0.1
    
    X2.sample.j<-X.sample[,-j]*X.sample[,-j]
    Y2.sample.j<-Y.sample[,-j]*Y.sample[,-j]
    Z2.sample.j<-Z.sample[,-j]*Z.sample[,-j]
      
    alpha.x<-0
    alpha.y<-0
    alpha.z<-0
    beta.x<-rep(0,(p-1))
    beta.y<-rep(0,(p-1))
    beta.z<-rep(0,(p-1))
    tao.y<-0
    tao.z<-0
    sigma<-1
    
    k<-0
    
    while(1){
      k<-k+1
      
      para.original<-c(alpha.x,alpha.y,alpha.z,beta.x,beta.y,beta.z,tao.y,tao.z,sigma)
      
      l<-X.sample[,-j]%*%beta.x+Y.sample[,-j]%*%beta.y+Z.sample[,-j]%*%beta.z
      
      qj<-alpha.x+l
      muj<-alpha.y+tao.y*l
      lj<-alpha.z+tao.z*l
      
      d_ljtaozj<-l
      d_mujtaoyj<-l
      
      d_Lqj<--exp(qj)/(1+exp(qj))+X.sample[,j]
      d_Lqj[which(is.na(d_Lqj))]<--1+X.sample[which(is.na(d_Lqj)),j]
      d_Lmuj<-X.sample[,j]*(Y.sample[,j]-muj)/sigma^2
      d_Llj<-X.sample[,j]*(Z.sample[,j]-exp(lj))
      
      delta_alpha_xj<-sum(d_Lqj)/n.sample+2*lambda*alpha.x/n.sample
      delta_alpha_yj<-sum(d_Lmuj)/n.sample+2*lambda*alpha.x/n.sample
      delta_alpha_zj<-sum(d_Llj)/n.sample+2*lambda*alpha.x/n.sample
      
      delta_beta_xj<-t(X.sample[,-j])%*%(d_Lqj+tao.y*d_Lmuj+tao.z*d_Llj)/n.sample+2*lambda*beta.x/n.sample
      delta_beta_yj<-t(Y.sample[,-j])%*%(d_Lqj+tao.y*d_Lmuj+tao.z*d_Llj)/n.sample+2*lambda*beta.x/n.sample
      delta_beta_zj<-t(Z.sample[,-j])%*%(d_Lqj+tao.y*d_Lmuj+tao.z*d_Llj)/n.sample+2*lambda*beta.x/n.sample
      
      delta_tao_yj<-sum(d_Lmuj*d_mujtaoyj)/n.sample+2*lambda*tao.y/n.sample
      delta_tao_zj<-sum(d_Llj*d_ljtaozj)/n.sample+2*lambda*tao.z/n.sample
      
      delta_sigma<-sum(-X.sample[,j]/sigma+X.sample[,j]*(Y.sample[,j]-muj)^2/sigma^3)/n.sample+2*lambda*sigma/n.sample
      
      d2_L_qj<--exp(qj)/(1+exp(qj))^2
      d2_L_qj[which(is.na(d2_L_qj))]<-0
      d2_L_muj<--X.sample[,j]/sigma^2
      d2_L_lj<--X.sample[,j]*exp(lj)
      
      delta2_alpha_xj<-sum(d2_L_qj)/n.sample+2*lambda/n.sample
      delta2_alpha_yj<-sum(d2_L_muj)/n.sample+2*lambda/n.sample
      delta2_alpha_zj<-sum(d2_L_lj)/n.sample+2*lambda/n.sample
      
      delta2_beta_xj<-(t(X2.sample.j)%*%(d2_L_qj+tao.y^2*d2_L_muj+tao.z^2*d2_L_lj))/n.sample+2*lambda/n.sample
      delta2_beta_yj<-(t(Y2.sample.j)%*%(d2_L_qj+tao.y^2*d2_L_muj+tao.z^2*d2_L_lj))/n.sample+2*lambda/n.sample
      delta2_beta_zj<-(t(Z2.sample.j)%*%(d2_L_qj+tao.y^2*d2_L_muj+tao.z^2*d2_L_lj))/n.sample+2*lambda/n.sample
      
      delta2_tao_yj<-sum(d2_L_muj*d_mujtaoyj^2)/n.sample+2*lambda/n.sample
      delta2_tao_zj<-sum(d2_L_lj*d_ljtaozj^2)/n.sample+2*lambda/n.sample
      
      delta2_sigma<-sum(X.sample[,j]/sigma^2-3*X.sample[,j]*(Y.sample[,j]-muj)^2/sigma^4)/n.sample+2*lambda/n.sample
      
      alpha.x<-alpha.x-t*delta_alpha_xj/delta2_alpha_xj
      alpha.y<-alpha.y-t*delta_alpha_yj/delta2_alpha_yj
      alpha.z<-alpha.z-t*delta_alpha_zj/delta2_alpha_zj
      
      beta.x<-beta.x-t*delta_beta_xj/delta2_beta_xj
      beta.y<-beta.y-t*delta_beta_yj/delta2_beta_yj
      beta.z<-beta.z-t*delta_beta_zj/delta2_beta_zj
      
      tao.y<-tao.y-t*delta_tao_yj/delta2_tao_yj
      tao.z<-tao.z-t*delta_tao_zj/delta2_tao_zj
      
      sigma<-sigma-t*delta_sigma/delta2_sigma
      
      para.update<-c(alpha.x,alpha.y,alpha.z,beta.x,beta.y,beta.z,tao.y,tao.z,sigma)
      
      if(all(abs((para.update-para.original)/para.original)<thr)){
        break
      }else if(k>10000){
        break
      }
    }  
    
    rm(list=(c("l","qj","muj","lj",ls(pattern="^d_"),ls(pattern="^d2_"),ls(pattern="^delta_"),ls(pattern="^delta2_"),ls(pattern="^para."))))
    list(alpha.x,alpha.y,alpha.z,beta.x,beta.y,beta.z,tao.y,tao.z,sigma)
  }
  return(tmp.list1)
}

para.initial.list<-myFun2(X.sample, Y.sample, Z.sample, n.sample, p)

q2<-Sys.time()
q2-q1

save(X, Y, Z, sample, para.initial.list, file=paste("xyz_single_setting",set,"_rep",rep,"_initial.RData", sep=""))

##################################################
q1<-Sys.time()

lambda_seq<-c(seq(0.011,0.02,0.001))

myFun2<-function(X.sample, Y.sample, Z.sample, n.sample, p, lambda_seq, para.initial.list){
  tmp.list2<-foreach(lambda=lambda_seq) %:%
    foreach(j=1:p) %dopar%{
      
      .GlobalEnv$tmp.list<-para.initial.list
      
      thr<-0.0001
      t<-0.01
      
      alpha.x<-tmp.list[[j]][[1]]
      alpha.y<-tmp.list[[j]][[2]]
      alpha.z<-tmp.list[[j]][[3]]
      
      beta.x<-tmp.list[[j]][[4]]
      beta.y<-tmp.list[[j]][[5]]
      beta.z<-tmp.list[[j]][[6]]
      
      tao.y<-tmp.list[[j]][[7]]
      tao.z<-tmp.list[[j]][[8]]
      
      sigma<-tmp.list[[j]][[9]]
      
      k<-0
      
      while(1){    
        
        k<-k+1
        
        para.original<-c(alpha.x,alpha.y,alpha.z,beta.x,beta.y,beta.z,tao.y,tao.z,sigma)
        
        l<-X.sample[,-j]%*%beta.x+Y.sample[,-j]%*%beta.y+Z.sample[,-j]%*%beta.z
        
        qj<-alpha.x+l
        muj<-alpha.y+tao.y*l
        lj<-alpha.z+tao.z*l
        
        d_ljtaozj<-l      
        d_mujtaoyj<-l
        
        d_Lqj<--exp(qj)/(1+exp(qj))+X.sample[,j]
        d_Lqj[which(is.na(d_Lqj))]<--1+X.sample[which(is.na(d_Lqj)),j]
        d_Lmuj<-X.sample[,j]*(Y.sample[,j]-muj)/sigma^2
        d_Llj<-X.sample[,j]*(Z.sample[,j]-exp(lj))
        
        delta_alpha_xj<-sum(d_Lqj)/n.sample
        delta_alpha_yj<-sum(d_Lmuj)/n.sample
        delta_alpha_zj<-sum(d_Llj)/n.sample
        
        delta_beta_xj<-t(X.sample[,-j])%*%(d_Lqj+tao.y*d_Lmuj+tao.z*d_Llj)/n.sample
        delta_beta_yj<-t(Y.sample[,-j])%*%(d_Lqj+tao.y*d_Lmuj+tao.z*d_Llj)/n.sample
        delta_beta_zj<-t(Z.sample[,-j])%*%(d_Lqj+tao.y*d_Lmuj+tao.z*d_Llj)/n.sample
        
        delta_tao_yj<-sum(d_Lmuj*d_mujtaoyj)/n.sample
        delta_tao_zj<-sum(d_Llj*d_ljtaozj)/n.sample
        
        delta_sigma<-sum(-X.sample[,j]/sigma+X.sample[,j]*(Y.sample[,j]-muj)^2/sigma^3)/n.sample
        
        alpha.x<-alpha.x+t*delta_alpha_xj
        alpha.y<-alpha.y+t*delta_alpha_yj
        alpha.z<-alpha.z+t*delta_alpha_zj
        
        beta.x<-beta.x+t*delta_beta_xj
        beta.y<-beta.y+t*delta_beta_yj
        beta.z<-beta.z+t*delta_beta_zj
        
        tao.y<-tao.y+t*delta_tao_yj
        tao.z<-tao.z+t*delta_tao_zj
        
        sigma<-sigma+t*delta_sigma
        
        beta<-rbind(t(beta.x),t(beta.y),t(beta.z))
        beta<-apply(beta,2,function(m){L<-sqrt(sum(m^2))
        if (L>lambda*t){return((L-lambda*t)*m/L)
        }else(return(c(0,0,0)))})
        
        beta.x<-beta[1,]
        beta.y<-beta[2,]
        beta.z<-beta[3,]
        
        para.update<-c(alpha.x,alpha.y,alpha.z,beta.x,beta.y,beta.z,tao.y,tao.z,sigma)


        
        if(all(abs(para.update-para.original)<thr)){
          break
        }else if(k>30000){
          break
        }
      }
      
      rm(list=(c("l","qj","muj","lj",ls(pattern="^d_"),ls(pattern="^d2_"),ls(pattern="^delta_"),ls(pattern="^delta2_"),ls(pattern="^para."))))
      list(beta.x,beta.y,beta.z)
    }
  return(tmp.list3)
}

para.roc.list<-myFun(X.sample, Y.sample, Z.sample, n.sample, p, lambda_seq, para.initial.list)

q2<-Sys.time()
q2-q1

closeCluster(cl)

######################################

save(X, Y, Z, sample, para.true.list, para.initial.list, para.roc.list, file=paste("xyz_single_setting",set,"_rep",rep,".RData", sep=""))

mpi.quit()
