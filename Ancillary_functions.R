logistic1<-function(x){
  return(1/(1+exp(-x)))
}

FOC<-function(x){return(1/(exp(x)+2+exp(-x)))}

is.zero<-function(x){return(ifelse(abs(x)>0.00000001,1,0))}

struct<-function(d,thr){
  for (i in 1:dim(d)[1]){
    for(j in 1:dim(d)[2]){
      d[i,j]<-ifelse(abs(d[i,j])>thr,1,0)
    }
  }
  
  d<-is.zero(d+t(d))
  
  for (i in 1:dim(d)[1]){
    for(j in 1:dim(d)[2]){
      d[i,j]<-ifelse((i<j) & (d[i,j]>0.001),1,0)
    }
  }
  
  return(d)
}

rt<-function(true_mat,est_mat){
  TPR<-sum(true_mat==1 & true_mat==est_mat)/sum(true_mat==1)
  FPR<-sum(est_mat==1 & true_mat!=est_mat)/(dim(est_mat)[1]*(dim(est_mat)[1]-1)/2-sum(true_mat==1))
  return(c(TPR,FPR))
}

est_ad<-function(est_mat,way='union'){
  nd<-dim(est_mat)[1]
  mat_sign<-is.zero(est_mat)
  if (way=='union'){
    mat_sign<-is.zero(mat_sign+t(mat_sign))
  }else{mat_sign<-mat_sign*t(mat_sign)}
  for (i in 1:nd){
    for (j in 1:nd){
      if (i>=j){mat_sign[i,j]=0
      }
    }
  }
  return(mat_sign)
}

insert<-function(x_new,x,loc){
  if (loc==1){
    y<-append(x_new,x)
  }else if (loc==length(x)+1){
    y<-append(x,x_new)
  }else{y<-c(x[1:loc-1],x_new,x[loc:length(x)])}
  return(y)
}

nb_mat<-function(p){
  a<-matrix(0,nrow=p,ncol=p)
  for (i in 1:(p-1)){
   
    a[i,i+1]<-1
    
  }
  return(a)
}

tri_mat<-function(p,ac){
  mat<-diag(p)
  for (i in 1:p){
    for (j in 1:p){
      if (abs(i-j)==1){
        mat[i,j]=ac
      }
    }
  }
  return(mat)
}

r_auto_mvnorm<-function(p=5,acf=0.1,n=10,mu=0,sigma_=1){
  Sigma=diag(p)
  for (i in 1:p){
    for (j in 1:p){
      Sigma[i,j]<-acf^(abs(i-j))
    }
  }
  return(mvrnorm(n,rnorm(p,mu),sigma_*Sigma))
}

acf_mat<-function(p,ac){
  Sigma=diag(p)
  for (i in 1:p){
    for (j in 1:p){
      Sigma[i,j]<-ac^(abs(i-j))
    }
  }
  return(Sigma)
}


tri_strape_mat<-function(p,len=3,ac=0.5){
  Sigma=diag(p)
  for (i in 1:p){
    for (j in 1:p){
      if (abs(i-j)<=len & abs(i-j)!=0){Sigma[i,j]<-ac}
      
    }
  }
  return(Sigma)
}

block_construct<-function(block_num=5,p_list=rep(4,5),ac_list=rep(0.1,5),form=acf_mat){
  p_total<-sum(p_list)#total num of nodes
  
  beta<-diag(sum(p_list))
  
  beta[1:p_list[1],1:p_list[1]]<-form(p=p_list[1],ac=ac_list[1])
  
  for (i in 2:block_num){
    pre_<-sum(p_list[1:(i-1)])
    beta[(pre_+1):(pre_+p_list[i]),(pre_+1):(pre_+p_list[i])]<-form(p=p_list[i],ac=ac_list[i])
  }
  
  return(beta)
}