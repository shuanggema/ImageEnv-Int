### counting tp/fp
tpfp<-function(a, b){
  tp <- sum( a !=0 & b!=0)
  fp <- sum( a !=0 & b==0)
  return(list(tp=tp,fp=fp))
}

### estimated y hat
gf<-function(cbeta_env,cbeta_main,cbeta_inter,
             cenviro,cmain,cinter){
  main_ef<-cenviro%*%cbeta_env+ cmain%*%cbeta_main

  inter_eff<-lapply(1:J,function(i,a,b){a[[i]]%*%(b[i,]*cbeta_main)},a=cinter,b=cbeta_inter)
  return(main_ef+Reduce("+",inter_eff))
}

#### objective function 
qf<-function(cbeta_main,cbeta_inter, ## vector
             gf_value,y,
             lam1,lam2){
  mse_sq<-mean(sum((y-gf_value)^2))
  pen<-lam1*sum(abs(cbeta_main))+lam2*sum(abs(cbeta_inter))
  
  return(sum(pen)+mse_sq)
}

### soft_thresholding
sthres<-function(a,b){
  if(abs(a)<=b){
    return(0) 
  }else{
    return(sign(a)*(abs(a)-b))
  }
}

######## algorithm ########
joint_model<-function(env,cmain,cinter,y,lam1,lam2){
  #### starting values 
  beta_env_old<-rep(0,J)
  beta_main_old<-rep(0,K)
  beta_inter_old<-matrix(0,nrow=J,ncol=K)
  yfit<-rep(0,n)
  q<-c()
  loop=0
  d=1
  while(d>10^-4 && loop<201 && d<10^4){
    loop=loop+1
    q_old<-qf(beta_main_old,beta_inter_old,yfit,y,lam1,lam2)
    ## update alpha
    ytilda<-y-yfit+env%*%beta_env_old
    beta_env_old<-solve(t(env)%*%env)%*%t(env)%*%ytilda
    yfit<-y-ytilda+env%*%beta_env_old

    ## update gamma_k = eta k
    xtilda<-cmain+Reduce("+",lapply(1:J,function(i,a,b){a[[i]]%*%diag(b[i,])},a=cinter,b=beta_inter_old))
    beta_s<-beta_main_old
    for(i in 1:K){
      ytilda<-y-yfit+xtilda%*%beta_s
      noi<-ytilda-xtilda[,-i]%*%beta_s[-i]
      beta_s[i]<-sthres((t(xtilda[,i])%*%xtilda[,i])^(-1)*xtilda[,i]%*%noi,(t(xtilda[,i])%*%xtilda[,i])^(-1)*lam1)
      yfit<-y-ytilda+xtilda%*%beta_s
    }
  
    beta_main_old<-beta_s
    ## update tau = theta jk
    xtilda<-lapply(cinter,function(x){x%*%diag(beta_s)})
    tau<-beta_inter_old
    id_nonz<-which(beta_s!=0)
    print( length(id_nonz))
    for(l in 1:J){
      tau_nol<-tau[-l,]
      xtilda_nol<-xtilda[-l]
      for(j in id_nonz){
        ytilda<-y-yfit+Reduce("+",lapply(1:J,function(i,a,b){a[[i]]%*%b[i,]},a=xtilda,b=tau))
        noi<-ytilda-Reduce("+",lapply(1:(J-1),function(i,a,b){a[[i]]%*%b[i,]},a=xtilda_nol,b=tau_nol))

        tau[l,j]<-sthres((t(xtilda[[l]][,j])%*%xtilda[[l]][,j])^(-1)*xtilda[[l]][,j]%*%noi,(t(xtilda[[l]][,j])%*%xtilda[[l]][,j])^(-1)*lam2)
        yfit<-y-ytilda+Reduce("+",lapply(1:J,function(i,a,b){a[[i]]%*%b[i,]},a=xtilda,b=tau))
      }
    }
    beta_inter_old<-tau
    q_new<-qf(beta_main_old,beta_inter_old,yfit,y,lam1,lam2)
    d<-abs(q_new-q_old)/abs(q_old)
    print(d)
    q[loop]<-q_new
  }
  if(d>10^4){
    return(print("not converge"))
  }else{
    ## return the estimated para and seq of objective 
    return(list(beta_main=beta_main_old,beta_inter=beta_inter_old,
                beta_env=beta_env_old,yhat=yfit,obj=q))
  }
}
#### tuning 
results<-function(fit,cbeta_main,cbeta_inter,y){
  
  main_tpfp<-do.call(c,tpfp(fit_$beta_main,cbeta_main))
  inter_tpfp<--do.call(c,tpfp(fit_$beta_inter,cbeta_main))
  
  mse<-mean(sum((fit$yhat-y)^2))
  df_new<-sum(main_tpfp)+sum(inter_tpfp)+J
  
  ebic<-n*log(mse)+df_new*log(n)+2*1*log(choose(K*J+J,df_new)) 
  
  
  results_v<-c(main_tpfp,inter_tpfp,
               df_new-J,
               mse,
               ebic)
names(results_v)<-c("m:tp","m:fp","i:tp","i:fp",
                         "df","mse","ebic")
  return(results_v)
}

