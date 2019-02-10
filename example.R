rm(list=ls())

n=200
K=772 
J=4
# simulate X, E, interactions
set.seed(K)
X<-matrix(rnorm(n*K),ncol=K)
set.seed(J)
E<-matrix(rnorm(n*J),ncol=J)
inter<-lapply(1:J, function(x){
  E[,x]*X
})
# simulate Y 
set.seed(J)
beta_E_sim<-runif(J)
set.seed(K)
beta_X_sim<-c(runif(5),rep(0,K-5))
beta_inter_sim<-do.call(rbind,lapply(1:J, function(x){
  set.seed(K+x)
 c(runif(2),rep(0,K-2))
}))

Y<-E%*%beta_E_sim + X%*%beta_X_sim +
  Reduce("+",lapply(1:J,function(x){
  inter[[x]]%*%(beta_inter_sim[x,]*beta_X_sim)
}))

##### marginal analysis #####
p_raw<-matrix(0,ncol=K,nrow=J+1)

for(k in 1:K){
  fit<-lm(Y~E+X[,k]+E*X[,k])
  temp<-summary(fit)
  p_raw[,k]<-temp$coefficients[6:10,4]
}
p_adjust<-apply(pval,2,function(x){
  p.adjust(x,method="fdr")
})

##### joint analysis #####
source("~/Path/joint_analysis_functions.R")
inter_s<-lapply(inter,function(x){scale(x)})
Y_s<-scale(Y)
fit_<-joint_model(E,X,inter_s,Y_s,15,1.5)
print(results(fit_,beta_X_sim,beta_inter_sim,Y_s))
