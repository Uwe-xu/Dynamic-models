# This version includs Hafner and Rombout's work

# manual generating 
rm(list=ls(all=TRUE))
#library(matlib)
library(fGarch)
library(xtable)
# require("fBasics")
# require("timeDate")
# require("timeSeries")
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
#library(rugarch)

n_vector=c(250,500,1000)
#n_vector=c(250)
t_df_vector=c(3,4,6,9,20)
#t_df_vector=c(6)
m=3
n_start=1000
NumIt = 1000
sigma=1

k=1
theta=c(omega=1,alpha1=0.6,beta1=0.3)
source("Functions_Garch(1,1).R")

# setting initial values

#F_i is for generating, and F_i2 is for estimation
F_i=F_i2=F1

estimates_Mfit<- matrix(0, nrow= NumIt, ncol= m)
estimates_OPIV<- matrix(0, nrow= NumIt, ncol= m)
starter_sse = proc.time()
time_index=starter_sse

  
  for (d_s in 1: length(t_df_vector)){
  t_df=t_df_vector[d_s]
  result_MISE=matrix(0, nrow= 1+2*length(n_vector), ncol= 2*m+4)
  result_MISE[1,1:3]=theta
  for (n_s in 1: length(n_vector)){
    n=n_vector[n_s]
  
  print(paste("n=",n,"and df=",t_df))

ones=matrix(1,n-1,1, byrow = FALSE)
for (j in 1: NumIt){
  #process generating
  source("generating_Garch(1,1).R")
# MLE 
  init=theta
  m_mit  <- try(optim(init, M_fit,method = "L-BFGS-B",lower=c(0.01,0,0.01),upper=c(5,0.99,0.99)),silent=TRUE)
  if(inherits(m_mit, "try-error")) {
    m_mit  <- optim(init, M_fit)
  } 
  estimates_Mfit[j,]=m_mit$par
  init=m_mit$par
  
# Optimal IV
  #calculating k3,k4
  source("k3_k4_Garch(1,1).R")
  oiv  <- try(optim(init, IVLS,method = "L-BFGS-B",lower=c(0.01,0,0.01),upper=c(5,0.99,0.99)),silent=TRUE)
  if(inherits(oiv, "try-error")) {
    oiv  <- optim(init, IVLS)
  } 
  estimates_OPIV[j,]=oiv$par

# counter 
  if (j/100==floor(j/100))
  {print(c(j,(proc.time() -time_index)[3]))
    time_index=proc.time() }
}
filename <- paste("n=", n,"_df=",t_df,".csv", sep = "")

ones=matrix(1,NumIt,1, byrow = FALSE)

true=theta
MLE_mean=colMeans(estimates_Mfit)
GEE_mean=colMeans(estimates_OPIV)
MLE_mse=colMeans((estimates_Mfit-ones%*%theta)^2)*n
GEE_mse=colMeans((estimates_OPIV-ones%*%theta)^2)*n
start=2+(n_s-1)*2
end=3+(n_s-1)*2
result_MISE[start:end,1:3]= rbind(MLE_mean,GEE_mean)
 
result_MISE[start:end,4]=9999
  rbind(MLE_mean,GEE_mean)
result_MISE[start:end,5:7]=rbind(MLE_mse,GEE_mse)
result_MISE[start:end,8]=8888
result_MISE[start:end,9]=n
result_MISE[start:end,10]=t_df

time_sse = (proc.time() -starter_sse)[3]
time_sse
  }
  colnames(result_MISE) <- c("w","a","b","space","mse_w","mse_a","mse_b","space","n","df")
  names_vector<- c("MLE","GEE")
  rownames(result_MISE)<- c("true",rep(names_vector, times = length(n_vector)))
  result.mat = round(result_MISE,digits=4)
  #save(result.mat, file="sim_gls_DGP3-MSE.Rdata")
  result.mat
  latex_table <- xtable(result.mat,digits = 4)
  filename=paste("result_of_dof=", t_df, ".csv", sep = "")
  
  write.csv(result.mat, file = filename)
  filename=paste("result_of_dof=", t_df, ".tex", sep = "")
  capture.output(print(latex_table, type = "latex", include.rownames = TRUE), file = filename)
  }
time_sse = (proc.time() -starter_sse)[3]
time_sse

