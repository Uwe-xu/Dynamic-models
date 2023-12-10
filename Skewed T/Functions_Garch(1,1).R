#initial values
F1<- function(w,a,b){ 
  init=abs(w/(1-a-b))
  return(init)
}

F2<- function(w,a,b){ 
  init=abs(w)
  return(init)
}

# estiamtion functions
M_fit <- function(x_ini) {  
  w <- x_ini[1]
  a <- x_ini[2]
  b <- x_ini[3]
  
  sig2=e2=matrix(0,n+1,1, byrow = FALSE)
  sig2[1]=e2[1]=F_i2(w,a,b)
  e2[2:(n+1)]=Y^2
  
  for (f_i in 1:n){
    sig2[1+f_i]=abs(w+a*e2[f_i]+b*sig2[f_i]) 
  }
  
  e2_1=e2[1:n]
  e2_2=e2[2:(n+1)]
  sig2_1=sig2[1:n]
  sig2_2=sig2[2:(n+1)]
  
  E1=(e2_2-sig2_2)/(sig2_2)^2
  E2=E1*e2_1
  E3=E1*sig2_1
  
  sum(E1)^2+sum(E2)^2+sum(E3)^2
}
IVLS <- function(x_ini) {
  
  w <- x_ini[1]
  a <- x_ini[2]
  b <- x_ini[3]
  
  sig2=e2=e1=matrix(0,n+1,1, byrow = FALSE)
  sig2[1]=e2[1]=abs(w/(1-a-b))
  sig2[1]=e2[1]=F_i2(w,a,b)
  e1[1]=e2[1]^0.5
  e2[2:(n+1)]=Y^2
  for (f_i in 1:n){
    sig2[1+f_i]=abs(w+a*e2[f_i]+b*sig2[f_i])
  }
  e2_1=e2[1:n]
  e2_2=e2[2:(n+1)]
  sig2_1=sig2[1:n]
  sig2_2=sig2[2:(n+1)]
  
  E1=(e2_2-sig2_2)/(sig2_2)^2-k3/(sig2_2)^1.5*Y
  E2=E1*e2_1
  E3=E1*sig2_1
  
  sum(E1)^2+sum(E2)^2+sum(E3)^2
}

#Hafner and Rombouts 2013
  # derivative of normal kernal.
d_n<- function(x){ 
     K=(2*pi)^(-0.5)*exp(-x^2/2)*(-1)*x
    return(K)
}

# HR
SP_HR<- function(x_ini) {
  w <- x_ini[1]
  a <- x_ini[2]
  b <- x_ini[3]
  
  sig2=e2=e1=matrix(0,n+1,1, byrow = FALSE)
  sig2[1]=e2[1]=F_i2(w,a,b)
  e1[1]=e2[1]^0.5
  e1[2:(n+1)]=Y
  e2[2:(n+1)]=Y^2
  
  for (f_i in 1:n){
    sig2[1+f_i]=abs(w+a*e2[f_i]+b*sig2[f_i]) 
  }
  
  e1_2=e1[2:(n+1)]
  e2_1=e2[1:n]
  e2_2=e2[2:(n+1)]
  
  sig2_1=sig2[1:n]
  sig2_2=sig2[2:(n+1)]
  # this is v_t in HR
  v=e1_2/(sig2_2^0.5)
  # bandwidth temporarily
  h=(range(v)[2]-range(v)[1])*n^(bw)
  X=cbind(1,e2_1,sig2_1)
  v_m=v%*%matrix(1,1,n, byrow = FALSE)
  temp=(v_m-t(v_m))/h
  test_kk=dnorm(temp)
  #level estimator
  g_v=1/(h*n)*colSums(dnorm(temp))
  #derivative estimator
  dg_v=1/(h^2*n)*colSums(d_n(temp))
  
  #pSi 
  psi=-(1+dg_v*v/g_v)
 
  W=0.5*X*(1/sig2_2)
  l=W*psi
  M_pF=c(0,2)
  F_t=cbind(v,(v^2-1))
  M_FF=t(F_t)%*%F_t/n
  
  P=(psi-F_t%*%t(M_pF%*%solve(M_FF)))%*%colMeans(W)
  l_star=l-P
  
  sum(colSums(l_star)^2)
}
SP_HR_ss<- function(x_ini) {
  #x_ini=init
  w <- x_ini[1]
  a <- x_ini[2]
  b <- x_ini[3]
  
  sig2=e2=e1=matrix(0,n+1,1, byrow = FALSE)
  #initial values
  sig2[1]=e2[1]=F_i2(w,a,b)
  e1[1]=e2[1]^0.5
  e1[2:(n+1)]=Y
  e2[2:(n+1)]=Y^2
  
  for (f_i in 1:n){
    sig2[1+f_i]=abs(w+a*e2[f_i]+b*sig2[f_i]) 
  }
  
  e1_2=e1[2:(n+1)]
  e2_1=e2[1:n]
  e2_2=e2[2:(n+1)]
  
  sig2_1=sig2[1:n]
  sig2_2=sig2[2:(n+1)]
  # this is v_t in HR
  v_0=e1_2/(sig2_2^0.5)
  X=cbind(1,e2_1,sig2_1)
  W0=0.5*X*(1/sig2_2)
  # sample splitting
  N=n/2
  
  #sample 1
  v=v_0[1:N]
  # bandwidth temporarily
  h=(range(v)[2]-range(v)[1])*N^(bw)
  v_m=v%*%matrix(1,1,N, byrow = FALSE)
  
  temp=(v_m-t(v_m))/h
  #test_kk=dnorm(temp)
  #level estimator
  g_v=1/(h*N)*colSums(dnorm(temp))
  #derivative estimator
  dg_v=1/(h^2*N)*colSums(d_n(temp))
  
  #pSi 
  # first we try the version without sample spliting
  # NO truncating
  psi_1=-(1+dg_v*v/g_v)
  W_1=W0[1:N,]

  M_pF=c(0,2)
  F_t_1=cbind(v,(v^2-1))
  M_FF_1=t(F_t_1)%*%F_t_1/N
  
  #sample 2
  v=v_0[(N+1):n]
  # bandwidth temporarily
  h=(range(v)[2]-range(v)[1])*N^(bw)

  v_m=v%*%matrix(1,1,N, byrow = FALSE)
  
  temp=(v_m-t(v_m))/h
  #test_kk=dnorm(temp)
  g_v=1/(h*N)*colSums(dnorm(temp))
  #derivative estimator
  dg_v=1/(h^2*N)*colSums(d_n(temp))
  
  #pSi 
  # first we try the version without sample spliting
  # NO truncating
  psi_2=-(1+dg_v*v/g_v)
  W_2=W0[(N+1):n,]
  
 M_pF=c(0,2)
  F_t_2=cbind(v,(v^2-1))
  M_FF_2=t(F_t_2)%*%F_t_2/N
  #construct the score function
  l_1=W_1*psi_2
  l_2=W_2*psi_1

  P_1=(psi_2-F_t_1%*%t(M_pF%*%solve(M_FF_1)))%*%colMeans(W_1)
  P_2=(psi_1-F_t_2%*%t(M_pF%*%solve(M_FF_2)))%*%colMeans(W_2)

  l_star_1=l_1-P_1
  l_star_2=l_2-P_2
  
  S1=sum(colSums(l_star_1)^2)
  S2=sum(colSums(l_star_2)^2)
  
  S1+S2
}

# one step HR
SP_HR_1s<- function(x_ini) {
  #x_ini=init
  w <- x_ini[1]
  a <- x_ini[2]
  b <- x_ini[3]
  
  sig2=e2=e1=matrix(0,n+1,1, byrow = FALSE)
  sig2[1]=e2[1]=F_i2(w,a,b)
  e1[1]=e2[1]^0.5
  e1[2:(n+1)]=Y
  e2[2:(n+1)]=Y^2
  
  for (f_i in 1:n){
    sig2[1+f_i]=abs(w+a*e2[f_i]+b*sig2[f_i]) 
  }
  
  e1_2=e1[2:(n+1)]
  e2_1=e2[1:n]
  e2_2=e2[2:(n+1)]
  
  sig2_1=sig2[1:n]
  sig2_2=sig2[2:(n+1)]
  # this is v_t in HR
  v=e1_2/(sig2_2^0.5)
  # bandwidth temporarily
  h=(range(v)[2]-range(v)[1])*n^(bw)
  X=cbind(1,e2_1,sig2_1)
  #kernel estimatiion
  #matrix of v 
  v_m=v%*%matrix(1,1,n, byrow = FALSE)
  temp=(v_m-t(v_m))/h
  test_kk=dnorm(temp)
  #level estimator
  g_v=1/(h*n)*colSums(dnorm(temp))
  #derivative estimator
  dg_v=1/(h^2*n)*colSums(d_n(temp))
  
  #pSi 
  psi=-(1+dg_v*v/g_v) #equation (3)
  
  W=0.5*X*(1/sig2_2)
  
  l=W*psi
  M_pF=c(0,1)
  F_t=cbind(v,(v^2-1))
  M_FF=t(F_t)%*%F_t/n
  
  P=(psi-F_t%*%t(M_pF%*%solve(M_FF)))%*%colMeans(W)
  l_star=l-P
  
  J=t(l_star)%*%l_star/n
  
  Theta=x_ini+solve(J)%*%colMeans(l_star)
  return(Theta)
  }
