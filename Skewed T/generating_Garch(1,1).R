sim = 103+j
set.seed(sim)
#data generating
#eta=rnorm(n+n_start,0,1)
#eta=rsnorm(n+n_start, mean = 0, sd = 1, xi = 2)

t_var=t_df/(t_df-2)
#eta=rt(n+n_start, df = t_df)* sqrt(1/ t_var)
eta=rsstd(n+n_start,mean=0, sd=1, nu = t_df, xi = 2)

w0 <- theta[1]
a0 <- theta[2]
b0 <- theta[3]

Y=sig2_0=eet_0=matrix(0,n+n_start,1, byrow = FALSE)
#sig2_0[1]=abs(w0/(1-a0-b0))
sig2_0[1]=F_i(w0,a0,b0)
#sig2_0[1]=abs(w0)
eet_0[1]=sig2_0[1]^0.5*eta[1]

for (f_i in 1:(n+n_start-1)){
  sig2_0[1+f_i]=w0+a0*eet_0[f_i]^2+b0*sig2_0[f_i]
  eet_0[1+f_i]=sig2_0[1+f_i]^0.5*eta[f_i+1]
  Y[1+f_i]=eet_0[1+f_i]
}
Y=Y[(1+n_start):(n+n_start)]