w <- init[1]
a <- init[2]
b <- init[3]
sig2=eet=matrix(0,n+1,1, byrow = FALSE)
#sig2[1]=abs(w/(1-a-b))
sig2[1]=F_i(w,a,b)
#sig2[1]=abs(w)
eet[1]=sig2[1]^0.5

for (f_i in 1:n){
  sig2[1+f_i]=abs(w+a*eet[f_i]^2+b*sig2[f_i])
  eet[1+f_i]=Y[f_i]
}

eet=eet[2:(n+1)]/sig2[2:(n+1)]^0.5
k3 = mean(eet^3)
eta4 = mean((eet^2-1)^2)
