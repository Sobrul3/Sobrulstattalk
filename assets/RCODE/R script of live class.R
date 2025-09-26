Y1=c(15,17,15,13,20,15,15,13,14,17,15,17,15,18,18,15,18,10,18,18,13,16,11,16,16,18)
Y2=c(24,32,29,10,26,26,26,22,30,30,26,28,29,32,31,26,33,19,30,34,30,16,25,26,23,34)
Y3=c(14,26,23,16,28,21,22,22,17,27,20,24,24,28,27,21,26,17,29,26,24,16,23,16,21,24)
data= data.frame(Y1,Y2,Y3)
mu0=c(14,25,20)
n=nrow(data)
p=ncol(data)
ybar=colMeans(data)
#calculate sample variance
sigma=var(data)
sigma_inverse=solve(sigma)
T2=n*(ybar- mu0)%*%sigma_inverse%*%(ybar-mu0)
T2
p_value=pchisq(q=T2,df=p,lower.tail = FALSE)
p_value

Fval=qf(p=0.05,df1=p,df2=(n-p),lower.tail = FALSE)
Fval
T2_tab=((n-1)*p/(n-p))*Fval
T2_tab
T2
#similtanous confidence interval
wd=sqrt((n-1)*p/(n-p)*Fval)*sqrt(diag(sigma)/n)
wd
cls=cbind(ybar-wd,ybar+wd)
cls
