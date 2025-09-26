
y1<-c(15,17,15, 13,20,15,15,13,14,17,15,17,15,18,18,15,18,10,18,18,13,16,11,16,16,18)
y2<-c(24,32 ,29,10, 26,26,26 ,22,30,30,26,28,29,32,31 ,26,33,19,30,34,30,16,25,26,23,34)
y3<-c(14,26,23,16,28,21 ,22,22,17 ,27 ,20,24,24,28,27,21 ,26,17 , 29,26,24,16,23,16,21,24)

length(y1)
length(y2)
length(y3)

# hot<-matrix(c(y1,y2,y3),ncol=3)
hot<-data.frame(y1,y2,y3)
hot

mu0<-c(14,25,20)
n<-nrow(hot)
p<-ncol(hot)
bar.y<-colMeans(hot)
bar.y

Sigma<-var(hot)
Sigma

# Population variance
n1<-length(hot)
Sigma*(n1-1)/n1

Sigma.inv<-solve(Sigma)
Sigma.inv

T2<-n*(bar.y-mu0)%*%Sigma.inv%*%(bar.y-mu0)
T2

# Decision by p-value
p.value<-pchisq(q=T2, df=p, lower.tail=FALSE)
p.value

ndf<-p
ndf
ddf<-(n-p)
ddf
fval<-qf(p=0.05,df1=ndf,df2=ddf,lower.tail=FALSE)
fval

#Decision by critical value
T2.tab<-((n-1)*p/(n-p))*fval
T2.tab



# If you are interest to solve by creating own function
mu<-c(14,25,20)
H1S<-function(x,mu){
  n<-nrow(hot)
  p<-ncol(hot)
  stopifnot(n>3*p)
  bar.y<-colMeans(hot)
  Sigma.inv<-solve(var(hot))
  T2<-n*(bar.y-mu)%*%Sigma.inv%*%(bar.y-mu)
  p.value<-pchisq(q=T2,df=p,lower.tail=FALSE)
  fval<-qf(p=0.05,df1=p,df2=n-p,lower.tail=FALSE)
  T2.tab<-((n-1)*p/(n-p))*fval
  return(list(Statistic=T2, Critical.value=T2.tab, p.value=p.value))
}
hot<-matrix(c(y1,y2,y3),ncol=3)
mu<-c(14,25,20)
H1S(hot,mu)



# Confidence Interval (CI) for mu1, mu2 and mu3
alpha<-0.05
me<-sqrt(((n-1)*p/(n-p))*qf(1-alpha,p,n-p)*diag(Sigma)/n)
cis<-cbind(bar.y-me,bar.y+me)
cat("95% simultaneous confidence interval for mu1, mu2 and mu3","\n")
cis

# Confidence Interval for the difference mu1-mu3
alpha<-0.05
med<-sqrt(((n-1)*p/(n-p))*qf(1-alpha,p,n-p)*(Sigma[1,1]*Sigma[3,3]-2*Sigma[1,3])/n)
cid<-cbind((bar.y[1]-bar.y[3])-med, (bar.y[1]-bar.y[3])+med)
cat("95% simultaneous confidence interval for mu1-mu3","\n")
cid



## Calculation Using
#install.packages("MVTests")
library(MVTests)
#install.packages("rmarkdown")
#library("rmarkdown")

hot<-data.frame(y1,y2,y3)
hot
summary(hot)
