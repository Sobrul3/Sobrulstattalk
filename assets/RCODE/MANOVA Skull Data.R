
rm(list=ls()) # clear the memory
# Egyptian Skull Data (Page: 349/369, Table-6.13)

skull<-read.table("F:/Statistics semester wise book & sheet/7th semester/STA 4106/RCODE/T6-13.dat")

colnames(skull)<-c("MaxBreath","BasHeight","BasLength","NasHeight","TunePeriod")

head(skull)

x<-as.matrix(skull[,1:4])
time<-skull[,5]

n1<-sum(time==1)
n2<-sum(time==2)
n3<-sum(time==3)
n<-n1+n2+n3
ens<-c(n1,n2,n3)
p<-4
g<-3

s1<-cov(x[time==1,])
s2<-cov(x[time==2,])
s3<-cov(x[time==3,])

W<-(n1-1)*s1+(n2-1)*s2+(n3-1)*s3
spooled<-W/(n-g)

# Test equality of the Covariance Matrices or Box's M test
M<-(n-g)*log(det(spooled))-((n1-1)*log(det(s1))+(n2-1)*log(det(s2))+(n3-1)*log(det(s3)))
u<-(sum(1/(ens-1))-1/sum(ens-1))*(2*p*p+3*p-1)/(6*(p+1)*(g-1))
c<-(1-u)*m
df<-p*(p+1)*(g-1)/2
pval<-pchisq(c,df,lower.tail = FALSE)                                          # 0.3942866

xb=colMeans(x)
t1h<-colMeans(x[time==1,])-xb
t2h<-colMeans(x[time==2,])-xb
t3h<-colMeans(x[time==3,])-xb
th<-cbind(t1h,t2h,t3h)

B<-th%*%diag(ens)%*%t(th)
wilks<-det(W)/det(W+B)                                                         # 0.8301027
bart<--(n-1-((p+g)/2))*log(wilks)
asy.pval<-pchisq(bart,p*(g-1),lower.tail = FALSE)                              # 0.04353074


# 'matplot' plot all the columns of one matrix against all the columns of another.
# Here we plot the 4 columns of tauhat against the one column c(1,2,3,4).

matplot(1:4, th, type='l',xaxp=c(1,4,3),lty=c(1,2,4),col=c(1,2,1),xlab="Variables",ylab="Treatment Effects")
legend("topright",c("tau1","tau2","tau3"),lty=c(1,2,4),col=c(1,2,1))

fit<-manova(x~as.factor(time))   # Note: time is a vector of labels - 'factor'
summary(fit)
summary(fit, test = "Wilks")

#Critical value F with df 2p and 2*(n-p-2)
#This is the Distribution of wilks' Lambda because p=4, g=3 
crit.value <- qf (0.95, 8, 168)
crit.value

# Conclusion: with F test statistic = 2.049 > 1.9939, and a p-value of 0.04358<0.05 = ??,
# we reject our HO and conclude that the time effect differences exist.
# There is a difference of male Egyptian skulls for three different time periods.

summary(fit, test = "Pillai")
summary(fit, test = "Hotelling-Lawley")
summary(fit, test = "Roy") 
summary.aov(fit)

#Bonferroni intervals on all differences

alpha = 0.05    # Some are significant at alpha = 0.05 m = p*g*(g-1)/2
m<-p*g*(g-1)/2
q = qt(alpha/(2*m), n-g, lower.tail = FALSE)
confint<-array(dim=c(g,g,p,3),dimnames = list(NULL, NULL, NULL, c("lower", "point", "upper")))

for (k in 1:g) {
  for (l in 1:g) {
    halfwidth = q*sqrt((diag (w)/(n-g))*(1/ens [k] + 1/ens [1]))
    mid = th [, k] - th[,1]
    lower = mid - halfwidth
    upper = mid + halfwidth
    confint [k, 1,, ] = cbind(lower, mid, upper)
  }
}
for(k in 2:g){
  for (l in 1:(k-1)) {
    mat=round(confint[k,l,,],3)
    sig=sign(mat[,1]*mat[,3])
    sig[sig==-1]=''
    sig[sig==1]='*'
    mat=cbind(mat,sig)
    cat(paste(100*(1-alpha),"% Bonferrani CI on tau",k,"minus tau",l,"are:","\n"))
    prmatrix(mat,mat,quote = FALSE)
    cat("\n")
  }
}

# pair comparison
g<-3
p<-4

#All Tune Periods are 30
n1<-length (which((skull$TunePeriod ==1)))
n2<-length (which((skull$TunePeriod ==2)))
n3<-length (which((skull$TunePeriod ==3)))
n<-n1+n2+n3

#Finding means across the columns and we don't want Tune Period so we remove it!
x1b<-colMeans(skull[skull$TunePeriod==1, -5])
x2b<-colMeans(skull[skull$TunePeriod==2,-5])
x3b<-colMeans(skull[skull$TunePeriod==3, -5])
xb<-(n1*x1b+n2*x2b+n3*x3b)/(n1+n2+n3)

#Finding SE
S1<-cov(skull[skull$TunePeriod==1, -5])
S2<-cov(skull[skull$TunePeriod==2, -5])
S3<-cov(skull[skull$TunePeriod==3, -5])
W<-(n1-1)*S1+(n2-1)*S2+(n3-1)*S3

alpha<-0.05
crit.b<-qt(1-alpha/(p*g*(g-1)),n-g)

for(i in 1:p){
  LCI12<-(x1b[i]-x2b[i])-crit.b*sqrt(W[i,i]/(n-g)*(1/n1+1/n2))
  UCI12<-(x1b[i]-x2b[i])+crit.b*sqrt(W[i,i]/(n-g)*(1/n1+1/n2))
  cat("mu1[",i,"]-mu2[",i,"]=(",LCI12,",",UCI12,")\n",sep="")
  
  #\mu_{11}-\mu_{31}
  LCI13<-(x1b[i]-x3b[i])-crit.b*sqrt(W[i,i]/(n-g)*(1/n1+1/n3))
  UCI13<-(x1b[i]-x3b[i])+crit.b*sqrt(W[i,i]/(n-g)*(1/n1+1/n3))
  cat("mu1[",i,"]-mu3[",i,"]=(",LCI13,",",UCI13,")\n",sep="")
  
  #\mu_{21}-\mu_{31}
  LCI23<-(x2b[i]-x3b[i])-crit.b*sqrt(W[i,i]/(n-g)*(1/n2+1/n3))
  UCI23<-(x2b[i]-x3b[i])+crit.b*sqrt(W[i,i]/(n-g)*(1/n2+1/n3))
  cat("mu2[",i,"]-mu3[",i,"]=(",LCI23,",",UCI23,")\n",sep="")
}

#Intervals 1-3 are the 95% confidence intervals for Max Breath during the 3 different time periods
#Intervals 4-6 are the 95% confidence intervals for Bas Height during the 3 different time periods
#Intervals 7-9 are the 95% confidence intervals for Bas Length during the 3 different time periods 
#Intervals 10-12 are the 95% confidence intervals for Nas Height during the 3 different time periods

# Test Normality
#Time Period 1
S1inv<-solve(S1) # Inverse the matrix s1. 
skull1<- skull[skull$TunePeriod==1, -5]
datachisq<-diag(t(t (skull1)-x1b)%*%S1inv%*%(t(skull1)-x1b))
qqplot(qchisq(ppoints(500),df=p),datachisq, main="",xlab="Theoretical Quantiles", ylab="Sample Quantiles")

#Time Period 2
S2inv<-solve(S2) # Inverse the matrix s2. 
skull2<- skull[skull$TunePeriod==2, -5]
datachisq<-diag(t(t (skull2)-x2b)%*%S2inv%*%(t(skull2)-x2b))
qqplot(qchisq(ppoints(500),df=p),datachisq, main="",xlab="Theoretical Quantiles", ylab="Sample Quantiles")

#Time Period 3
S3inv<-solve(S3) # Inverse the matrix s3. 
skull3<- skull[skull$TunePeriod==3, -5]
datachisq<-diag(t(t (skull3)-x3b)%*%S3inv%*%(t(skull3)-x3b))
qqplot(qchisq(ppoints(500),df=p),datachisq, main="",xlab="Theoretical Quantiles", ylab="Sample Quantiles")

# abline(a = 0, b = 1, col = "red", lwd = 2)
