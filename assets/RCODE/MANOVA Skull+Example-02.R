rm(list = ls())  # clear the memory
###########################
## Egyptian Skull Data (Page: 349/369, Table-6.13)

data <- read.table("C:/Users/User/Downloads/datasets/T6-13.dat")
X <- as.matrix(data[,1:4])
time <- data[,5]

n1 <- sum(time==1)
n2 <- sum(time==2)
n3 <- sum(time==3)
n <- n1 + n2 + n3
ens <- c(n1,n2,n3)
p=4
g = 3

S1 = cov(X[time==1,])
S2 = cov(X[time==2,])
S3 = cov(X[time==3,])
W = (n1-1)*S1 + (n2-1)*S2 + (n3-1)*S3
Spooled = W/(n-g)

# Test equality of the Covariance Matrices 
u = (sum(1/(ens-1)) - 1/sum(ens-1))*(2*p*p + 3*p - 1)/(6*(p+1)*(g-1))
M = (n-g)*2(det(Spooled)) - 
  ((n1-1)*log(det(S1)) + (n2-1)*log(det(S2)) + (n3-1)*log(det(S3)))
C = (1-u)*M
dfC = p*(p+1)*(g-1)/2
pval = pchisq(C, dfC, lower.tail = 0)
pval # 0.3942866


xbar = colMeans(X)
tauhat1 = colMeans(X[time==1,]) - xbar
tauhat2 = colMeans(X[time==2,]) - xbar
tauhat3 = colMeans(X[time==3,]) - xbar
tauhat = cbind(tauhat1, tauhat2, tauhat3)

#Test overall treatment effects using Bartlett's approximation:
B = tauhat%*%diag(ens)%*%t(tauhat)
Wilks = det(W)/det(B+W) #  0.8301027
bart = -(n-1-((p+g)/2))*log(Wilks)
asym.pval = pchisq(bart, p*(g-1), lower.tail = 0)
asym.pval # 0.04353074

# 'matplot' plots all the columns of one matrix against all the columns of another.  
# Here we plot the 4 columns of tauhat against the one column c(1,2,3,4):

matplot(1:4, tauhat, type = 'l', xaxp = c(1,4,3), 
        lty = c(1,2,4), col = c(1,2,1), xlab = "component", ylab = "treatment effects")
legend("topright", c("tau1", "tau2", "tau3"), lty = c(1,2,4), col = c(1,2,1))


fit <- manova(X ~ as.factor(time)) #Note: time is a vector of labels - 'factors'
summary(fit)
summary(fit, test = "Wilks") # p = 0.04358
summary(fit, test = "Roy") # p = 0.005278
summary(fit, test = "Hotelling-Lawley") # p = 0.03896
summary.aov(fit)

#Bonferroni intervals on all differences
alpha = 0.05 # Some are significant at alpha = 0.05
m = p*g*(g-1)/2
q = qt(alpha/(2*m), n-g, lower.tail = 0)
confint = array(dim=c(g,g,p,3),
                dimnames = list(NULL, NULL, NULL, c("lower", "point", "upper")))
for(k in 1:g) {
  for(l in 1:g) {
    halfwidth = q*sqrt((diag(W)/(n-g))*(1/ens[k] + 1/ens[l]))
    mid = tauhat[,k] - tauhat[,l]
    lower = mid - halfwidth
    upper = mid + halfwidth
    confint[k,l, , ] = cbind(lower,mid,upper)
  }
}

for(k in 2:g) { 
  for(l in 1:(k-1)) {
  mat = round(confint[k,l, ,],3)
  sig = sign(mat[,1]*mat[,3])
  sig[sig==-1] = ' '
  sig[sig==1] = '  *'
  mat = cbind(mat, sig)
  cat(paste(100*(1-alpha),"% Bonferroni confidence intervals on tau", k, "minus tau", l, "are:","\n"))
  prmatrix(mat, quote = 0)
  cat("\n")
}}


## SObrul's Stat Talk
#----------------------------------------------------------------------#-----------------------------------#

n1<-271
n2<-138
n3<-107
n<-n1+n2+n3
n
p<-4
g<-3
x1b<-matrix(c(2.066,0.480,0.082,0.360),ncol=1)
x1b
x2b<-matrix(c(2.167,0.596,0.124,0.418),ncol=1)
x2b
x3b<-matrix(c(2.273,0.521,0.125,0.383),ncol=1)
x3b

s1<-matrix(c(0.291,-0.001,0.002,0.010,-0.001,0.011,0.000,0.003,0.002,0.000,0.001,0.000,0.010,0.003,0.000,0.010),ncol=4)
s1
s2<-matrix(c(0.561,0.011,0.001,0.037,0.011,0.025,0.004,0.007,0.001,0.004,0.005,0.002,0.037,0.007,0.002,0.019),ncol=4)
s2
s3<-matrix(c(0.261,0.030,0.003,0.018,0.030,0.017,-0.000,0.006,0.003,-0.000,0.004,0.001,0.018,0.006,0.001,0.013),ncol=4)
s3

w<-(n1-1)*s1+(n2-1)*s2+(n3-1)*s3
w
xb<-((n1*x1b)+(n2*x2b)+(n3*x3b))/n
xb
b<-(n1)*(x1b-xb)%*%t(x1b-xb)+(n2)*(x2b-xb)%*%t(x2b-xb)+(n3)*(x3b-xb)%*%t(x3b-xb)
b
wilk<-det(w)/det(b+w)
wilk

# for p>=1 and g=3
fcal<-((n-p-2)/p)*((1-sqrt(wilk))/sqrt(wilk))
ftab<-qf(0.01,2*p,2*(n-p-2),lower.tail = FALSE)
cat("Calculated value =",fcal,"\n")
cat("Tabulated vaue =",ftab, "\n")

# for large n=516 the exact test
chcal<--(n-1-(p+g)/2)*log(wilk)
chtab<-qchisq(0.01,p*(g-1),lower.tail = FALSE)
cat("Calculated value =",chcal,"\n")
cat("Tabulated vaue =",chtab, "\n")

t1h<-x1b-xb
t2h<-x2b-xb
t3h<-x3b-xb
a<-t1h[3,1]-t3h[3,1]
b<-t1h[3,1]-t2h[3,1]
c<-t2h[3,1]-t3h[3,1]
crv<-qt((0.05/(p*g*(g-1))),n-g,lower.tail = FALSE)
d<-sqrt((w[3,3]/(n-g))*((1/n1)+(1/n3)))
e<-sqrt((w[3,3]/(n-g))*((1/n1)+(1/n2)))
f<-sqrt((w[3,3]/(n-g))*((1/n2)+(1/n3)))

cat("CI for tau13-tau33 =",a-crv*d,a+crv*d,"\n","\n")
cat("CI for tau13-tau23 =",b-crv*e,b+crv*e,"\n","\n")
cat("CI for tau23-tau33 =",c-crv*f,c+crv*f,"\n","\n")

sp<-w/((n1-1)+(n2-1)+(n3-1))
m<-((n1-1)+(n2-1)+(n3-1))*log(det(sp))-((n1-1)*log(det(s1))+(n2-1)*log(det(s2))+(n3-1)*log(det(s3)))
m
u<-(1/(n1-1)+1/(n2-1)+1/(n3-1)-1/((n1-1)+(n2-1)+(n3-1)))*((2*p*p+3*p-1)/(6*(p+1)*(g-1)))
u
chcal<-(1-u)*m
chtab<-qchisq(0.05,(p*(p+1)*(g-1)/2),lower.tail = FALSE)
cat("Calculated value =",chcal,"\n","\n")
cat("Tabulated vaue =",chtab, "\n")



