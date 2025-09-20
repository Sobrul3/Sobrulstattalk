x0<-rep(1,8)
x1<-c(23,23,30,30,25,25,30,30)
x2<-c(40,42,58,45,55,42,48,51)
y1<-c(41.5,33.8,27.7,21.7,19.9,15,12.2,19.3)
y2<-c(45.9,53.3,57.5,58.8,60.6,58,60.6,58)
y3<-c(11.2,11.2,12.7,16,16.2,22.6,24.5,21.3)

x<-matrix(c(x0,x1,x2),8)
x

y<-matrix(c(y1,y2,y3),8)
y

b1<-solve(t(x)%*%x)%*%t(x)%*%y1
b1
b2<-solve(t(x)%*%x)%*%t(x)%*%y2
b2
b3<-solve(t(x)%*%x)%*%t(x)%*%y3
b3

b<-matrix(c(b1,b2,b3),3)
b
# Verify
B<-solve(t(x)%*%x)%*%t(x)%*%y
B

# Y hat
yh<-x%*%b
yh

# Residuals
eh<-(y-yh)
eh

# Sum of Squares
ssr<-t(yh)%*%yh
ssr

sse<-t(eh)%*%eh
sse

sst<-t(y)%*%y
sst
ssr+sse

# Properties check
t(yh)%*%eh
t(x)%*%eh

# Solve by Regression model
m<-lm(y~cbind(x1,x2))
summary(m)
coef(m)
residuals(m)
