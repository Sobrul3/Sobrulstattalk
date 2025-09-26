
z0<-rep(1,5)
z1<-c(0,1,2,3,4)
y1<-c(1,4,3,8,9)
y2<-c(-1,-1,2,3,2)

y<-matrix(c(y1,y2),ncol=2)
y

z<-matrix(c(z0,z1),ncol=2)
z
tz.z<-t(z)%*%z
tz.z

tz.z.inv<-solve(tz.z)
tz.z.inv

tz.y1<-t(z)%*%y1
tz.y1

b1<-solve(t(z)%*%z)%*%t(z)%*%y1
b1

tz.y2<-t(z)%*%y2
tz.y2

b2<-solve(t(z)%*%z)%*%t(z)%*%y2
b2

b<-matrix(c(b1,b2),ncol=2)
b

y.hat<-z%*%b
y.hat

e.hat<-(y-y.hat)
e.hat


# Properties Check
sum(e.hat)
t(y.hat)%*%e.hat
t(z)%*%e.hat

ssr<-t(y.hat)%*%y.hat
ssr

sse<-t(e.hat)%*%e.hat
sse

sst<-t(y)%*%y
sst
ssr+sse


# By Regreesion model
m<-lm(y~z1)
summary(m)

