#Problem-01

x<-c(1250,60,1230,60,900,415,1230,415,1920)
MulCol<-function(x,n,k){
  e<-matrix(x,ncol=k)
  r<-matrix(rep(0,k*k),ncol=k)
  for(i in 1:k){
    for(j in 1:k){
  r[i,j]<-(e[i,j]/sqrt(e[i,i]*e[j,j]))
  j<-j+1
    }
    i<-i+1
  }
  xtx<-matrix(r,ncol=k)
  print(xtx)
  dxtx<-det(xtx)
  print(dxtx)
  ch<--(n-1-((2*k+5)/6))*log(dxtx)
  v<-(k*(k-1)/2)
  cat("Calculated value = ",ch,"\n")
  tab<-qchisq(0.05,v,lower.tail = FALSE)
  cat("Tabulated value = ",tab,"\n")
  R1<-(r[1,2]^2+r[1,3]^2-2*r[1,2]*r[1,3]*r[2,3])/(1-r[2,3]^2)
  R2<-(r[1,2]^2+r[2,3]^2-2*r[1,2]*r[1,3]*r[2,3])/(1-r[1,3]^2)
  R3<-(r[1,3]^2+r[2,3]^2-2*r[1,2]*r[1,3]*r[2,3])/(1-r[1,2]^2)
  vf1<-1/(1-R1)
  vf2<-1/(1-R2)
  vf3<-1/(1-R3)
  print(vf1)
  print(vf2)
  print(vf3)
}

##-------------------------------------------------------------------------------------------------------------------##
# Problem-02

ob<-2000:2009
ob<-as.character(ob)
y<-c(6.1,6.3,6.5,7.2,7.5,7.6,8.0,8.5,9.0,9.4)
x1<-c(45,46.3,47.5,48.8,49.3,51.6,56,60.2,64.5,66)
x2<-c(5.8,4.5,5.5,6.4,7.2,8.5,10.7,25,17.6,22)
x3<-c(106,98,108,102,96,97,100,98,93,103)
AA<-(x1-mean(x1))^2
BB<-(x2-mean(x2))^2
CC<-(x3-mean(x3))^2
AB<-(x1-mean(x1))*(x2-mean(x2))
AC<-(x1-mean(x1))*(x3-mean(x3))
BC<-(x2-mean(x2))*(x3-mean(x3))

dat<-data.frame("Observations"=ob,"A^2"=AA,"B^2"=BB,"C^2"=CC,"AB"=AB,"AC"=AC,"BC"=BC)
dat

#----------------------------------------------------------------#
ob<-c(2000:2009)
ob<-as.character(ob)
ob<-c(ob,"Total")
y<-c(6.1,6.3,6.5,7.2,7.5,7.6,8.0,8.5,9.0,9.4)
x1<-c(45,46.3,47.5,48.8,49.3,51.6,56,60.2,64.5,66)
x2<-c(5.8,4.5,5.5,6.4,7.2,8.5,10.7,25,17.6,22)
x3<-c(106,98,108,102,96,97,100,98,93,103)
AA<-c((x1-mean(x1))^2,sum((x1-mean(x1))^2))
BB<-c((x2-mean(x2))^2,sum((x2-mean(x2))^2))
CC<-c((x3-mean(x3))^2,sum((x3-mean(x3))^2))
AB<-c((x1-mean(x1))*(x2-mean(x2)),sum((x1-mean(x1))*(x2-mean(x2))))
AC<-c((x1-mean(x1))*(x3-mean(x3)),sum((x1-mean(x1))*(x3-mean(x3))))
BC<-c((x2-mean(x2))*(x3-mean(x3)),sum((x2-mean(x2))*(x3-mean(x3))))

dat<-data.frame("Observations"=ob,"A^2"=AA,"B^2"=BB,"C^2"=CC,"AB"=AB,"AC"=AC,"BC"=BC)
dat
#--------------------------------------------------------------#

r12<-sum(AB)/sqrt(sum(AA)*sum(BB))
r13<-sum(AC)/sqrt(sum(AA)*sum(CC))
r23<-sum(BC)/sqrt(sum(BB)*sum(CC))
xtx<-matrix(c(1,r12,r13,r12,1,r23,r13,r23,1),ncol=3)
xtx
dxtx<-det(xtx)
dxtx
ch<--(10-1-((2*3+5)/6))*log(dxtx)
v<-(3*(3-1)/2)
cat("Calculated value = ",ch,"\n")
tab<-qchisq(0.05,v,lower.tail = FALSE)
cat("Tabulated value = ",tab,"\n")

R1<-(r12^2+r13^2-2*r12*r13*r23)/(1-r23^2)
R1
R2<-(r12^2+r23^2-2*r12*r13*r23)/(1-r13^2)
R2
R3<-(r13^2+r23^2-2*r12*r13*r23)/(1-r12^2)
R3
R<-c(R1,R2,R3)
for(i in 1:3){
  s<-(R[i]/2)/((1-R[i])/7)
  t<-qf(0.05,2,7,lower.tail = FALSE)
  cat("Calculated Value ","=",i,":",s,"\n")
  cat("Tabulated Value =",t,"\n","\n")
}

rr12<-(r12-(r13*r23))^2/((1-r23^2)*(1-r13^2))
rr13<-(r13-(r12*r23))^2/((1-r23^2)*(1-r12^2))
rr23<-(r23-(r12*r13))^2/((1-r13^2)*(1-r12^2))
r<-c(rr12,rr13,rr23)
for(i in 1:3){
  s<-sqrt((r[i]*7)/(1-r[i]))
  t<-qt(0.05,7,lower.tail = FALSE)
  cat("Calculated Value ","=",i,":",s,"\n")
  cat("Tabulated Value =",t,"\n","\n")
}


##--------------------------------------------------------------------------------------------------------------------##
# Problem-03

ob<-c(seq(1959,1980,1))
ob<-as.character(ob)
RC<-c(58.5,59.9,61.7,63.9,65.3,67.8,69.3,71.8,73.7,76.5,77.6,79,80.5,82.9,84.7,83.7,84.5,87,88.1,89.7,90,89.7)
P<-c(47.2,48,49.8,52.1,54.1,54.6,58.6,61,62.3,64.5,64.8,66.2,68.8,71,73.1,72.2,74.8,77.2,78.4,79.5,79.7,79.8)

dat<-data.frame("Observations"=ob,"Real Compensation"=RC,"Productivity"=P,"P^2"=P*P,"RC.P"=RC*P)
dat
bh<-((length(P)*sum(RC*P))-(sum(RC)*sum(P)))/(length(P)*sum(P^2)-(sum(P))^2)
ah<-mean(RC)-bh*mean(P)
RCh<-ah+(bh*P)
e<-RC-RCh
ee1<-e[-length(e)]*e[-1]

dat<-data.frame("Observations"=ob,"Real Compensation"=RC,"Productivity"=P,"P^2"=P*P,"RC.P"=RC*P,"RC hat"=RCh,"Error"=round(e,2),"Error*Error1"=round(ee1,2))
dat

plot(P,e,xlab="Observation",ylab="Residual",pch=16,col="blue")
lines(P[order(P)],e[order(P)],lty=2,col="red")

d<-2*(1-(sum(ee1)/sum(e*e)))


##----------------------------------------------------------------------------------------------------------------------------------------------------##
# Problem-04

x<-c(49,55,55,70,53,70,55,62,62,80,73,92,92,73,66,73,78,92,78,90,93,74,95,81,94,97,84,98,99,57)
y<-c(65.4,56,55.9,49,46.5,46.2,45.4,59.2,53.3,43.4,41.1,40.9,40.9,40.4,39.6,39.3,38.9,38.8,38.2,42.2,40.9,40.7,40,39.3,38.8,38.4,38.4,38.4,46.9,36.3)

dat<-data.frame("X"=x, "Y"=y)
dat

plot(x,y,xlab="Engine Horsepower",ylab="Miles Per Gallon",pch=16,col="green")
abline(lm(y~x), col="red",lwd=2)

xod<-order(x)
x11<-xod[1:11]
x22<-xod[20:30]

x1<-x[x11]
y1<-y[x11]
x1
y1

dat<-data.frame("X"=x1,"Y"=y1,"X^2"=x1^2,"XY"=x1*y1)
dat

bh<-((length(x1)*sum(y1*x1))-(sum(x1)*sum(y1)))/(length(x1)*sum(x1^2)-(sum(x1))^2)
ah<-mean(y1)-bh*mean(x1)
y1h<-ah+(bh*x1)
e<-y1-y1h

dat<-data.frame("X"=x1,"Y"=y1,"X^2"=x1^2,"XY"=x1*y1,"Y hat"=y1h,"Error"=round(e,2),"Error^2"=round(e*e,2))
dat

sse1<-sum(e*e)


x2<-x[x22]
y2<-y[x22]
x2
y2

dat<-data.frame("X"=x2,"Y"=y2,"X^2"=x2^2,"XY"=x2*y2)
dat

bh<-((length(x2)*sum(y2*x2))-(sum(x2)*sum(y2)))/(length(x2)*sum(x2^2)-(sum(x2))^2)
ah<-mean(y2)-bh*mean(x2)
y2h<-ah+(bh*x2)
e<-y2-y2h

dat<-data.frame("X"=x2,"Y"=y2,"X^2"=x2^2,"XY"=x2*y2,"Y hat"=y2h,"Error"=round(e,2),"Error^2"=round(e*e,2))
dat

sse2<-sum(e*e)

v<-((length(x)-8)/2)-2
fcal<-sse2/sse1
cat("Calculated Value =", fcal,"\n")
ftab<-qf(0.05,v,v,lower.tail = FALSE)
cat("Tabulated Value =", ftab,"\n")


dat<-data.frame("X"=x,"Y"=y,"X^2"=x^2,"XY"=x*y)
dat

bh<-((length(x)*sum(y*x))-(sum(x)*sum(y)))/(length(x)*sum(x^2)-(sum(x))^2)
ah<-mean(y)-bh*mean(x)
yh<-ah+(bh*x)
e<-y-yh
e<-round(e,2)
d2<-(rank(x)-rank(abs(e)))^2

dat<-data.frame("X"=x,"Y"=y,"X^2"=x^2,"XY"=x*y,"Y hat"=yh,"Error"=e,"|Error|"=abs(e),"Rank X"=rank(x),"Rank Error"=rank(abs(e)),"D^2"=d2)
dat

rs<-1-(6*(sum(d2))/(length(x)*((length(x))^2-1)))

rcal<-(rs*sqrt(length(x)-2))/sqrt(1-rs^2)
rtab<-qt(0.025,length(x)-2, lower.tail = FALSE)

cat("Absolute Calculated value =", abs(rcal),"\n")
cat("Tabulated value =", rtab,"\n")

#SOBRUL'S STAT TALK
