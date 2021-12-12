#####################################################
################ Illustration figure 1 ##############
#####################################################
library(latex2exp)
library(dSTEM)
##
k1 = 0.08
k2 = 0.15
gamma=8
v = 100
c = 4 # bandwidth of smoothing kernel
d1mu = function(x) a/gamma*dnorm((v-x)/gamma)+(k1-k2)*pnorm((v-x)/gamma)+k2
d2mu = function(x) (a*(v-x)+(k2-k1)*gamma^2)/gamma^3*dnorm((v-x)/gamma)
par(mfrow=c(3,3),mar=c(2,2.5,2,0.5))
#### Type I change Point
a = 0
b = a - (k2-k1)*v
xmax = 200
ymax = k2*xmax+b
x = 1:xmax
curve(k1*x,0,floor(xmax/2),xlim = c(0,xmax),ylim=c(0,ymax+1),lwd=2,xlab="",ylab="",main=TeX("$\\mu(t)$"))
curve(k2*x+b,floor(xmax/2),xmax,add=TRUE,lwd=2,xlab="",ylab="")
points(x=v,y=k1*v,pch=16)
points(x=v,y=k2*v+b,pch=1)
abline(v=v,col="red",lty=2)
y = d1mu(x)
plot(loess(y~x,span=0.6),type="l",lwd=2,xlab="",ylab="",main=TeX("$\\mu_{\\gamma}^{\\prime}(t)$"))
abline(v=v,col="red",lty=2,lwd=1)
points(x = c(v-c*gamma,v+c*gamma),y=y[c(v-c*gamma,v+c*gamma)],pch =16 ,col="blue")
y = d2mu(x)
plot(loess(y~x,span=0.6),type="l",lwd=2,xlab="",ylab="",main=TeX("$\\mu_{\\gamma}^{\\prime\\prime}(t)$"))
abline(v=v,col="red",lty=2,lwd=1)
points(x = c(v-c*gamma,v+c*gamma),y=y[c(v-c*gamma,v+c*gamma)],pch =16 ,col="blue")
#### Type II Change Point
a = 3
b = a - (k2-k1)*v
xmax = 200
ymax = k2*xmax+b
x = 1:xmax
curve(k1*x,0,floor(xmax/2),xlim = c(0,xmax),ylim=c(0,ymax+1),lwd=2,xlab="",ylab="")
curve(k2*x+b,floor(xmax/2),xmax,add=TRUE,lwd=2,xlab="",ylab="")
points(x=v,y=k1*v,pch=16)
points(x=v,y=k2*v+b,pch=1)
abline(v=v,col="red",lty=2)
y = d1mu(x)
plot(loess(y~x,span=0.6),type="l",lwd=2,xlab="",ylab="")
abline(v=v,col="red",lty=2,lwd=1)
points(x = c(v-c*gamma,v+c*gamma),y=y[c(v-c*gamma,v+c*gamma)],pch =16 ,col="blue")
y = d2mu(x)
plot(loess(y~x,span=0.6),type="l",lwd=2,xlab="",ylab="")
abline(v=v,col="red",lty=2,lwd=1)
points(x = c(v-c*gamma,v+c*gamma),y=y[c(v-c*gamma,v+c*gamma)],pch =16 ,col="blue")
#### Mixture of Type I and Type II Change Points
a = 0
b = a - (k2-k1)*v
xmax1 = 200 ; xmax2 = 300
ymax = k2*xmax1+b
x = 1:xmax2
a3 = -3 ; k3 = 0 ; v3 = xmax1
b3 = a3 - (k3-k2)*v3 +b
curve(k1*x,0,floor(xmax1/2),xlim = c(0,xmax2),ylim=c(0,ymax+1),lwd=2,xlab="",ylab="")
curve(k2*x+b,floor(xmax1/2),xmax1,add=TRUE,lwd=2,xlab="",ylab="")
points(x=v,y=k1*v,pch=16)
abline(v=v,col="red",lty=2)
curve(k3*x+b3,xmax1,xmax2,add=TRUE,lwd=2,xlab="",ylab="")
abline(v=v3,col="red",lty=2)
points(x=v3,y=k2*v3+b,pch=16)
points(x=v3,y=k3*v3+b3,pch=1)
#
y1 = d1mu(x[1:floor(xmax2/2)])
plot(loess(y1~x[1:floor(xmax2/2)],span=0.6),type="l",lwd=2,xlab="",ylab="",xlim=c(0,xmax2),ylim=c(-0.10,0.15))
abline(v=v,col="red",lty=2,lwd=1)
points(x = c(v-c*gamma,v+c*gamma),y=y1[c(v-c*gamma,v+c*gamma)],pch =16 ,col="blue")
k1=k2;k2=k3;a=a3;v=v3
y2 = d1mu(x[floor(xmax2/2):xmax2])
lines(loess(y2~x[floor(xmax2/2):xmax2],span=0.6),lwd=2,xlab="",ylab="")
abline(v=v3,col="red",lty=2,lwd=1)
points(x = c(v3-c*gamma,v3+c*gamma),y=y2[c(v3-c*gamma,v3+c*gamma)-floor(xmax2/2)],pch =16 ,col="blue")
#
k1 = 0.08;k2 = 0.15;v = 100;a = 0
y1 = d2mu(x[1:floor(xmax2/2)])
plot(loess(y1~x[1:floor(xmax2/2)],span=0.6),type="l",lwd=2,xlab="",ylab="",xlim=c(0,xmax2),ylim=c(-0.016,0.007))
abline(v=v,col="red",lty=2,lwd=1)
points(x = c(v-c*gamma,v+c*gamma),y=y1[c(v-c*gamma,v+c*gamma)],pch =16 ,col="blue")
k1=k2;k2=k3;a=a3;v=v3
y2 = d2mu(x[floor(xmax2/2):xmax2])
lines(loess(y2~x[floor(xmax2/2):xmax2],span=0.6),lwd=2,xlab="",ylab="")
abline(v=v3,col="red",lty=2,lwd=1)
points(x = c(v3-c*gamma,v3+c*gamma),y=y2[c(v3-c*gamma,v3+c*gamma)-floor(xmax2/2)],pch =16 ,col="blue")

