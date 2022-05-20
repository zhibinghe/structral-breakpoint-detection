#####################################################
################ Illustration figure 1 ##############
#####################################################
library(latex2exp)
library(dSTEM)
set.seed(2021)
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
#####################################################
############# Illustration figure 2 ###########
#####################################################
#### Type I change points
l = 1200
h = seq(150,by=150,length.out=6) 
jump = rep(0,7)
beta1 = c(2,-1,2.5,-3,-0.2,2.5)/50
beta1 = c(beta1,-sum(beta1*(c(h[1],diff(h))))/(l-tail(h,1)))
signal = gen.signal(l,h,jump,beta1)
#signal = gen.signal(l,h,jump,beta1,rep=10)
noise = rnorm(length(signal),0,1)
data = (signal+noise)
#modeling
gamma = 20
dy = diff(smth.gau(data,gamma),na.rm=T)
ddy = diff(dy[!is.na(dy)])  
model1_2 = cpTest(x=ddy,order=2,gamma=gamma,alpha=0.05)
##
par(mfrow = c(2,2),mar=c(4,4,2,2))
b = 8
##
col1 = "blue"; col2 = "darkorange"
plot(data,col="grey",xlab="",ylab="",bty="n")
title(xlab="t",ylab="y(t)",mgp=c(2.5,2.5,0),cex.lab=1.0)
lines(signal,lwd=2,col="red")
plot(ddy,type="l",lwd=2,xlab="",ylab="",bty="n")
title(xlab="t",ylab=latex2exp::TeX("$y''_{\\gamma}(t)"),mgp=c(2.5,2.5,0),cex.lab=1.0)
abline(h=0,lty=1,lwd=2,col="lightgrey")
peak = which.peaks(ddy); val = which.peaks(ddy,decreasing = T)
points(peak,ddy[peak],col=col1,pch=16,bg=col1)
points(val,ddy[val],col=col2,pch=16,bg=col2)
tpeak = model1_2$peak; tval = model1_2$vall
points(tpeak,ddy[tpeak],col=col1,pch=24,bg=col1)
points(tval,ddy[tval],col=col2,pch=25,bg=col2)
Tmax = h[which(h %in% which.peaks(signal,decreasing=T))]; Tmin = h[which(h %in% which.peaks(signal))]
for(i in 1:length(Tmax)) rect(Tmax[i]-b,0,Tmax[i]+b,1,col="cyan",lwd=1.5,density=0)
for(i in 1:length(Tmin)) rect(Tmin[i]-b,0,Tmin[i]+b,-1,col="pink",lwd=1.5,density=0)

#### Type II chnage points
l = 1200
h = seq(150,by=150,length.out=6) 
jump = c(0,1.5,2,-2.2,1.8,2,-1.5)*2
beta1 = c(2,-1,1.5,-2.5,-0.2,2.5,-1)/50
signal = gen.signal(l,h,jump,beta1)
noise = rnorm(length(signal),0,1)
data = (signal+noise)
# modeling
gamma=10
dy = diff(smth.gau(data,gamma),na.rm=T)
ddy = diff(dy[!is.na(dy)])  
# estimated location and slope
model2_2 = cpTest(x=ddy,order=2,gamma=gamma,alpha=0.1)
breaks = est.pair(model2_2$vall,model2_2$peak,gamma)$cp
slope = est.slope(data,breaks)
model2_1 = cpTest(x=dy,order=1,alpha=0.05,gamma=gamma,breaks=breaks,slope=slope)
##
col3="green"; col4 = "red"
plot(data,col="grey",xlab="",ylab="",bty="n")
title(xlab="t",ylab="y(t)",mgp=c(2.5,2.5,0),cex.lab=1.0)
#
segment = split(1:l,cut(1:l,breaks=c(1,h,l),include.lowest = T))
for(i in 1:length(segment)) lines(x=segment[[i]],y=signal[segment[[i]]],lwd=2,col="red")
#
plot(dy,type="l",lwd=2,xlab="",ylab="",bty="n")
title(xlab="t",ylab=latex2exp::TeX("$y'_{\\gamma}(t)"),mgp=c(2.5,2.5,0),cex.lab=1.0)
abline(h=0,lty=1,lwd=2,col="lightgrey")
peak = which.peaks(dy); val = which.peaks(dy,decreasing = T)
points(peak,dy[peak],col=col3,pch=16,bg=col3)
points(val,dy[val],col=col4,pch=16,bg=col4)
tpeak = model2_1$peak; tval = model2_1$vall
points(tpeak,dy[tpeak],col=col3,pch=24,bg=col3)
points(tval,dy[tval],col=col4,pch=25,bg=col4)
Tmin = h[c(3,6)]; Tmax = h[c(1,2,4,5)]
for(i in 1:length(Tmax)) rect(Tmax[i]-b,0,Tmax[i]+b,1,col="cyan",lwd=1.5,density=0)
for(i in 1:length(Tmin)) rect(Tmin[i]-b,0,Tmin[i]+b,-1,col="pink",lwd=1.5,density=0)
for(i in 1:length(segment)) lines(x=segment[[i]],y=rep(beta1[i],length(segment[[i]])),
                                  lwd=2,col="darkgrey")
##### Mixture of type I and type II
l = 1500
## contains only 8 (6) change points
h = seq(150,by=150,length.out=8)
jump = c(0,0,2.5,0,3,0,0,2.8,-2.5)
beta1 = c(1,-1,1.5,-2,-0.2,1.2,-0.5,1.2,0)/50
signal = gen.signal(l,h,jump,beta1)
noise = rnorm(length(signal),0,1)
data = (signal+noise)
## Modeling
gamma = 15
dy = diff(smth.gau(data,gamma),na.rm=T)
ddy = diff(dy[!is.na(dy)])  
model2 = cpTest(x=ddy,order=2,gamma=gamma,alpha=0.1)
breaks = est.pair(model2$vall,model2$peak,gamma)$cp
slope = est.slope(data,breaks)
model1 = cpTest(x=dy,order=1,alpha=0.05,gamma=gamma,breaks=breaks,slope=slope)
model2.new = cpTest(x=ddy,order=2,gamma=gamma,alpha=0.1,untest=c(model1$peak,model1$vall))
##
par(mfrow = c(1,3),mar=c(4,4,2,2))
plot(data,col="grey",xlab="",ylab="",bty="n")
title(xlab="t",ylab="y(t)",mgp=c(2.5,2.5,0),cex.lab=1.0)
#
segment = split(1:l,cut(1:l,breaks=c(1,h[c(2,4,7,8)],l),include.lowest = T))
for(i in 1:length(segment)) lines(x=segment[[i]],y=signal[segment[[i]]],lwd=2,col="red")
#
b=15
plot(dy,type="l",lwd=2,xlab="",ylab="",bty="n")
title(xlab="t",ylab=latex2exp::TeX("$y'_{\\gamma}(t)"),mgp=c(2.5,2.5,0),cex.lab=1.0)
abline(h=0,lty=1,lwd=2,col="lightgrey")
peak = which.peaks(dy); val = which.peaks(dy,decreasing = T)
points(peak,dy[peak],col=col3,pch=16,bg=col3)
points(val,dy[val],col=col4,pch=16,bg=col4)
tpeak = model1$peak; tval = model1$vall
points(tpeak,dy[tpeak],col=col3,pch=24,bg=col3)
points(tval,dy[tval],col=col4,pch=25,bg=col4)
Tmin = h[8]; Tmax = h[c(2,4,7)]
for(i in 1:length(Tmax)) rect(Tmax[i]-b,0,Tmax[i]+b,1,col="cyan",lwd=1.5,density=0)
for(i in 1:length(Tmin)) rect(Tmin[i]-b,0,Tmin[i]+b,-1,col="pink",lwd=1.5,density=0)
#
segment = split(1:l,cut(1:l,breaks=c(1,h,l),include.lowest = T))
for(i in 1:length(segment)) lines(x=segment[[i]],y=rep(beta1[i],length(segment[[i]])),
                                  lwd=2,col="darkgrey")
#
plot(ddy,type="l",lwd=2,xlab="",ylab="",bty="n")
title(xlab="t",ylab=latex2exp::TeX("$y''_{\\gamma}(t)"),mgp=c(2.5,2.5,0),cex.lab=1.0)
abline(h=0,lty=1,lwd=2,col="lightgrey")
peak = which.peaks(ddy); val = which.peaks(ddy,decreasing = T)
points(peak,ddy[peak],col=col1,pch=16,bg=col1)
points(val,ddy[val],col=col2,pch=16,bg=col2)
tpeak = model2.new$peak; tval = model2.new$vall
points(tpeak,ddy[tpeak],col=col1,pch=24,bg=col1)
points(tval,ddy[tval],col=col2,pch=25,bg=col2)
Tmin = h[c(1,3,6)]; Tmax = h[5]
for(i in 1:length(Tmax)) rect(Tmax[i]-b,0,Tmax[i]+b,1,col="cyan",lwd=1.5,density=0)
for(i in 1:length(Tmin)) rect(Tmin[i]-b,0,Tmin[i]+b,-1,col="pink",lwd=1.5,density=0)
Tmax = c(model1$peak,model1$vall)
for(i in 1:length(Tmax)) rect(Tmax[i]-2*gamma,0,Tmax[i]+1.5*gamma,1,col="grey",lwd=1.5,density=10)
for(i in 1:length(Tmax)) rect(Tmax[i]-2*gamma,0,Tmax[i]+1.5*gamma,-1,col="grey",lwd=1.5,density=10)




