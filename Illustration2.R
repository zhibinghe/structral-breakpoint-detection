#####################################################
############# Main Idea Illustration ################
############# Illustration figure 2 and 2 ###########
#####################################################
library(dSTEM) 
set.seed(2021)

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



