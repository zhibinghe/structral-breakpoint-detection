#############################################
############ Computational Analysis##########
########### Time VS L #######################
library(dSTEM)
library(M)
##### Type I 
l = 1200
h = seq(150,by=150,length.out=6) 
jump = rep(0,7)
beta1 = c(2,-1,2.5,-3,-0.2,2.5)/50
beta1 = c(beta1,-sum(beta1*(c(h[1],diff(h))))/(l-tail(h,1)))
#modeling
Rep = 10^c(0,1,2,3,4)
Time = rep(NA,length=length(Rep)) 
for(i in 1:length(Rep)){
  signal = gen.signal(l,h,jump,beta1,rep=Rep[i])
  noise = rnorm(length(signal),0,1)
  data = (signal+noise)
  gamma = 20
  t0 = Sys.time()
  dy = diff(smth.gau(data,gamma),na.rm=T)
  ddy = diff(dy[!is.na(dy)])  
  model2 = cpTest(x=ddy,order=2,gamma=gamma,alpha=0.05)
  Time[i] = Sys.time() - t0
}
Time = cbind(l=l*Rep,time=Time)
plot(log10(Time[,2]),xaxt="n",yaxt="n",xlab="L",ylab="Time(s)",type="l",lwd=2)
points(log10(Time[,2]),pch=21,bg="red")
axis(side=1, at=1:5, labels=c(expression(10^3),expression(10^4),expression(10^5),expression(10^6),
                              expression(10^7)))
axis(side=2, at=-2:2, labels=c(0.01,0.1,1,10,100))

##### Type II
l = 1200
h = seq(150,by=150,length.out=6) 
jump = c(0,1.5,2,-2.2,1.8,2,-1.5)*2
beta1 = c(2,-1,1.5,-2.5,-0.2,2.5,-1)/50
# modeling
Rep = 10^c(0,1,2,3,4)
Time = rep(NA,length=length(Rep)) 
for(i in 1:length(Rep)){
  signal = gen.signal(l,h,jump,beta1,rep=Rep[i])
  noise = rnorm(length(signal),0,1)
  data = (signal+noise)
  gamma=10
  t0 = Sys.time()
  dy = diff(smth.gau(data,gamma),na.rm=T)
  ddy = diff(dy[!is.na(dy)])  
  # estimated location and slope
  model2 = cpTest(x=ddy,order=2,gamma=gamma,alpha=0.1)
  breaks = est.pair(model2$vall,model2$peak,gamma)$cp
  slope = est.slope(data,breaks)
  model1 = cpTest(x=dy,order=1,alpha=0.05,gamma=gamma,breaks=breaks,slope=slope)
  Time[i] = Sys.time() - t0
}
