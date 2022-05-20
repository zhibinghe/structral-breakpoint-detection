#####################################################
############## Simulation ###########################
###### Comparason with other methods ################
#####################################################
source("main.R")
library(foreach)
################ PART I #############################
############# Short term data sequence ##############

x = signal.info(type="I")
signal = gen.signal(l=x$n,h=x$cpt,jump = x$jump.size,b1 = x$slope)
data = signal + rnorm(length(signal))

# run 100 times
cpt.dstem = unlist(comp.detection(data,method="dstem","I",gamma=10))
cpt.not = comp.detection(data,method="not","I",M=0.1*length(data))
cpt.nsp = comp.detection(data,method="nsp","I",M=0.1*length(data))
cpt.bp = comp.detection(data,method="bp","I",M=20)



####################################################
gamma=15
#### piecewise linear function
t0 = Sys.time()
l = 1200
h = seq(150,by=150,length.out=6) 
jump = rep(0,7)
beta1 = c(2,-1,2.5,-3,-0.2,2.5)/25
beta1 = c(beta1,-sum(beta1*(c(h[1],diff(h))))/(l-tail(h,1)))
signal = gen.signal(l,h,jump,beta1)
## Modeling
R = 100 #repeatation
D.linebreak = function(i){
  noise = rnorm(length(signal),0,1) # white noise
  data = signal+noise
  # modelling
  dy = diff(smth.gau(data,gamma),na.rm=T)
  ddy = diff(dy[!is.na(dy)])
  model2 = cpTest(x=ddy,order=2,gamma=gamma,alpha=0.05) # alpha = 0.05
  return(sort(c(model2$peak,model2$vall)))
}
# parallel computing
c = 5 # number of cpu cores
cl =parallel::makeCluster(c)
doParallel::registerDoParallel(cl)

out_linebreak_dstem = foreach(iterators::icount(R),.packages="dSTEM",.errorhandling="pass") %dopar% {
  D.linebreak(0)
}
time_linebreak_dstem = Sys.time()-t0

################################
## piecewise constant function
t0 = Sys.time()
l = 1200
h = seq(150,by=150,length.out=6) 
jump = c(0,1.5,2,2.2,1.8,2,1.5)
beta1 = rep(0,length(h)+1)
signal = gen.signal(l,h,jump,beta1)
#modeling
R = 100 #repeatation
D.stepjump = function(i){
  noise = rnorm(length(signal),0,1) # white noise
  data = signal+noise
  # modelling
  dy = diff(smth.gau(data,gamma),na.rm=T)
  ddy = diff(dy[!is.na(dy)])
  model1 = cpTest(x=dy,order=1,alpha=0.05,gamma=gamma,is.constant=T) # alpha = 0.05
  return(sort(c(model1$peak,model1$vall)))
}
out_stepjump_dstem = foreach(iterators::icount(R),.packages="dSTEM",.errorhandling="pass") %dopar% {
  D.stepjump(0)
}
time_stepjump_dstem = Sys.time()-t0
#############################
## Normal linear breaks
t0 = Sys.time()
l = 1200
h = seq(150,by=150,length.out=6) 
jump = c(0,1.5,2,2.2,1.8,2,1.5)*5
beta1 = c(2,-1,2.5,-3,-0.2,2.5,-0.5)/25
signal = gen.signal(l,h,jump,beta1)
#modeling
R = 100 #repeatation
D.linejump = function(i){
  noise = rnorm(length(signal),0,1) # white noise
  data = signal+noise
  dy = diff(smth.gau(data,gamma),na.rm=T)
  ddy = diff(dy[!is.na(dy)])
  # modelling with true slopes
  #model1 = cpTest(x=dy,order=1,alpha=0.05,gamma=gamma,breaks=h,slope=beta1)
  # estimated location and slope
  model2 = cpTest(x=ddy,order=2,gamma=gamma,alpha=0.1)
  breaks = est.pair(model2$vall,model2$peak,gamma)$cp
  slope = est.slope(data,breaks)
  model1 = cpTest(x=dy,order=1,alpha=0.05,gamma=gamma,breaks=breaks,slope=slope)
  return(sort(c(model1$peak,model1$vall)))
}
out_linejump_dstem = foreach(iterators::icount(R),.packages="dSTEM",.errorhandling="pass") %dopar% {
  D.linejump(0)
}
time_linejump_dstem = Sys.time()-t0
parallel::stopCluster(cl)
