########################################################
#### FDR (Power) vs. gamma for all types #################
########################################################
library(foreach)
library(MASS)
library(dSTEM)

#### Type I Change Point
Gamma = 3:15
## modeling
R = 2000 # repeatation
D.linebreak = function(i,signal,gamma){
  noise = rnorm(length(signal),0,1) # white noise
  data = (signal+noise)
  dy = diff(smth.gau(data,gamma),na.rm=T)
  ddy = diff(dy[!is.na(dy)])
  model2 = cpTest(x=ddy,order=2,gamma=gamma,alpha=0.05) # alpha = 0.05
  fdr_pow = Fdr(uh=c(model2$peak),th=h[-length(h)],b=10) 
  return(fdr_pow)
}
## parallel computing
c = 40 # number of cpu cores
cl = parallel::makeForkCluster(c)
doParallel::registerDoParallel(cl)
t0 = Sys.time()
out.linebreak = foreach::foreach(i = 1:length(Gamma),.packages=c("dSTEM","MASS"),.errorhandling = "pass") %dopar%{
  gamma = Gamma[i]
  l = 15000
  h = seq(150,l,150)
  jump = rep(0,length(h)+1)
  diffb = 1 # dk = {0.15 0.2 0.3}
  beta1 = seq(from=0,by=diffb,length = length(h)+1)
  signal = gen.signal(l,h,jump,beta1)
  do.call("rbind",lapply(1:R,D.linebreak,signal=signal,gamma=gamma))
}
t.linebreak = Sys.time() - t0

#### Type II Change Point (piecewise step function)
Gamma = 3:15
## modeling
R = 2000 # repeatation
D.stepjump = function(i,signal,gamma){
  noise = rnorm(length(signal),0,1) # white noise
  data = (signal+noise)
  dy = diff(smth.gau(data,gamma),na.rm=T)
  ddy = diff(dy[!is.na(dy)])
  model1 = cpTest(x=dy,order=1,alpha=0.05,gamma=gamma,is.constant=T) # alpha = 0.05
  fdr_pow = Fdr(uh=c(model1$peak),th=h[-length(h)],b=10)
  return(fdr_pow)
}
## parallel computing
t0 = Sys.time()
out.stepjump = foreach::foreach(i=1:length(Gamma),.packages=c("dSTEM","MASS"),.errorhandling = "pass") %dopar%{
  gamma = Gamma[i]
  l = 15000
  h = seq(150,l,150)
  jump = rep(2,length(h)+1) # jump = {1 1.5 2}
  beta1 = rep(0,length(h)+1)
  signal = gen.signal(l,h,jump,beta1)
  do.call("rbind",lapply(1:R,D.stepjump,signal=signal,gamma=gamma))
}
t.stepjump = Sys.time() - t0

#### Type II Change Point (piecewise linear function)
Gamma = 3:15
## modeling
R = 2000 # repeatation 
D.linejump = function(i,signal,gamma){
  noise = rnorm(length(signal),0,1) # white noise
  data = (signal+noise)
  dy = diff(smth.gau(data,gamma),na.rm=T)
  ddy = diff(dy[!is.na(dy)])
  model2 = cpTest(x=ddy,order=2,gamma=gamma,alpha=0.1)
  breaks = est.pair(model2$vall,model2$peak,gamma)$cp
  slope = est.slope(data,breaks)
  model1 = cpTest(x=dy,order=1,alpha=0.05,gamma=gamma,breaks=breaks,slope=slope)
  fdr_pow = Fdr(uh=c(model1$peak),th=h[-length(h)],b=10) 
  return(fdr_pow)
}
## parallel computing
t0 = Sys.time()
out.linejump = foreach::foreach(i = 1:length(Gamma),.packages=c("dSTEM","MASS"),errorhandling="pass") %dopar%{
  gamma = Gamma[i]
  l = 15000
  h = seq(150,l,150)
  jump = rep(3,length(h)+1) # jump = {2 2.5 3}
  addb = 0.05  # k_j + k_{j+1} = {0.02 0.05 0.1}   
  beta1 = c(rep(c(0,addb),floor(length(h)/2)),0) 
  signal = gen.signal(l,h,jump,beta1)
  do.call("rbind",lapply(1:R,D.linejump,signal=signal,gamma=gamma))
}
t2 = Sys.time() - t0

#### Mixture of Mixture of Type I and Type II Change Point  (ratio = 0.5)
Gamma = 3:15
## modeling
R = 2000 # repeatation 
D.mixture = function(i,signal,gamma){
  noise = rnorm(length(signal),0,1) # white noise
  data = (signal+noise)
  dy = diff(smth.gau(data,gamma),na.rm=T)
  ddy = diff(dy[!is.na(dy)])
  model2 = cpTest(x=ddy,order=2,gamma=gamma,alpha=0.1)
  breaks = est.pair(model2$vall,model2$peak,gamma)$cp
  if(length(breaks)==0) breaks = floor(length(data)/2)
  slope = est.slope(data,breaks)
  model1 = cpTest(x=dy,order=1,alpha=0.05,gamma=gamma,breaks=breaks,slope=slope)
  fdr_pow_jump = Fdr(uh=c(model1$peak,model1$vall),th=h[(nline+1):(length(h)-1)],b=gamma) 
  jump_seg = unique(c(sapply(c(model1$peak,model1$vall),function(x) (x-2*gamma):(x+2*gamma))))
  model2_new = cpTest(x=ddy,order=2,gamma=gamma,alpha=0.1,untest=jump_seg)
  fdr_pow_linebreak = Fdr(uh=model2_new$vall,th=h[1:nline],b=10) 
  return(cbind(fdr_pow_jump,fdr_pow_linebreak))
}
## parallel computing
t0 = Sys.time()
out.mixture = foreach::foreach(i=1:length(Gamma),.packages=c("dSTEM","MASS"),.errorhandling="pass") %dopar%{
  gamma = Gamma[i]
  ratio = 0.5 #{0.3 0.5 0.7} 
  l = 15000*2
  h = seq(150,l,150)
  njump = floor(ratio*length(h))
  nline = length(h) - njump
  jump = c(rep(0,nline+1),rep(3,njump))
  diffb = -1     # dk = {0.15 0.2 0.3}
  addb = 0.05    # addb = 0.05
  beta1 = c(seq(from=0,by=diffb,length = nline+1),rep(c(0,addb),floor(njump/2)))
  if(length(beta1)-length(h) < 1) beta1 = c(beta1,0)
  signal = gen.signal(l,h,jump,beta1)
  do.call("rbind",lapply(1:R,D.linejump,signal=signal,gamma=gamma))
}
t2 = Sys.time() - t0


