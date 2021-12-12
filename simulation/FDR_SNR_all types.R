########################################################
### FDR (Power) vs. SNR(slope change) for all types ####
########################################################
library(foreach)
library(MASS)
library(dSTEM)

#### Calculate slope difference based on SNR 
inverse.snr = function(Snr,order,gamma,addslope){
  # order: the order of derivative
  # addslope: k_j + k_{j+1}, only required for order = 1
  '%notin%' <- Negate('%in%')
  if (order %notin% c(1,2)) stop("order must be 1 or 2")
  if (missing(addslope) & order==1) stop("addslope is required when order is 1")
  if (order==1) {
    jump = (Snr - addslope*pi^(1/4)*gamma^(3/2))*pi^(1/4)/(sqrt(2*gamma))
    return(jump)
  }
  else {
    minusslope = Snr*sqrt(3)*pi^(1/4)/(2*gamma^(3/2))
    return(minusslope)
  }
}
#### Type I Change Point
gamma = 10
## modeling
R = 2000 # repetition
D.linebreak = function(i,signal){
  noise = rnorm(length(signal),0,1)
  data = (signal+noise)
  # modeling
  dy = diff(smth.gau(data,gamma),na.rm=T)
  ddy = diff(dy[!is.na(dy)])
  model2 = cpTest(x=ddy,order=2,gamma=gamma,alpha=0.05) # alpha = 0.05
  fdr_pow = Fdr(uh=c(model2$peak),th=h[-length(h)],b=gamma) # b = gamma
  return(fdr_pow)
}
## parallel computing
t0 = Sys.time()
c = 40 # number of cpu cores
cl = parallel::makeForkCluster(c)
doParallel::registerDoParallel(cl)
bc = inverse.snr(Snr = c(seq(3,10,0.5),seq(11,30,1)),2,gamma)
out.linebreak = foreach::foreach(i = 1:length(bc),.packages=c("dSTEM","MASS")) %dopar%{
  l = 15000
  h = seq(150,l,150)
  jump = rep(0,length(h)+1)
  diffb = bc[i] # diffb = {0.15 0.2 0.3}
  beta1 = seq(from=0,by=diffb,length.out = length(h)+1)
  signal = gen.signal(l,h,jump,beta1)
  do.call("rbind",lapply(1:R,D.linebreak,signal=signal))
}
t.linebreak = Sys.time() - t0

#### Type II Change Point (piecewise step function)
gamma = 10
## modeling
R = 2000 # repeatation 
D.stepjump = function(i,signal){
  noise = rnorm(length(signal),0,1) 
  data = (signal+noise)
  # modeling
  dy = diff(smth.gau(data,gamma),na.rm=T)
  ddy = diff(dy[!is.na(dy)])
  model1 = cpTest(x=dy,order=1,alpha=0.05,gamma=gamma,is.constant=T) # alpha = 0.05
  fdr_pow = Fdr(uh=c(model1$peak),th=h[-length(h)],b=gamma) # b = gamma
  return(fdr_pow)
}
## parallel computing
bc = inverse.snr(Snr = c(seq(3,10,0.5),seq(11,30,1)),1,gamma,0)
t0 = Sys.time()
out.stepjump = foreach::foreach(i = 1:length(bc),.packages=c("dSTEM","MASS")) %dopar%{
  l = 15000
  h = seq(150,l,150)
  jump = rep(bc[i],length(h)+1) # jump = {1 1.5 2}
  beta1 = rep(0,length(h)+1)
  signal = gen.signal(l,h,jump,beta1)
  do.call("rbind",lapply(1:R,D.stepjump,signal=signal))
}
t.stepjump = Sys.time() - t0

#### Type II Change Point (piecewise linear function)
gamma = 10
## modeling
R = 2000 # repeatation 
D.linejump = function(i,signal){
  noise = rnorm(length(signal),0,1) 
  data = (signal+noise)
  dy = diff(smth.gau(data,gamma),na.rm=T)
  ddy = diff(dy[!is.na(dy)])
  model2 = cpTest(x=ddy,order=2,gamma=gamma,alpha=0.1)
  breaks = est.pair(model2$vall,model2$peak,gamma)$cp
  slope = est.slope(data,breaks)
  model1 = cpTest(x=dy,order=1,alpha=0.05,gamma=gamma,breaks=breaks,slope=slope)
  fdr_pow = Fdr(uh=c(model1$peak),th=h[-length(h)],b=gamma) # b = gamma
  return(fdr_pow)
}
## parallel computing
bc = inverse.snr(Snr = c(seq(3,10,0.5),seq(11,30,1)),1,gamma,0.05)
t0 = Sys.time()
out.linejump = foreach::foreach(i = 1:length(bc),.packages=c("dSTEM","MASS")) %dopar%{
  l = 15000
  h = seq(150,l,150)
  jump = rep(bc[i],length(h)+1) # jump = {2 2.5 3}
  addb = 0.05  # k_j + k_{j+1} = {0.02 0.05 0.1}  
  beta1 = c(rep(c(0,addb),floor(length(h)/2)),0) 
  signal = gen.signal(l,h,jump,beta1)
  do.call("rbind",lapply(1:R,D.linejump,signal=signal))
}
t.linejump = Sys.time() - t0

#### Mixture of Type I and Type II Change Point  (ratio = 0.5)
gamma = 10
## modeling
R = 2000 # repeatation 
D.mixture = function(i,signal){
  noise = rnorm(length(signal),0,1) # white noise
  data = (signal+noise)
  dy = diff(smth.gau(data,gamma),na.rm=T)
  ddy = diff(dy[!is.na(dy)])
  model2 = cpTest(x=ddy,order=2,gamma=gamma,alpha=0.1)
  breaks = est.pair(model2$vall,model2$peak,gamma)$cp
  if(length(breaks)==0) breaks = floor(length(data)/2)
  slope = est.slope(data,breaks)
  model1 = cpTest(x=dy,order=1,alpha=0.05,gamma=gamma,breaks=breaks,slope=slope)
  fdr_pow_jump = Fdr(uh=c(model1$peak,model1$vall),th=h[(nline+1):(length(h)-1)],b=gamma) # b = gamma
  jump_seg = unique(c(sapply(c(model1$peak,model1$vall),function(x) (x-2*gamma):(x+2*gamma))))
  model2_new = cpTest(x=ddy,order=2,gamma=gamma,alpha=0.1,untest=jump_seg)
  fdr_pow_break = Fdr(uh=model2_new$vall,th=h[1:nline],b=gamma) # b = gamma
  return(cbind(fdr_pow_jump,fdr_pow_break))
}
## parallel computing
bc = inverse.snr(Snr = c(seq(3,10,0.5),seq(11,30,1)),2,gamma)
ac = inverse.snr(Snr = c(seq(3,10,0.5),seq(11,30,1)),1,gamma,0.05)
t0 = Sys.time()
out.mixture = foreach::foreach(i = 1:length(bc),.packages=c("dSTEM","MASS"),.errorhandling="pass") %dopar%{
  ratio = 0.5 #{0.3 0.5 0.7} 
  l = 15000*2
  h = seq(150,l,150)
  njump = floor(ratio*length(h))
  nbreak = length(h) - njump
  jump = c(rep(0,nbreak+1),rep(ac[i],njump))
  diffb = -bc[i]# dk = {0.15 0.2 0.3}
  addb = 0.05    # addb = 0.05
  beta1 = c(seq(from=0,by=diffb,length = nbreak+1),rep(c(0,addb),floor(njump/2)))
  if(length(beta1)-length(h) < 1) beta1 = c(beta1,0)
  signal = gen.signal(l,h,jump,beta1)
  do.call("rbind",lapply(1:R,D.mixture,signal=signal))
}
t.mixture = Sys.time() - t0

