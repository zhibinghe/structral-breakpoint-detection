########################################################
################ FDR/Power vs. gamma  ##################
########################################################
source("main.R")
########################################################
######## Type I Change Point
Gamma = 3:25
R = 100 # repeatation
## modeling
D.type1 = function(signal,th,gamma,R,b = 10){
  f = function(i){
    noise = rnorm(length(signal),0,1) # white noise
    model = dstem(signal + noise, "I", gamma = gamma, alpha = 0.05)
    Fdr(uh=c(model$peak,model$vall),th=th,b=b) }
  colMeans(do.call(rbind,lapply(1:R,f)))
}
## parallel computing
n.cpu = 4 # number of cpu cores
cl = parallel::makeForkCluster(n.cpu)
doParallel::registerDoParallel(cl)
out.gamma = foreach::foreach(i = 1:length(Gamma),.combine = "rbind", .packages=c("dSTEM","MASS"),
                             .errorhandling = "pass") %dopar%{
  l = 15000
  h = seq(150,l,150)
  jump = rep(0,length(h)+1)
  diffb = 0.15 # dk = {0.15 0.2 0.3}
  beta1 = seq(from=0,by=diffb,length = length(h)+1)
  signal = gen.signal(l,h,jump,beta1)
  D.type1(signal,th=h[-length(h)],gamma=Gamma[i],R = R)
}
#########################################################
########## Type II-step Change Point 
D.type2step = function(signal,th,gamma,R,b = 10){
  f = function(i){
    noise = rnorm(length(signal),0,1) # white noise
    model = dstem(signal + noise, "II-step", gamma = gamma, alpha = 0.05)
    Fdr(uh=c(model$peak,model$vall),th=th,b=b) }
  colMeans(do.call(rbind,lapply(1:R,f)))
}
## parallel computing
out.gamma = foreach::foreach(i=1:length(Gamma),.combine = "rbind",.packages=c("dSTEM","MASS"),
                             .errorhandling = "pass") %dopar%{
  l = 15000
  h = seq(150,l,150)
  jump = rep(2,length(h)+1) # jump = {1 1.5 2}
  beta1 = rep(0,length(h)+1)
  signal = gen.signal(l,h,jump,beta1)
  D.type2step(signal,th=h[-length(h)],gamma=Gamma[i],R = R)
 }                            
##########################################################
############ Type II-linear Change Point
D.type2line = function(signal,th,gamma,R,b = 10){
  f = function(i){
    noise = rnorm(length(signal),0,1) # white noise
    model = dstem(signal + noise, "II-linear", gamma = gamma, alpha = 0.05)
    Fdr(uh=c(model$peak,model$vall),th=th,b=b) }
  colMeans(do.call(rbind,lapply(1:R,f)))
}
## parallel computing
out.gamma = foreach::foreach(i=1:length(Gamma),.combine = "rbind",.packages=c("dSTEM","MASS"),
                             .errorhandling = "pass") %dopar%{
  l = 15000
  h = seq(150,l,150)
  jump = rep(2,length(h)+1) # jump = {1 1.5 2}
  beta1 = rep(0,length(h)+1)
  signal = gen.signal(l,h,jump,beta1)
  D.type2linear(signal,th=h[-length(h)],gamma=Gamma[i],R = R)
}   
##########################################################
######## Mixture of Type I and Type II Change Point 
D.type.mix = function(signal,th1,th2,th,gamma,R,b = 10){
  f = function(i){
    noise = rnorm(length(signal),0,1) # white noise
    model = dstem(signal + noise, "mixture", gamma = gamma, alpha = 0.05)
    cbind(Fdr(uh=c(model$type1$peak,model$type1$vall),th=th1,b=b),
          Fdr(uh=c(model$type2$peak,model$type2$vall),th=th2,b=b),
          Fdr(uh=as.vector(unlist(model)),th=th,b=b))
  }
  colMeans(do.call(rbind,lapply(1:R,f)))
}
## parallel computing
out.gamma = foreach::foreach(i=1:length(Gamma),.combine = "rbind",.packages=c("dSTEM","MASS"),
                             .errorhandling = "pass") %dopar%{
  ratio = 0.5 #{0.3 0.5 0.7}
  l = 15000 * 2
  h = seq(150, l, 150)
  njump = floor(ratio * length(h))
  nline = length(h) - njump
  jump = c(rep(0, nline + 1), rep(3, njump))
  diffb = 0.2    # dk = {0.15 0.2 0.3}
  addb = 0.05    # addb = 0.05
  beta1 = c(seq(from = 0,by = diffb,length = nline + 1), rep(c(0, addb), floor(njump / 2)))
  if (length(beta1) - length(h) < 1) beta1 = c(beta1, 0)
  signal = gen.signal(l, h, jump, beta1)
  D.type.mix(signal,th1 =h[1:nline],th2 = h[(nline+1):(length(h)-1)], 
             th=h[-length(h)],gamma=Gamma[i],R = R)
                             }
parallel::stopCluster(cl)

