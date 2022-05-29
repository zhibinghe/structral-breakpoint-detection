#####################################################
############## Simulation ###########################
###### Comparason with other methods ################
#####################################################
source("main.R")
################ PART I #############################
############# Short term data sequence ##############
#### Type I
x = signal.info(type="I") # change type to "II-step", "II-linear"
signal = gen.signal(l=x$n,h=x$cpt,jump = x$jump.size,b1 = x$slope)
# x$cpt # true change point locations
# start parallel computing
n.cpu = 4; R = 100 # repetition
cl = parallel::makeCluster(4) # cpu cores
doParallel::registerDoParallel(cl)
# 
cpt = NULL; time = NULL
method = c("dstem","not","nsp")
for (m in method){
  tic()
  if (m == "dstem") {
    method.cpt = foreach(iterators::icount(R),.packages="dSTEM",.errorhandling="pass") %dopar% {
      as.vector(unlist(comp.detection(signal + rnorm(x$n),method="dstem",x$name,gamma=20)))}}
  if (m == "not")
    method.cpt = foreach(iterators::icount(R),.packages="not",.errorhandling="pass") %dopar% {
      comp.detection(signal + rnorm(x$n),method="not",x$name,M=0.1* x$n)}
  if (m == "nsp")
    method.cpt = foreach(iterators::icount(R),.packages="nsp",.errorhandling="pass") %dopar% {
      comp.detection(signal + rnorm(x$n),method="nsp",x$name,M=0.1* x$n)}
  if (m == "bp")
    method.cpt = foreach(iterators::icount(R),.packages="strucchange",.errorhandling="pass") %dopar% {
      comp.detection(signal + rnorm(x$n),method="bp",x$name,M=20)}
  time.run = toc(quiet=TRUE); time.run = time.run$toc - time.run$tic
  cpt = append(cpt,method.cpt)
  time = append(time,time.run)
}
names(time) = method
cpt = split(cpt,as.factor(rep(1:length(method),each=R)))
parallel::stopCluster(cl)
#### Mixture 
#### only dstem can discriminate type I and type II change points
x = signal.info(type="mixture") # change type to "II-step", "II-linear"
signal = gen.signal(l=x$n,h=x$cpt,jump = x$jump.size,b1 = x$slope)
# start parallel computing
n.cpu = 4; R = 100 # repetition
cl = parallel::makeCluster(4) # cpu cores
doParallel::registerDoParallel(cl)
# 
cpt = foreach(iterators::icount(R),.packages="dSTEM",.errorhandling="pass") %dopar% {
      as.vector(unlist(comp.detection(signal + rnorm(x$n),
                                      method="dstem",x$name)))}
 parallel::stopCluster(cl) 
################ PART II #############################
############# Long term data sequence ##############
#### Type I
x = signal.info(type="I")
signal = gen.signal(l=x$n,h=x$cpt,jump = x$jump.size,b1 = x$slope,rep=10)
n = length(signal)
# seq(c(x$cpt,n),n*rep,n) # true change point locations
# start parallel computing
cl = parallel::makeCluster(n.cpu) # cpu cores
doParallel::registerDoParallel(cl)
# 
cpt = NULL; time = NULL
method = c("dstem","not","nsp")
for (m in method){
  tic()
  if (m == "dstem") {
    method.cpt = foreach(iterators::icount(R),.packages="dSTEM",.errorhandling="pass") %dopar% {
      as.vector(unlist(comp.detection(signal + rnorm(n),method="dstem",x$name,gamma=20)))}}
  if (m == "not")
    method.cpt = foreach(iterators::icount(R),.packages="not",.errorhandling="pass") %dopar% {
      comp.detection(signal + rnorm(n),method="not",x$name,M=0.1*n)}
  if (m == "nsp")
    method.cpt = foreach(iterators::icount(R),.packages="nsp",.errorhandling="pass") %dopar% {
      comp.detection(signal + rnorm(n),method="nsp",x$name,M=0.1*n)}
  if (m == "bp")
    method.cpt = foreach(iterators::icount(R),.packages="strucchange",.errorhandling="pass") %dopar% {
      comp.detection(signal + rnorm(n),method="bp",x$name,M=100)}
  time.run = toc(quiet=TRUE); time.run = time.run$toc - time.run$tic
  cpt = append(cpt,method.cpt)
  time = append(time,time.run)
}
names(time) = method
cpt = split(cpt,as.factor(rep(1:length(method),each=R)))
parallel::stopCluster(cl)
###################################################################
######################## Coverage Rate ############################
capture.rate = function(x,th,b,L){
  # x: estimated change points
  # th: true locations
  # b: location tolerance
  # lower,upper: number of b
  seg.capture = function(x,th,b,lower,upper,L){
    f = function(x,lower,upper){
      lowerb = lower*b + 1 
      if(lower*b %% 1 == 0 & lower!=0) c(seq(ceiling(x-upper*b),x-lowerb),seq(x+lowerb,floor(x+upper*b)))
      else c(seq(ceiling(x-upper*b),floor(x-lower*b)),seq(ceiling(x+lower*b),floor(x+upper*b)))
    }
    if(upper == "Inf") seg = setdiff(1:L,unique(unlist(lapply(th,f,lower=0,upper=lower))))
    else seg = unique(unlist(lapply(th,f,lower=lower,upper=upper)))
    round(mean(sapply(x, function(x) sum(x %in% seg)))/length(th),4)
  }
  f = function(t) cbind(seg.capture(t,th=th,b=b,lower=0,upper=1/3,L=L),
                        seg.capture(t,th=th,b=b,lower=1/3,upper=1,L=L),
                        seg.capture(t,th=th,b=b,lower=1,upper=2,L=L),
                        seg.capture(t,th=th,b=b,lower=2,upper=4,L=L),
                        seg.capture(t,th=th,b=b,lower=4,upper=Inf,L=L))
  table = as.data.frame(f(x)); colnames(table) = c("0--1/3","1/3--1","1--2","2--4",">4")
  return(table)
}
# coverage rate table
b = 20 # b = gamma
do.call(rbind,lapply(cpt,capture.rate,th = x$cpt,b=20,L=x$n))
