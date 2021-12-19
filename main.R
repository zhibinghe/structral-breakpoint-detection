####************ source main functions
set.seed(202100)
##
library(dSTEM)
library(not)
library(strucchange)
source("NSP_for_Github_v4.R")
# signal setting 
signal.info = function(type = c("I","II-step","II-linear","mix")) {
  type = match.arg(type)
  # changes in slope
  if(type == "I") {
    name = "Type I"
    cpt.type = "pcwsLinContMean"
    n = 1200
    cpt = seq(150,by=150,length.out=6)
    jump.size = rep(0,7)
    beta1 = c(2,-1,2.5,-3,-0.2,2.5)/25
    # To make the endpoint be the same with start-point
    slope = c(beta1,-sum(beta1*(c(cpt[1],diff(cpt))))/(l-tail(cpt,1)))
  }
  # changes in intercept
  else if(type == "II-step") {
    name = "Type II (stepwise)"
    cpt.type = "pcwsConstMean"
    n = 1200
    cpt = seq(150,by=150,length.out=6) 
    jump.size = c(0,1.5,2,2.2,1.8,2,1.5)*2
    slope = rep(0,length(cpt)+1)
  }
  # changes in intercept and slope
  else if(type == "II-linear") {
    name = "Type II (linear)"
    cpt.type = "pcwsLinMean"
    n = 1200
    cpt = seq(150,by=150,length.out=6) 
    jump.size = c(0,1.5,2,-2.2,1.8,2,-1.5)*3
    slope = c(2,-1,2.5,-3,-0.2,2.5,-0.5)/25
  }
  # mixture of two types
  else {
    name = "mixture"
    cpt.type = "mixture"
    n = 1500
    cpt = seq(150,by=150,length.out=8)
    jump.size = c(0,0,2.5,0,3,0,0,2.8,-2.5)*3
    slope = c(1,-1,1.5,-2,-0.2,1.2,-0.5,1.2,0)/25
  }
  list(name=name,cpt.type=cpt.type,n=n,cpt=cpt,jump.size=jump.size,slope=slope)
}

#### Calculate slope change based on SNR 
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
#### dstem method
## if type is 'mixture', the output is a list of type I and type II change points,
## otherwise, it is a vector indicating the change points
dstem = function(data,type = c("I","II-step","II-linear","mixture"),gamma,level=0.05){
  type = match.arg(type)
  dy = diff(smth.gau(data,gamma))
  ddy = diff(dy)
  if (type == "I") {
    est = cpTest(x=ddy,order=2,gamma=gamma,alpha=level)
    out = list(vall = est$vall,peak = est$peak)
  }
  else if (type == "II-step") {
    est = cpTest(x=dy,order=1,alpha=level,gamma=gamma,is.constant=T)
    out = list(vall = est$vall,peak = est$peak)
  }
  else if (type == "II-linear") {
    model2 = cpTest(x=ddy,order=2,gamma=gamma,alpha=2*level)
    breaks = est.pair(model2$vall,model2$peak,gamma)$cp
    slope = est.slope(data,breaks)
    est = cpTest(x=dy,order=1,alpha=level,gamma=gamma,breaks=breaks,slope=slope)
    out = list(vall = est$vall, peak = est$peak)
  }
  # mixture  
  else {
    model2 = cpTest(x=ddy,order=2,gamma=gamma,alpha=2*level)
    breaks = est.pair(model2$vall,model2$peak,gamma)$cp
    if(length(breaks)==0) breaks = floor(length(data)/2)
    slope = est.slope(data,breaks)
    model1 = cpTest(x=dy,order=1,alpha=level,gamma=gamma,breaks=breaks,slope=slope)
    jump_seg = unique(c(sapply(c(model1$peak,model1$vall),function(x) floor(x-2*gamma):ceiling(x+2*gamma))))
    est = cpTest(x=ddy,order=2,gamma=gamma,alpha=level,untest=jump_seg)
    out = list(type1 = list(vall=est$vall,peak=est$peak), type2 = list(vall=model1$vall,peak=model1$peak))
  }
  return(out)
}
#### Method comparation
comp.detection = function(data,method = c("dstem","not","nsp","bp"),
                          type = c("I","II-step","II-linear","mixture"), level=0.05,gamma=NULL,M=NULL) {
  # M: initial number of change points
  method = match.arg(method)
  type = match.arg(type)
  if (method == "dstem") {
    out = dstem(data,type,gamma,level)
  }
  # 
  if (method != "dstem" & type == "mixture") 
    stop("Cannot deal with mixture signal")
  # 
  if (method == "not") {
    if (type == "I") 
      out = features(not(data,contrast="pcwsLinContMean"),q.max=M)$cpt
    else if (type == "II-step") 
      out = features(not(data,contrast="pcwsConstMean"),q.max=M)$cpt
    else
      out = features(not(data,contrast="pcwsLinMean"),q.max=M)$cpt
  }
  #
  if (method == "nsp") {
    degree = ifelse (type == "II-step", 0, 1)
    out = round(rowMeans(nsp_poly(data,M=M,sigma=1,alpha=0.05,deg=deg)$intervals[,-3]))
  }
  #
  if (method == "bp") {
    out = breakpoints(data~seq(1,length(data)),h=M)$breakpoints
  }
  return(out)
}

