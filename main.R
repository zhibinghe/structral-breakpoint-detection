######################################################
############# Main File #############################
######################################################
set.seed(202100)
# devtools::install_github("zhibinghe/dSTEM")
library(foreach)
library(doParallel)
library(tictoc)
library(dSTEM)
library(not)
library(strucchange)
library(nsp)
#### signal simulation settings 
signal.info = function(type = c("I","II-step","II-linear","mixture")) {
  type = match.arg(type)
  # changes in slope
  if(type == "I") {
    name = "I"
    cpt.type = "pcwsLinContMean"
    n = 1200
    cpt = seq(150,by=150,length.out=6)
    jump.size = rep(0,7)
    beta1 = c(2,-1,2.5,-3,-0.2,2.5)/25
    # To make the endpoint be the same with start-point
    slope = c(beta1,-sum(beta1*(c(cpt[1],diff(cpt))))/(n-tail(cpt,1)))
  }
  # changes in intercept
  else if(type == "II-step") {
    name = "II-step"
    cpt.type = "pcwsConstMean"
    n = 1200
    cpt = seq(150,by=150,length.out=6) 
    jump.size = c(0,1.5,2,2.2,1.8,2,1.5)*2
    slope = rep(0,length(cpt)+1)
  }
  # changes in intercept and slope
  else if(type == "II-linear") {
    name = "II-linear"
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

#### Method comparation
## NOT, NSP, BP methods 
comp.detection = function(data,method = c("dstem","not","nsp","bp"),
                          type = c("I","II-step","II-linear","mixture"),
                          level=0.05,gamma=20,M=NULL) {
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
    out = round(rowMeans(nsp_poly(data,M=M,sigma=1,alpha=level,
                                  deg=degree)$intervals[,-3]))
  }
  #
  if (method == "bp") {
    out = breakpoints(data~seq(1,length(data)),h=M)$breakpoints
  }
  return(out)
}

