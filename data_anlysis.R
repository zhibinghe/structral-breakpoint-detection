source("main.R")
#path = getwd()
#setwd(paste0(path,"/real_data"))

#### HST-stock price
hst = read.csv("HST_stock.csv",header=T)
hst$Date <- as.Date(hst$Date, "%m/%d/%Y") 
hsts = subset(hst,Date >= "2018-01-01") # start with 2018-1-1
gamma = 10
hst_dstem = dstem(hsts$Close,"mixture",gamma)
##
plot(Close ~ Date, data = hsts,type="l",ylab="Close Price",xlab="")
points(hsts[unlist(hst_dstem$type1),c("Date","Close")],pch=24,bg="blue")
points(hsts[unlist(hst_dstem$type2),c("Date","Close")],pch=25,bg="blue")
##
hst_not = features(not(hsts$Close,contrast="pcwsLinMean"),q.max=20)$cpt
points(hsts[hst_not,c("Date","Close")],pch=21,bg="red")
##
hst_nsp = round(rowMeans(nsp_poly(hsts$Close,sigma=1,alpha=0.05,deg=1)$intervals[,-3]))
points(hsts[hst_nsp,c("Date","Close")],pch=22,bg="green")
legend("topright",legend=c("dSTEM (type I)","dSTEM (type II)","NOT","NSP"),pch=c(24,25,21,22),
       pt.bg=c("blue","blue","red","green"),bty="n")

#### Covid-19-associated deaths in UK
cvduk <- read.csv("data_2020-Jul-23.csv")
# new Death per day
death = rev(cvduk$newDeathsByPublishDate)
# Anscombe transform 
death = stats::filter(death, rep(1, 7)/7, sides=1) 
death = 2 * sqrt(death[-(1:6)] + 3/8)
ts.plot(death, xlab="Time (days starting from March 12th 2020)", ylab = "",main="Covid19 Death in UK")
## dSTEM
gamma = 6
cvduk_dstem = dstem(death,"mixture",gamma)
points(unlist(cvduk_dstem$type1),death[unlist(cvduk_dstem$type1)],pch=24,bg="blue")
points(unlist(cvduk_dstem$type2),death[unlist(cvduk_dstem$type2)],pch=25,bg="blue")
## NOT
cvduk_not = features(not(death,contrast="pcwsLinMean"),q.max=20)$cpt
points(cvduk_not,death[cvduk_not],pch=21,bg="red")
## NSP
cvduk_nsp = round(rowMeans(nsp_poly(death, deg=1)$intervals[,-3]))
points(cvduk_nsp,death[cvduk_nsp],pch=22,bg="green")
legend("topright",legend=c("dSTEM (type I)","dSTEM (type II)","NOT","NSP"),pch=c(24,25,21,22),
       pt.bg=c("blue","blue","red","green"),bty="n")

#### covid19 data in USA
cvdus = read.csv("covid_usa.csv",header=T)
death = stats::filter(diff(cvdus$deaths), rep(1, 7)/7, sides=1) 
death = 2 * sqrt(death[-(1:6)] + 3/8)
ts.plot(death, xlab="Time (days starting from March 26th 2020)", ylab = "",main="Covid19 Death in USA")
## dSTEM
gamma = 6
cvdus_dstem = dstem(death,"mixture",gamma)
points(unlist(cvdus_dstem$type1),death[unlist(cvdus_dstem$type1)],pch=24,bg="blue")
points(unlist(cvdus_dstem$type2),death[unlist(cvdus_dstem$type2)],pch=25,bg="blue")
## NOT
cvdus_not = features(not(death,contrast="pcwsLinMean"),q.max=20)$cpt
points(cvdus_not,death[cvdus_not],pch=21,bg="red")
## NSP
cvdus_nsp = round(rowMeans(nsp_poly(death, deg=1)$intervals[,-3]))
points(cvdus_nsp,death[cvdus_nsp],pch=22,bg="green")
legend("bottomright",legend=c("dSTEM (type I)","dSTEM (type II)","NOT","NSP"),pch=c(24,25,21,22),
       pt.bg=c("blue","blue","red","green"),bty="n")
