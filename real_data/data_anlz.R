####*************************************
#### HST-stock price
hst = read.csv("HST_stock.csv",header=T)
hst$Date <- as.Date(hst$Date, "%m/%d/%Y") 
hsts = subset(hst,Date >= "2018-01-01") # start with 2018-1-1
gamma=10
model1 = dstem(hsts$Close,"II-linear",gamma)
model2 = dstem(hsts$Close,"mixture",gamma)
plot(Close ~ Date, data = hsts,type="l") 
points(x=hsts[model1$peak,"Date"],y=hsts[model1$peak,"Close"],pch=16,col="green")
points(x=hsts[model1$vall,"Date"],y=hsts[model1$vall,"Close"],pch=19,col="red")


###################################################
#### Covid-19-associated deaths in UK
cvd <- read.csv("data_2020-Jul-23.csv")
# new Death per day
z = rev(cvd$newDeathsByPublishDate)
# Anscombe transform 
zz = stats::filter(z, rep(1, 7)/7, sides=1) 
zz.ans = 2 * sqrt(zz[8:140] + 3/8)
# Modeling
gamma=5
dy = diff(smth.gaussian(zz.ans,window=10*gamma,alpha=5,tails=T),na.rm=T)
ddy = diff(dy[!is.na(dy)])  
model2 = cpTest(x=ddy,order=2,gamma=gamma,alpha=0.05)
# NSP
temp = nsp_poly(zz.ans, deg=1) 
ts.plot(zz.ans, xlab="Time (days starting from March 12th 2020)", ylab = "",main="Covid19 Death in UK")
draw.rects(temp$intervals, range(zz.ans), 20, "grey")
draw.rects(f.int(c(model2$peak,model2$vall),gamma),range(zz.ans),col="green",density=20)
nps_loc = round(rowMeans(temp$intervals[,-3]))
stem_loc = c(model2$peak,model2$vall)
points(nps_loc,zz.ans[nps_loc],pch=15,col="red")
points(stem_loc,zz.ans[stem_loc],pch=15,col="green")

# strucchange
x = 1:length(zz.ans)
bai = strucchange::breakpoints(zz.ans~x,h=10)
bai_loc = bai$breakpoints
points(bai_loc,zz.ans[bai_loc]+0.5,pch=15,col="black")
# segmented
temp=segmented(lm(zz.ans~x),seg.Z=~x, npsi=5, control=seg.control(display=FALSE))
seg_loc = as.vector(round(temp$psi[,2]))
points(seg_loc,zz.ans[seg_loc]+0.8,pch=15,col="blue")
#### covid19 data in USA
dat = covid19("USA",level=1,start="2020-03-20",end="2020-12-31")
#dat = read.csv("covid_usa.csv")
death = diff(dat$deaths)
zz = stats::filter(death, rep(1, 7)/7, sides=1) 
zz.ans = 2 * sqrt(zz[8:140] + 3/8)
# Modeling
gamma=6
dy = diff(smth.gaussian(zz.ans,window=10*gamma,alpha=5,tails=T),na.rm=T)
ddy = diff(dy[!is.na(dy)])  
model2 = cpTest(x=ddy,order=2,gamma=gamma,alpha=0.05)
# NSP
temp = nsp_poly(zz.ans, deg=1) 
nps_loc = round(rowMeans(temp$intervals[,-3]))
ts.plot(zz.ans, xlab="Time (days starting from March 21th 2020)", ylab = "",main="Covid19 Death in USA")
#draw.rects(temp$intervals, range(zz.ans), 20, "grey")
#draw.rects(f.int(c(model2$peak,model2$vall),gamma),range(zz.ans),col="green",density=20)
stem_loc = c(model2$peak,model2$vall)
points(x=stem_loc,y=zz.ans[stem_loc],pch=15,col="green")
points(nps_loc,zz.ans[nps_loc],pch=15,col="red")
#### 
x = 1:length(zz.ans)
bai = strucchange::breakpoints(zz.ans~x,h=10)
bai_loc = bai$breakpoints
points(bai_loc,zz.ans[bai_loc]+0.5,pch=15,col="black")

## segmented
temp=segmented(lm(zz.ans~x),seg.Z=~x, npsi=5, control=seg.control(display=FALSE))
seg_loc = as.vector(round(temp$psi[,2]))
points(seg_loc,zz.ans[seg_loc]+0.8,pch=15,col="blue")
