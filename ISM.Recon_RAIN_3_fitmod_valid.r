rm(list=ls())
setwd("/Users/emilygill/Documents/University of Colorado/PHD Research/6. SST Reconstruction/REDO Sep2015/")



# 1. Pull in the climatology for the full slab and SSTs for the full slab.
# 2. Read in both proxy data and organize and smooth
# 3. Using coordinates from #2, create a matrix of contemporary SST anomalies only where we have paleodata.
# 4. Scale paleodata using contemporary SSTs from that point
# 5. Pull in timeseries of total monsoon rainfall for each region of India.
# 6. Fit monsoon total rainfall for that particular region as a function of SSTs at all locations.
# 	Use something that chooses the best predictors.
# 7. Predict rainfall using paleo data.



###################################################################################################
##### 1. ## SOURCE FUNCTIONS AND READ IN LIBRARIES
source("BINxygrid.R")
source("functions_3.R")
source("read.recon.uk3.R")
source("read.recon.mgca3.R")
#source("initialize.R")
library(akima)			# for quilt.plot
library(MASS)			# for stepAIC
library(RColorBrewer)
library(locfit)			# for local polynomials
library(fields)			# for image.plot
myPalette1 <- colorRampPalette(rev(brewer.pal(9, "RdBu")), space="Lab")
myPalette2 <- colorRampPalette(rev(brewer.pal(9, "Greens")), space="Lab")
myPalette3 <- colorRampPalette(rev(brewer.pal(9, "Blues")), space="Lab")
myPalette4 <- colorRampPalette(rev(brewer.pal(9, "YlOrRd")), space="Lab")
myPalette5 <- colorRampPalette(rev(brewer.pal(9, "YlGnBu")), space="Lab")


proxy = "both"							# MgCa or Uk
region = "EqPac"						# This will always be full slab (I think)?? Maybe not.
season = "annual"

train1 = 1980
train2 = 2004
valid1 = 1901
valid2 = 1979

###################################################################################################
##### 2. ## FULL SLAB SST DATA AND CLIMATOLOGY. SCALE FULLSLAB DATA.
results = read.contemp(region = region, season = season) 	   
data.un <- results$data
xgrid = results$xgrid
nx = length(xgrid)
ygrid = results$ygrid
ny = length(ygrid)
grid <- results$grid
index1 <- results$ind
nsta <- results$nsta
clim <- results$clim
climdata = clim[index1]
data = sweep(data.un,2,climdata)




###################################################################################################
##### 3. ## READ IN PROXY DATA, ORGANIZE and SMOOTH.
if (proxy == "Uk"){
	sst.raw = read.recon.uk3(plot = TRUE)	
} else if (proxy == "MgCa"){
	sst.raw = read.recon.mgca3(plot = TRUE)
} else if (proxy == "both"){
	sst.raw.uk = read.recon.uk3(plot = TRUE)
	sst.raw.mgca = read.recon.mgca3(plot = TRUE)

	EP = c(sst.raw.uk$EPlist, sst.raw.mgca$EPlist)
	WP = c(sst.raw.uk$WPlist, sst.raw.mgca$WPlist)
	IO = c(sst.raw.uk$IOlist, sst.raw.mgca$IOlist)
	coord.EP = rbind(sst.raw.uk$EPcoord, sst.raw.mgca$EPcoord)
	coord.WP = rbind(sst.raw.uk$WPcoord, sst.raw.mgca$WPcoord)
	coord.IO = rbind(sst.raw.uk$IOcoord, sst.raw.mgca$IOcoord)

	sst.raw = list(EPlist=EP, WPlist=WP, IOlist=IO, 
		EPcoord=coord.EP, WPcoord=coord.WP, IOcoord=coord.IO)

	## Add something in that only selects unique locations??

}




###################################################################################################
##### 4. ## SMOOTH AND PLOT PALEO RECONSTRUCTIONS
sst.smooth = gcv.recon(sst.raw$EPlist, sst.raw$WPlist, sst.raw$IOlist, 
	EPdeg=2, WPdeg=2, IOdeg=2, base=TRUE)
#pEP = plot.recon(sst.raw$EPlist,sst.smooth$EP, map=TRUE, region="EP",sst.raw$EPcoord)
#pWP = plot.recon(sst.raw$WPlist,sst.smooth$WP, map=TRUE, region="WP",sst.raw$WPcoord)
#pIO = plot.recon(sst.raw$IOlist,sst.smooth$IO, map=TRUE, region="IO",sst.raw$IOcoord)




###################################################################################################
##### 5. ## CREATE A MATRIX OF CONTEMPORARY SSTS AT POINTS WHERE WE HAVE PALEO DATA
if (region == "EqPac"){
	reconpt <- rbind(sst.raw$EPcoord,sst.raw$WPcoord)
	if (proxy == "both"){
		doublemgca = c(14,15,17,18,28,30)		# when double, remove mg/ca
		doubleuk = c(2,6,8,10,23,24)			# when double, remove uk
		mixed = c(14,15,17,18,23,24) 			# only alk in east and mg/ca in west
		mixed2 = c(2,6,8,18,23,24) 				# taking the best resolution records
		mixed3 = c(2,6,8,9,18,21,23,24,26,29,30,31) # taking the best resolution records and removing the short records
	} else if (proxy == "Uk"){
		mixed3 = c(9,16,19)
	} else if (proxy == "MgCa"){
		mixed3 = c(6,9,10,11)
	}
	reconpt = reconpt[-mixed3,]	
} 

nreconpt <- dim(reconpt)[1] 
index <- matrix(seq(1:(nx*ny)),nx,ny)
reconsst <- c()
reconsst.a <- c()

for (i in 1:nreconpt){
	lat <- reconpt[i,1]
	xlat <- which(ygrid == lat)

	lon <- reconpt[i,2]
	xlon <- which(xgrid == lon)

	zz <- index[xlon,xlat]
	zsst <- data[,which(index1 == zz)]
	zclim <- climdata[which(index1 == zz)]

	print(which(index1 == zz))

	reconsst <- cbind(reconsst,zsst)				# (160 x nreconpt) Pulls all data at the point
	reconsst.a <- cbind(reconsst.a,zclim)			# (1 x nreconpt) Using the set climatology
}






###################################################################################################
##### 7. ## PULL IN THE RAINFALL DATA FOR VARIOUS REGIONS and calculate avg and total at each grid point
xgrid <- seq(66.5,100.5,by=1)
nx <- length(xgrid) # nrows = 35
ygrid <- seq(6.5,38.5,by=1)
ny <- length(ygrid) # ncols = 33

ntime = 37986
raindata <- readBin("datafiles/1901-2004_IMD_RAIN_daily.r4",
		what="numeric",n=(nx*ny*ntime),size=4,endian="swap")

rain <- BINxygrid(raindata,xgrid,ygrid,ntime,type="special", undefined=-999, condition="greater")
precip <- rain$data
index = rain$ind
leap = seq(1124,ntime,1460)	# leap years starting with Feb 29 1904
precip = precip[-leap,]		# 104 years of data

# Create average daily monsoon rainfall series, and total monsoon series.
ntime = dim(precip)[1]
Jan1s = seq(1,ntime,365)
Jun1s = seq(152,ntime,365)
Sep30s = seq(273,ntime,365)
Dec31s = seq(365,ntime,365)

precip.davg = c()
precip.sestot = c()
precip.sestot.ann = c()
for (i in 1:length(Jun1s)){
	xavg = apply(precip[Jun1s[i]:Sep30s[i],],2,mean)
	xtot = apply(precip[Jun1s[i]:Sep30s[i],],2,sum)
	xtota = apply(precip[Jan1s[i]:Dec31s[i],],2,sum)

	precip.davg = rbind(precip.davg,xavg)
	precip.sestot = rbind(precip.sestot,xtot)
	precip.sestot.ann = rbind(precip.sestot.ann,xtota)
}

# avgann.perc = (apply(precip.sestot,2,mean))*100/(apply(precip.sestot.ann,2,mean))

# zfull = rep(NaN,(nx*ny))
# zfull[index] = avgann.perc
# zmat = matrix(zfull,nrow=nx,ncol=ny)
# image.plot(xgrid,ygrid,zmat,
# 	ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),99), axes=FALSE, zlim=c(0,100), 
# 	ylab="",xlab="",col=rev(myPalette5(100)), horizontal=TRUE, cex=1.5
# 	)
# axis(1,cex.axis=1.5)
# axis(2,cex.axis=1.5)
# mapnames <- map("world2", 
# 	xlim=c(min(ygrid), max(xgrid)), ylim=c(min(ygrid), max(ygrid)), 
# 	boundary=TRUE, interior=TRUE, col="grey50", add=TRUE
# )
# contour(xgrid, ygrid, zmat, add = TRUE, levels = seq(0,100,10), lwd = 1.5, col="white", labcex=1.25)





###################################################################################################
##### 6. ## SCALE ALL THE PALEO DATA (A COUPLE WAYS) and PICK YEAR TO RECONSTRUCT
if (region == "EqPac"){
	smoodat = c(sst.smooth$EP, sst.smooth$WP)
	if (proxy == "both"){
		smoodat = smoodat[-mixed3]
	}
} 

# sst.scale = scale.recon(smoodat, nreconpt, reconsst.mu, reconsst.sig)
sst.scale = scale.recon(smoodat, nreconpt, reconsst.a)
mat.unsc = sst.scale$unscaled			# Not scaled, just put into a matrix
# mat.sc = sst.scale$scaled				# Put into matrix and scaled using itself
mat.anom = sst.scale$anom				# Put into matrix and scaled using the climatology

mat.base = sst.scale$base
mat.baselim = sst.scale$baselim
# Pick which scaled series you want to use
mat = mat.baselim	


yrsBP = 0:10

###################################################################################################
##### 7. ## PCA ON RAINFALL (seasonal total)
# means  = apply(precip.sestot,2,mean)
# a.precip = sweep(precip.sestot,2,means)
# sc.precip = a.precip

sc.precip = scale(precip.davg)
contemp.sd = attr(sc.precip,'scaled:scale')
contemp.mu = attr(sc.precip,'scaled:center')

aa = which(1901:2004 == train1)
bb = which(1901:2004 == train2)
xdata = sc.precip[aa:bb,]

zs.rain <- var(xdata)
zsvd.rain <- svd(zs.rain)
pcs.rain <- t(t(zsvd.rain$u) %*% t(xdata))
lambdas.rain <- (zsvd.rain$d/sum(zsvd.rain$d))


# zpred.davg = matrix(NA,11,6)
# zpred.davg.se = matrix(NA,11,6)
# zpred.sestot = matrix(NA,11,6)
# zpred.sestot.se = matrix(NA,11,6)
# for (i in 1:length(yrsBP)){
yrBP = 10

missing = which(is.na(mat[yrBP+1,]))
if (length(missing) > 0){
	var = seq(1,nreconpt,1)[-missing]
} else {
	var = seq(1,nreconpt,1)
}

istart = which(1854:2013 == 1901)
istop = which(1854:2013 == 2004)
zreconsst = reconsst[istart:istop,var]


xreconsst = zreconsst[aa:bb,]
zs.p <- var(xreconsst)
zsvd.p <- svd(zs.p)
pcs.p <- t(t(zsvd.p$u) %*% t(xreconsst))
lambdas.p <- (zsvd.p$d/sum(zsvd.p$d))


###################################################################################################
##### 9. ## CCA TO PREDICT RAINFALL SPATIALLY
###################################################################################################
##### 7. ## PERFORM CCA to reconstruct U-WIND

cc = which(1901:2004 == valid1)
dd = which(1901:2004 == valid2)
valid = reconsst[cc:dd,]
valid.rain = sc.precip[cc:dd,]

zs.valid <- var(valid.rain)
zsvd.valid <- svd(zs.valid)
pcs.valid <- t(t(zsvd.valid$u) %*% t(valid.rain))
lambdas.valid <- (zsvd.valid$d/sum(zsvd.valid$d))


ncca = 4
Y = pcs.rain[,1:ncca]
X = pcs.p[,1:ncca]
N = length(Y[,1])

X1 = scale(X)
Xmean = attr(X1,'scaled:center')
Xsd = attr(X1,'scaled:scale')
Y1 = scale(Y)
Ymean = attr(Y1,'scaled:center')
Ysd = attr(Y1,'scaled:scale')

Qx1 = qr.Q(qr(X1))
Qy1 = qr.Q(qr(Y1))
T11 = qr.R(qr(X1))
T22 = qr.R(qr(Y1))

VV=t(Qx1) %*% Qy1
BB = svd(VV)$v
AA = svd(VV)$u 

BB = solve(T22) %*% svd(VV)$v * sqrt(N-1)
wm1=Y1 %*% BB

AA = solve(T11) %*% svd(VV)$u * sqrt(N-1)
vm1 = X1 %*% AA

cancorln = svd(VV)$d[1:ncca]   #canonical correlation
 
Fyy = var(Y1) %*% BB

betahat = solve(t(AA) %*% t(X1)%*% X1 %*% AA) %*% t(AA) %*% t(X1) %*% Y1

### Pull in paleo data and put in PC space.
xx.hol = valid
newdata = xx.hol %*% zsvd.p$u[,1:ncca]

### Predict past wind PC field.
pcmeans = matrix(colMeans(pcs.rain[,-1:-ncca]),nrow=1)
fpc.val = matrix(NA,length(valid1:valid2),ncol(pcs.rain))
scaled = newdata %*% AA %*% betahat
unscaled = t(t(scaled) * Ysd + Ymean)
fpc.val[,1:ncca] = unscaled
fpc.val[,-1:-ncca] = pcmeans
recon.val = fpc.val %*% t(zsvd.rain$u)
recon.val.unsc = t(t(recon.val) * contemp.sd + contemp.mu)


# ####################### CALIB STATISTICS #######################
### Statistics on just the signal data
nkeep = ncca
Ev = matrix(0, nrow = dim(pcs.rain)[2], ncol = dim(pcs.rain)[2])
Ev[,1:nkeep] = zsvd.valid$u[,1:nkeep] 

rain.sig = pcs.valid %*% t(Ev)
rain.sig.unsc = t(t(rain.sig) * contemp.sd + contemp.mu)

# Do stats on rain.sig.unsc vs. recon.val.unsc

R2_rain = c()
for (k in 1:length(index)){
	R2_rain = c(R2_rain, cor(recon.val.unsc[,k], rain.sig.unsc[,k])^2)
}

sumdiff = apply((rain.sig.unsc - recon.val.unsc)^2,2,sum)
sumYraw = apply(rain.sig.unsc^2,2,sum)
beta_rain = 1 - (sumdiff/sumYraw)

write(R2_rain, file="./postfiles/ISM/diagnostics/valid_r2_unscsig_rain.txt")
write(beta_rain, file="./postfiles/ISM/diagnostics/valid_b_unscsig_rain.txt")



##################################################
dev.new(width=7,height=7)
zfull = rep(NaN,(nx*ny))
zfull[index] = beta_rain
zmat = matrix(zfull,nrow=nx,ncol=ny)
image(xgrid,ygrid,zmat,
	ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)),
	xlab="",ylab="", yaxt="n", xaxt="n", main=expression(paste("Rainfall Verification ", beta)),
	cex.main=2, col=rev(myPalette4(100)), zlim=c(-1,1)
	)
axis(1,seq(30,100,10), cex.axis=2)
axis(2,cex.axis=2)
mapnames <- map("world2", 
	xlim=c(min(ygrid), 98), ylim=c(min(ygrid), max(ygrid)), 
	boundary=TRUE, interior=TRUE, col="black", add=TRUE
)
#contour(xgrid, ygrid, zmat, add = TRUE, levels = c(0), lwd = 2, labcex=1.25, col="black")

dev.new(width=7,height=7)
zfull = rep(NaN,(nx*ny))
zfull[index] = R2_rain
zmat = matrix(zfull,nrow=nx,ncol=ny)
image(xgrid,ygrid,zmat,
	ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)),
	xlab="",ylab="", yaxt="n", xaxt="n", main=expression(paste("Rainfall Verification ", R^{2})),
	cex.main=2, col=rev(myPalette2(100)), zlim=c(0,1)
	)
axis(1,seq(30,100,10), cex.axis=2)
axis(2,cex.axis=2)
mapnames <- map("world2", 
	xlim=c(min(ygrid), 98), ylim=c(min(ygrid), max(ygrid)), 
	boundary=TRUE, interior=TRUE, col="black", add=TRUE
)
#contour(xgrid, ygrid, zmat, add = TRUE, levels = seq(0,1,0.25), lwd = 2, labcex=1.25, col="black")







#### PLOT COMPARISON YEARS
iix = which(1901:2004 == 1988)
year = "1988-1989"
max = 4.5
min = -4.5

dev.new(height=10, width=3)
par(mfrow = c(2, 1), mar = c(2,2,2,2))
# Actual
zfull = rep(NaN,(nx*ny))
zfull[index] = sc.precip[iix,]*contemp.sd					# add in [index1] if youre plotting India rainfall or SST
zmat = matrix(zfull,nrow=nx,ncol=ny)
image(xgrid,ygrid,zmat,
	ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)),
	xlab="",ylab="", main=paste("Actual", year),
	cex.axis=1.1, yaxt="n", col=rev(myPalette1(100)), zlim=c(-11,11)
	)
axis(2, at = c(10,20,30))
mapnames <- map("world2", 
	xlim=c(min(ygrid), max(xgrid)), ylim=c(min(ygrid), max(ygrid)), 
	boundary=TRUE, interior=TRUE, add=TRUE, lwd=1.5, col="grey50"
)

# Reconstructed
zfull = rep(NaN,(nx*ny))
zfull[index] = recon.contemp[iix,]*contemp.sd  						# add in [index1] if youre plotting India rainfall or SST
zmat = matrix(zfull,nrow=nx,ncol=ny)
image.plot(xgrid,ygrid,zmat,
	ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)),
	xlab="",ylab="", main=paste("Reconstructed", year),
	cex.axis=1.1, yaxt="n", col=rev(myPalette1(100)), zlim=c(-11,11)
	)
axis(2, at = c(10,20,30))
mapnames <- map("world2", 
	xlim=c(min(ygrid), max(xgrid)), ylim=c(min(ygrid), max(ygrid)), 
	boundary=TRUE, interior=TRUE, add=TRUE, lwd=1.5, col="grey50"
)




