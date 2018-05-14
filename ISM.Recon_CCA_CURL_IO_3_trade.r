rm(list=ls())
setwd("/Users/emilygill/Documents/University of Colorado/PHD Research/6. SST Reconstruction/REDO Sep2015/")

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
library(CCA) 			# to perform canonical correlation analysis
library(seqinr)			# to shade colors and add transparency
myPalette1 <- colorRampPalette(rev(brewer.pal(9, "RdBu")), space="Lab")
myPalette2 <- colorRampPalette(rev(brewer.pal(9, "Greens")), space="Lab")
myPalette3 <- colorRampPalette(rev(brewer.pal(9, "YlOrRd")), space="Lab")
myPalette4 <- colorRampPalette(rev(brewer.pal(9, "BrBG")), space="Lab")
myPalette5 <- colorRampPalette(rev(brewer.pal(9, "PRGn")), space="Lab")
myPalette6 <- colorRampPalette(rev(brewer.pal(9, "Purples")), space="Lab")


###################################################################################################
##### *** PICK AN AREA TO DO THE RECONSTRUCTION OVER
region = "EqPac"						# EqPac or IO
proxy = "both"							# MgCa or Uk


###################################################################################################
##### 2. ## READ IN CONTEMPORARY SST DATA FOR FULLSLAB
results = read.contemp(region = "EqPac", season = "annual") 	   
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

# Data for overlapping wind period of 1949-2013
data.sst = data[96:160,]


##### 3. ## READ IN CONTEMPORARY SST DATA FOR INDIAN OCEAN ONLY
IOresults = read.contemp(region = "IO", season = "summer") 	   
IOdata.un <- IOresults$data
IOxgrid = IOresults$xgrid
IOnx = length(IOxgrid)
IOygrid = IOresults$ygrid
IOny = length(IOygrid)
IOgrid <- IOresults$grid
IOindex1 <- IOresults$ind
IOnsta <- IOresults$nsta
IOclim <- IOresults$clim
IOclimdata = IOclim[IOindex1]
IOdata = sweep(IOdata.un,2,IOclimdata)

IOdata.sst = IOdata[96:160,]

###################################################################################################
##### 5. ## READ IN PALEO DATA POINT RECONSTRUCTIONS RECONSTRUCTIONS
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

}

nEP = length(sst.raw$EPlist)
EP.smooth = list()
for (i in 1:nEP){
	N = dim(sst.raw$EPlist[[i]])[1]
	zfit = locfit(sst.raw$EPlist[[i]]$SST2 ~ sst.raw$EPlist[[i]]$age, deg=2, kern="bisq", scale=TRUE)

	if (round(sst.raw$EPlist[[i]]$age[N]) > 10){
		maxyr = 10
	} else {
		maxyr = round(EP[[i]]$age[N])
	}
	minyr = floor(EP[[i]]$age[1])

	aa = 10
	bb = 0
	xseq = seq(bb,aa,0.5)
	lseq = length(xseq)

	# yestorig = predict(zfit, seq(bb,aa,1))

	# yestbase500 = predict(zfit,seq(bb,aa,0.5))
	# xsmo500 = rep(NaN,21)
	# x.smooth500 = matrix(cbind(seq(bb,aa,0.5), yestbase500, xsmo500), 21, 3)

	yestbase250 = predict(zfit,xseq)
	xsmo250 = rep(NaN,lseq)
	xmin = which(xseq == minyr)
	xmax = which(xseq == maxyr)
	xsmo250[xmin:xmax] = yestbase250[xmin:xmax]
	x.smooth250 = matrix(cbind(xseq, yestbase250, xsmo250), lseq, 3)

	EP.smooth[[i]] = x.smooth250
}

nWP = length(sst.raw$WPlist)
WP.smooth = list()
for (i in 1:nWP){
	N = dim(sst.raw$WPlist[[i]])[1]
	zfit = locfit(sst.raw$WPlist[[i]]$SST2 ~ sst.raw$WPlist[[i]]$age, deg=2, kern="bisq", scale=TRUE)

	if (round(sst.raw$WPlist[[i]]$age[N]) > 10){
		maxyr = 10
	} else {
		maxyr = round(WP[[i]]$age[N])
	}
	minyr = floor(WP[[i]]$age[1])

	aa = 10
	bb = 0
	xseq = seq(bb,aa,0.5)

	yestbase250 = predict(zfit,xseq)
	xsmo250 = rep(NaN,lseq)
	xmin = which(xseq == minyr)
	xmax = which(xseq == maxyr)
	xsmo250[xmin:xmax] = yestbase250[xmin:xmax]
	x.smooth250 = matrix(cbind(xseq, yestbase250, xsmo250), lseq, 3)

	WP.smooth[[i]] = x.smooth250
	}



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
zzpt <- c()
for (i in 1:nreconpt){
	lat <- reconpt[i,1]
	xlat <- which(ygrid == lat)

	lon <- reconpt[i,2]
	xlon <- which(xgrid == lon)

	zz <- index[xlon,xlat]
	zzpt <- c(zzpt,which(index1 == zz))
	zsst <- data.sst[,which(index1 == zz)]
	zclim <- climdata[which(index1 == zz)]

	print(which(index1 == zz))

	reconsst <- cbind(reconsst,zsst)			# SSTs in matrix form (160 x nreconpt)
	reconsst.a <- cbind(reconsst.a,zclim)		# SST climatology at each point (1 x nreconpt)
}



###################################################################################################
##### 6. ## SCALE ALL THE PALEO DATA (A COUPLE WAYS) and PICK YEAR TO RECONSTRUCT
if (region == "EqPac"){
	smoodat = c(EP.smooth, WP.smooth)
	if (proxy == "both"){
		smoodat = smoodat[-mixed3]
	}
} 

mat = matrix(NaN, lseq, nreconpt)
for (i in 1:nreconpt){
	mat[,i] = smoodat[[i]][,3] - smoodat[[i]][1,2]
}

				

###################################################################################################
##### 3. ## READ IN CONTEMPORARY CURL DATA 
ntime = 65
nclim = 1

# FOR ARABIAN SEA ONLY
xgrid <- seq(35.625, 76.875, by=1.875)
nx <- length(xgrid) # nrows = 33
ygrid <- seq(-6.666573, 31.42808, by=1.904732)
ny <- length(ygrid) # ncols = 23
nsta1 = nx*ny

data.curl <- readBin("datafiles/1949-2013_curl_36.76E-5S.30N_JUNtoSEP_ncepncar.r4",what="numeric",n=(nx*ny*ntime),size=4,endian="swap")
clim.curl <- readBin("datafiles/1981-2010_curl.clim_36.76E-5S.30N_JUNtoSEP_ncepncar.r4",what="numeric",n=(nx*ny*ntime),size=4,endian="swap")


# Create an xygrid of lat and lon for each gridcell.
xygrid <- matrix(0, nrow=(nx*ny), ncol=2)
i = 0
for (iy in 1:ny){
    for (ix in 1:nx){
        i=i+1
        xygrid[i,1] = ygrid[iy]
        xygrid[i,2] = xgrid[ix] 
    }
}

# Set up for masking and determining the index.
data.curl <- array(data=data.curl, dim=c(nx,ny,ntime))
index2 <- 1:(nx*ny)
data1 <- data.curl[,,1]
nanindex <- which(data1 == "NaN")
index3 <- index2[-nanindex]

# write(index3, file="./postfiles/ISM/diagnostics/index-curl.txt")

# Put the data into a matrix of ntime x mlocation
nsta1 = length(index3)
BINdata.curl = matrix(NA, nrow=ntime, ncol=nsta1)
for (i in 1:ntime){
    xcurl = data.curl[,,i]
    x2curl = xcurl[index3]
    BINdata.curl[i,] = x2curl
}

# Read in climatology and scale data.
CLIMdata.curl = clim.curl[index3]
data.curl = sweep(BINdata.curl,2,CLIMdata.curl)


###################################################################################################
## A. (REGRESSING THE FIRST TREND PC SST AND TAKING RESIDUALS)
# Do PCA on 1949-2013 SSTs (data.sst)
zs.sst <- var(IOdata.sst)
zsvd.sst <- svd(zs.sst)
pcs.sst <- t(t(zsvd.sst$u) %*% t(IOdata.sst))
lambdas.sst <- (zsvd.sst$d/sum(zsvd.sst$d))

first = pcs.sst[,1]
data.detrend.curl = matrix(NA,65,nsta1)
for (i in 1:nsta1){
	zzfitcurl = lm(data.curl[,i] ~ first)
	data.detrend.curl[,i] = zzfitcurl$residuals	
}

scdata.detrend.curl = scale(data.detrend.curl)
contemp.sd.curl = attr(scdata.detrend.curl,'scaled:scale')
contemp.mu.curl = attr(scdata.detrend.curl,'scaled:center')

zs.curl <- var(scdata.detrend.curl)
zsvd.curl <- svd(zs.curl)
pcs.curl <- t(t(zsvd.curl$u) %*% t(scdata.detrend.curl))
lambdas.curl <- (zsvd.curl$d/sum(zsvd.curl$d))


###################################################################################################
##### 7. ## SET UP TO PULL OUT GUPTA CURL DATA
ilon = which(xgrid == 58.125)
ilat = which(ygrid == 18.094943)

zfull = seq(1,(nx*ny),1)
zmat = matrix(zfull,nx,ny)
igupta = zmat[ilon,ilat]
iindex = which(index3 == igupta)

missinglon = 18
missinglat = 15
imiss = zmat[missinglon,missinglat]
iindexmiss = which(index3 == imiss)


pilats = c(23.809139, rep(21.904407, 2), rep(19.999675,5), rep(18.094943,4), rep(16.190211,6))
pilons = c(67.500, 65.625, 67.500, 58.125, 60.000, 61.875, 63.750, 65.625, 58.125, 60.000, 61.875, 63.750, 54.375, 56.250, 58.125, 60.000, 61.875, 63.750)
xpilats = c()
xpilons = c()
ipos = c()
for (i in 1:length(pilats)){
	xxlat = which(round(ygrid,6) == pilats[i])
	xpilats = c(xpilats, xxlat)

	xxlon = which(round(xgrid,6) == pilons[i])
	xpilons = c(xpilons, xxlon)

	ipos = c(ipos, which(index3 == zmat[xxlon,xxlat]))
}

###################################################################################################
##### 7. ## Pick the year you want to reconstruct for Holocene
npc = 8 								# how many PCS to keep for CCA
#yrBP = 0
yrsBP = xseq
gupta.holavg.curl = c()
trade.holavg.curl = c()

results = matrix(NaN,length(index3),lseq)
for (i in 1:length(yrsBP)){
	yrBP = yrsBP[i]
	iyrBP = seq(1:lseq)[i]

	missing = which(is.na(mat[iyrBP,]))
	if (length(missing) > 0){
		var = seq(1,nreconpt,1)[-missing]
	} else {
		var = seq(1,nreconpt,1)
	}

	###################################################################################################
	##### 9. ## PCA ON PARTIAL CONTEMPORARY SST FIELD (AND OPTIONAL DIAGNOSTICS)
	zs.p <- var(reconsst[,var])
	zsvd.p <- svd(zs.p)
	pcs.p <- t(t(zsvd.p$u) %*% t(reconsst[,var]))
	lambdas.p <- (zsvd.p$d/sum(zsvd.p$d))


	###################################################################################################
	################################    CURL #######################################################
	###################################################################################################
	Y = pcs.curl[,1:npc]
	X = pcs.p[,1:npc]
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

	cancorln = svd(VV)$d[1:npc]   #canonical correlation
	 
	Fyy = var(Y1) %*% BB

	betahat = solve(t(AA) %*% t(X1)%*% X1 %*% AA) %*% t(AA) %*% t(X1) %*% Y1

	### Pull in paleo data and put in PC space.
	xx.hol = matrix(mat[iyrBP,var],1,length(mat[iyrBP,var]))
	newdata = xx.hol %*% zsvd.p$u[,1:npc]

	### Predict past wind PC field.
	pcmeans = matrix(colMeans(pcs.curl[,-1:-npc]),nrow=1)
	fpc.hol = matrix(NA,1,ncol(pcs.curl))
	scaled = newdata %*% AA %*% betahat												# Reconstructed Hol PCs (scaled)
	unscaled = (scaled * Ysd) + Ymean												# Reconstructed Hol PCs (unscaled: pc mean and sd have been added in)
	fpc.hol[1,1:npc] = unscaled
	fpc.hol[1,-1:-npc] = pcmeans
	recon.hol.curl = fpc.hol %*% t(zsvd.curl$u) 									# Reconstructed Hol curl
	recon.hol.curl.unsc = (recon.hol.curl * contemp.sd.curl) + contemp.mu.curl 		# Note: this is still anomaly space (may need to add CLIM)

	### Predict contemporary wind PC field.
	fpc = matrix(NA,nrow(pcs.curl),ncol(pcs.curl))
	scaled = X1 %*% AA %*% betahat													# Reconstructed Contemp PCs (scaled)
	unscaled = t(t(scaled) * Ysd + Ymean)											# Reconstructed Contemp PCs (unscaled)
	fpc[,1:npc] = unscaled
	fpc[,-1:-npc] = pcmeans[rep(1, nrow(pcs.curl)),]
	recon.contemp.curl = fpc %*% t(zsvd.curl$u)											# Reconstructed Contemp curl
	recon.contemp.curl.unsc = (recon.contemp.curl * contemp.sd.curl) + contemp.mu.curl 	# Reconstructed Contemp curl in anomaly space


	sum = CLIMdata.curl + recon.hol.curl.unsc
	perc = ((recon.hol.curl.unsc)/abs(CLIMdata.curl))*100
	results[,i] = perc

	####################### GUPTA CURL #######################
	gupta.holavg.curl = c(gupta.holavg.curl, recon.hol.curl.unsc[iindex])


	####################### TRADE ROUTE CURL #######################
	tradeavg = mean(recon.hol.curl.unsc[ipos])
	trade.holavg.curl = c(trade.holavg.curl, tradeavg)

}

ka5 = which(xseq == 5)
ka3 = which(xseq == 3)

#coordinates for arabian sea


for (i in ka3:ka5){
	
	xperc = results[,i]

	# MAGNITUDES AS %
	dev.new(width=5, height=5)
	zfull = rep(NaN,(nx*ny))
	zfull[index3] = xperc
	zmat = matrix(zfull,nrow=nx,ncol=ny)
	# image(xgrid,ygrid,zmat,
	# 	ylim=c(15,30), xlim=range(min(xgrid),max(xgrid)),
	# 	xlab="",ylab="", yaxt="n", xaxt="n", cex.axis=1.5, col=myPalette1(100), zlim=c(-20,20)
	# 	)
	image(xgrid,ygrid,zmat,
		ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)),
		xlab="",ylab="", yaxt="n", xaxt="n", cex.axis=1.5, col=myPalette1(100), zlim=c(-20,20)
		)


	image(xgrid,ygrid,zmat,ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)),
		xlab="",ylab="", yaxt="n", xaxt="n", cex.axis=1.5, col="firebrick4", zlim=c(20.1,max(perc)), add=TRUE
		)
	image(xgrid,ygrid,zmat,ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)),
		xlab="",ylab="", yaxt="n", xaxt="n", cex.axis=1.5, col="dodgerblue4", zlim=c(min(perc),-20.1), add=TRUE
		)
	axis(1,seq(30,100,10), cex.axis=1.5)
	axis(2,c(-10,0,10,20,30),cex.axis=1.5, labels=c("-10","0","10","20","30"))
	world(ylim=range(min(ygrid),max(ygrid)),
		xlim=range(min(xgrid),max(xgrid)),
		add=TRUE,
		fill=TRUE,
		border="darkgrey"
	)
	contour(xgrid, ygrid, zmat, add = TRUE, levels = seq(-20,20,10), lwd = 1.2, labcex=1, col="grey50", drawlabels=FALSE)
	# contour(xgrid, ygrid, zmat, add = TRUE, levels = 0, lwd = 1.2, labcex=1, col="black", drawlabels=FALSE)
	# points(57.36, 18.05, pch=21, col="black", bg="#009E73", cex=1.5)

}




####################################################################################
####################################################################################
########## GUPTA COMPARISON
## GUPTA. Holocene GB.
paleobull <- matrix(scan("/Users/emilygill/Documents/University of Colorado/PHD Research/4. Reconstruction/gupta2003.txt"),ncol=3,byrow=TRUE)
paleobull <- cbind(paleobull,scale(paleobull[,3]))

par(mar=c(5, 12, 4, 4) + 0.1)
x = paleobull[34:91,2]
y = paleobull[34:91,3]
zfit = locfit(y ~ x, deg=3, kern="bisq")
yest = predict(zfit,x)
newx = xseq[which(xseq == 2):lseq]*1000

col1 <- col2alpha("#0072B2",alpha=0.1)

plot(x, y, type="l", ylab="", xlab="",xlim=rev(range(x)), col="grey25", axes=F, ylim=c(0,30))
lines(x,yest,lty=2,col="grey25",xlim=rev(range(x)))
axis(1, cex.axis=1.25, at=c(10000,9000,8000,7000,6000,5000,4000,3000,2000,1000,0))
axis(2, ylim = c(range(y)), col="grey25", lwd=2, cex.axis=1.25)
mtext(2, text="G. bulloides % [Gupta et al. 2003]",line=2, cex=1.25, col="grey25")
mtext(1, text="ka BP", line=2, cex=1.25)

ODP.curl = gupta.holavg.curl[which(xseq == 2):lseq]
ODP.curl = trade.holavg.curl[which(xseq == 2):lseq]
sum = CLIMdata.curl[iindex] + ODP.curl
ODP.perc = ((sum - CLIMdata.curl[iindex])/CLIMdata.curl[iindex])*100
par(new = T)
#plot(newx,gupta.holavg.curl, type="l",axes = F, lwd=2, xlab = NA, ylab = NA, xlim=rev(range(newx)), col="#0072B2", ylim = c(-2,4))
#plot(newx,gupta.holavg.curl, type="l",axes = F, lwd=2, xlab = NA, ylab = NA, xlim=c(10877,0), col="#0072B2", ylim = c(-1,4))
plot(newx,ODP.perc, type="l",axes = F, lwd=2, xlab = NA, ylab = NA, xlim=rev(range(newx)), col="#0072B2", ylim = c(0,20))
#polygon(c(rev(newx),newx),c(rev(gupta.holavg.u - (1.95*mse.mag.gupta)), gupta.holavg.mag + (1.95*mse.gupta)),col=col1,border=col1)
#lines(newx, gupta.holavg.mag - mse.mag.gupta, xlim=rev(range(newx)), col="#0072B2", lty=2)
#lines(newx, gupta.holavg.mag + mse.mag.gupta, xlim=rev(range(newx)), col="#0072B2", lty=2)
axis(4, ylim=c(0,35),lwd=2,col="#0072B2", cex.axis=1.25)
mtext(4,text="Percentage Change in Mean Summer Curl", line=2,col="#0072B2", cex=1.25)




# MAGNITUDES AS %
dev.new(width=5, height=5)
zfull = rep(NaN,(nx*ny))
zfull[index3] = xperc
zmat = matrix(zfull,nrow=nx,ncol=ny)
# image(xgrid,ygrid,zmat,
# 	ylim=c(15,30), xlim=range(min(xgrid),max(xgrid)),
# 	xlab="",ylab="", yaxt="n", xaxt="n", cex.axis=1.5, col=myPalette1(100), zlim=c(-20,20)
# 	)
image(xgrid,ygrid,zmat,
	ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)),
	xlab="",ylab="", yaxt="n", xaxt="n", cex.axis=1.5, col=myPalette1(100), zlim=c(-20,20)
	)

## To test if I'm reading the curl data in correctly.
dev.new(width=9,height=7)
zfull = rep(NaN,(nx*ny))
zfull[index3] = apply(BINdata.curl[33:62,],2,mean)
zmat = matrix(zfull,nrow=nx,ncol=ny)
image.plot(xgrid,ygrid,zmat,
	ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)),
	xlab="",ylab="", main="",
	cex.axis=1.1,
	zlim=c(-70,70), col=myPalette4(100)
	)
axis(1, cex.axis=1.5)
axis(2, at = c(-10,30,10), cex.axis=1.5)
mapnames <- map("world2",
	xlim=c(min(ygrid), max(xgrid)), ylim=c(min(ygrid), max(ygrid)), 
	boundary=TRUE, interior=TRUE, add=TRUE, lwd=1.5, col="grey75", fill=TRUE
)
contour(xgrid, ygrid, zmat, add = TRUE, levels = seq(-60,60,10), lwd = 1.75, labcex=1.25, col="black")
## To plot the points for trade route:
points(pilons, pilats, pch=21, cex=1.25, bg="#0072B2")

## To plot the gupta point:
guptax = 58.125
guptay = 18.094943
points(guptax, guptay)

