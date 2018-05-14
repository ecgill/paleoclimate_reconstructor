rm(list=ls())
setwd("/Users/emilygill/Documents/University of Colorado/PHD Research/6. SST Reconstruction/REDO Sep2015/")


##### 1. ## SOURCE FUNCTIONS AND READ IN LIBRARIES
source("/Users/emilygill/Documents/University of Colorado/PHD Research/6. SST Reconstruction/BINxygrid.R")
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
myPalette1 <- colorRampPalette(rev(brewer.pal(9, "RdBu")), space="Lab")
myPalette2 <- colorRampPalette(rev(brewer.pal(9, "Greens")), space="Lab")
myPalette3 <- colorRampPalette(rev(brewer.pal(9, "YlOrRd")), space="Lab")
myPalette4 <- colorRampPalette(rev(brewer.pal(9, "BrBG")), space="Lab")
myPalette5 <- colorRampPalette(rev(brewer.pal(9, "PRGn")), space="Lab")

###################################################################################################
##### *** PICK AN AREA TO DO THE RECONSTRUCTION OVER
region = "EqPac"						# EqPac or IO
proxy = "both"							# MgCa or Uk
season = "annual"						# winter or summer


###################################################################################################
##### 2. ## READ IN CONTEMPORARY DATA 
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

# Data for overlapping wind period of 1949-2013
data.sst = data[96:160,]



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




###################################################################################################
##### 4. ## SMOOTH AND PLOT PALEO RECONSTRUCTIONS
sst.smooth = gcv.recon(sst.raw$EPlist, sst.raw$WPlist, sst.raw$IOlist, 
	EPdeg=2, WPdeg=2, IOdeg=2, base=TRUE)
# pEP = plot.recon(sst.raw$EPlist,sst.smooth$EP, map=TRUE, region="EP",sst.raw$EPcoord)
# pWP = plot.recon(sst.raw$WPlist,sst.smooth$WP, map=TRUE, region="WP",sst.raw$WPcoord)
# pIO = plot.recon(sst.raw$IOlist,sst.smooth$IO, map=TRUE, region="IO",sst.raw$IOcoord)




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
# else if (region == "IO"){
# 	reconpt <- rbind(sst.raw$IOcoord)
# }

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
	smoodat = c(sst.smooth$EP, sst.smooth$WP)
	if (proxy == "both"){
		smoodat = smoodat[-mixed3]
	}
} 
# else if (region == "IO"){
# 	smoodat = c(sst.smooth$IO)
# }

# sst.scale = scale.recon(smoodat, nreconpt, reconsst.mu, reconsst.sig)
sst.scale = scale.recon(smoodat, nreconpt, reconsst.a)
mat.unsc = sst.scale$unscaled			# Not scaled, just put into a matrix
# mat.sc = sst.scale$scaled				# Put into matrix and scaled using itself
mat.anom = sst.scale$anom				# Put into matrix and scaled using the climatology

mat.base = sst.scale$base
mat.baselim = sst.scale$baselim
# Pick which scaled series you want to use
mat = mat.baselim					






###################################################################################################
##### 3. ## READ IN CONTEMPORARY WIND DATA 
ntime = 65
nclim = 1

xgrid <- seq(100,300,by=2.5)
nx <- length(xgrid) # nrows = 109
ygrid <- seq(-10,10,by=2.5)
ny <- length(ygrid) # ncols = 21
nsta1 = nx*ny

data.u <- readBin("datafiles/1949-2013_U_100.300E-10S.10N_MAYtoAPR_ncepncar.r4",what="numeric",n=(nx*ny*ntime),size=4,endian="swap")	
clim.u <- readBin("datafiles/1981-2010_U.clim_100.300E-10S.10N_MAYtoAPR_ncepncar.r4",what="numeric",n=(nx*ny*nclim),size=4,endian="swap")
data.v <- readBin("datafiles/1949-2013_V_100.300E-10S.10N_MAYtoAPR_ncepncar.r4",what="numeric",n=(nx*ny*ntime),size=4,endian="swap")	
clim.v <- readBin("datafiles/1981-2010_V.clim_100.300E-10S.10N_MAYtoAPR_ncepncar.r4",what="numeric",n=(nx*ny*nclim),size=4,endian="swap")

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
data.u <- array(data=data.u, dim=c(nx,ny,ntime))
data.v <- array(data=data.v, dim=c(nx,ny,ntime))

index1 <- 1:(nx*ny)
data1 <- data.u[,,1]

# Put the data into a matrix of ntime x mlocation
BINdata.u = matrix(NA, nrow=ntime, ncol=nsta1)
BINdata.v = matrix(NA, nrow=ntime, ncol=nsta1)
for (i in 1:ntime){
    udata = data.u[,,i]
    data2u = udata[index1]
    BINdata.u[i,] = data2u

    vdata = data.v[,,i]
    data2v = vdata[index1]
    BINdata.v[i,] = data2v

}

# Read in climatology and scale data.
CLIMdata.u = clim.u[index1]
data.u = sweep(BINdata.u,2,CLIMdata.u)

CLIMdata.v = clim.v[index1]
data.v = sweep(BINdata.v,2,CLIMdata.v)


###################################################################################################
##### 4. ## PCA ON FULL CONTEMPORARY WIND FIELD (AND OPTIONAL DIAGNOSTICS)

## U-WIND
zs.u <- var(data.u)
zsvd.u <- svd(zs.u)
pcs.u <- t(t(zsvd.u$u) %*% t(data.u))
lambdas.u <- (zsvd.u$d/sum(zsvd.u$d))

## V-WIND
zs.v <- var(data.v)
zsvd.v <- svd(zs.v)
pcs.v <- t(t(zsvd.v$u) %*% t(data.v))
lambdas.v <- (zsvd.v$d/sum(zsvd.v$d))

# write.table(zsvd.u$u, 
# 	file=paste(paste("./postfiles/", region, "/CCA-PC6/diagnostics/", region, "_", proxy, "_", season, "_zvsd.u$u", sep=""),
# 	sep="\t"))
# write.table(pcs.u, 
# 	file=paste(paste("./postfiles/", region, "/CCA-PC6/diagnostics/", region, "_", proxy, "_", season, "_pcs.u", sep=""),
# 	sep="\t"))
# write(lambdas.u, 
# 	file=paste(paste("./postfiles/", region, "/CCA-PC6/diagnostics/", region, "_", proxy, "_", season, "_lambdas.u", sep=""),
# 	sep="\t"))

###################################################################################################
##### 7. ## Pick the year you want to reconstruct for Holocene
npc = 6 								# how many PCS to keep for CCA
#yrBP = 0
yrsBP = 0:10
for (i in 1:length(yrsBP)){
	print(i)
	yrBP = yrsBP[i]

	missing = which(is.na(mat[yrBP+1,]))
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
	##### 7. ## PERFORM CCA to reconstruct U-WIND
	Y = pcs.u[,1:npc]
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
	xx.hol = matrix(mat[yrBP+1,var],1,length(mat[yrBP+1,var]))
	newdata = xx.hol %*% zsvd.p$u[,1:npc]

	### Predict past wind PC field.
	pcmeans = matrix(colMeans(pcs.u[,-1:-npc]),nrow=1)
	fpc.hol = matrix(NA,1,ncol(pcs.u))
	scaled = newdata %*% AA %*% betahat
	unscaled = (scaled * Ysd) + Ymean
	fpc.hol[1,1:npc] = unscaled
	fpc.hol[1,-1:-npc] = pcmeans
	recon.hol.u = fpc.hol %*% t(zsvd.u$u)

	### Predict contemporary wind PC field.
	fpc = matrix(NA,nrow(pcs.u),ncol(pcs.u))
	scaled = X1 %*% AA %*% betahat
	unscaled = t(t(scaled) * Ysd + Ymean)
	fpc[,1:npc] = unscaled
	fpc[,-1:-npc] = pcmeans[rep(1, nrow(pcs.u)),]
	recon.contemp.u = fpc %*% t(zsvd.u$u)


	### Statistics
	sumdiff = apply((data.u - recon.contemp.u)^2,2,sum)
	sumYraw = apply(data.u^2,2,sum)
	beta = 1 - (sumdiff/sumYraw)

	cor <- c()
	for (j in 1:nsta1){
		xcor <- cor(recon.contemp.u[,j],data.u[,j])^2
		cor <- c(cor,xcor)
	}


	# ###################################################################################################
	# ##### 7. ## PERFORM CCA to reconstruct U-WIND
	# X = pcs.v[,1:npc]
	# Y = pcs.p[,1:npc]
	# N = length(Y[,1])

	# X1 = scale(X)
	# Xmean = attr(X1,'scaled:center')
	# Xsd = attr(X1,'scaled:scale')
	# Y1 = scale(Y)

	# Qx1 = qr.Q(qr(X1))
	# Qy1 = qr.Q(qr(Y1))
	# T11 = qr.R(qr(X1))
	# T22 = qr.R(qr(Y1))

	# VV=t(Qx1) %*% Qy1
	# BB = svd(VV)$v
	# AA = svd(VV)$u 

	# BB = solve(T22) %*% svd(VV)$v * sqrt(N-1)
	# wm1=Y1 %*% BB

	# AA = solve(T11) %*% svd(VV)$u * sqrt(N-1)
	# vm1 = X1 %*% AA

	# cancorln = svd(VV)$d[1:npc]   #canonical correlation
	 
	# Fyy = var(Y1) %*% BB

	# betahat = solve(t(AA) %*% t(X1)%*% X1 %*% AA) %*% t(AA) %*% t(X1) %*% Y1

	# ### Pull in paleo data and put in PC space.
	# xx.hol = matrix(mat[yrBP+1,var],1,length(mat[yrBP+1,var]))
	# newdata = xx.hol %*% zsvd.p$u[,1:npc]

	# ### Predict past wind PC field.
	# pcmeans = matrix(colMeans(pcs.v[,-1:-npc]),nrow=1)
	# fpc.hol = matrix(NA,1,ncol(pcs.v))
	# scaled = newdata %*% AA %*% betahat
	# unscaled = (scaled * Xsd) + Xmean
	# fpc.hol[1,1:npc] = unscaled
	# fpc.hol[1,-1:-npc] = pcmeans
	# recon.hol.v = fpc.hol %*% t(zsvd.v$u)

	# ### Predict contemporary wind PC field.
	# fpc = matrix(NA,nrow(pcs.v),ncol(pcs.v))
	# scaled = X1 %*% AA %*% betahat
	# unscaled = t(t(scaled) * Xsd + Xmean)
	# fpc[,1:npc] = unscaled
	# fpc[,-1:-npc] = pcmeans[rep(1, nrow(pcs.v)),]
	# recon.contemp.v = fpc %*% t(zsvd.v$u)


	write.table(recon.hol.u, 
		file=paste(paste("./postfiles/EqPac/CCA-PC", npc, "/", region, "_", proxy, "_", season, "_", npc, "npc", "_", yrBP, "ka_u.txt", sep=""),
		sep="\t"))
	# # write.table(recon.hol.v, 
	# # 	file=paste(paste("./post-files/EqPac/CCA-WINDS/PC", npc, "/", region, "_", proxy, "_", season, "_", npc, "npc", "_", yrBP, "ka_v.txt", sep=""),
	# # 	sep="\t"))

	if (i == 11){
		write.table(recon.contemp.u, 
			file=paste(paste("./postfiles/EqPac/CCA-PC", npc, "/", region, "_", proxy, "_", season, "_", npc, "npc", "_", yrBP, "ka_1949-2013.txt", sep=""),
			sep="\t"))

		write.table(data.u, 
			file=paste(paste("./postfiles/EqPac/CCA-PC", npc, "/", region, "_", proxy, "_", season, "_", npc, "npc", "_actual_1949-2013.txt", sep=""),
			sep="\t"))
	}

}

validate = plot.validate2(var="wind",npc=as.character(npc),region=region, proxy=proxy, 
	season=season, index1, xgrid, ygrid, reconpt, zzpt, mat)



###################################################################################################

dev.new(width=10, height=3)
par(mfrow = c(1, 1), mar = c(4, 4, 4, 1))
zfull = rep(NaN,(nx*ny))
zfull[index1] = recon.hol.u 						# add in [index1] if youre plotting India rainfall or SST
zmat = matrix(zfull,nrow=nx,ncol=ny)
image.plot(xgrid,ygrid,zmat,
	ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)),
	xlab="",ylab="", main=expression(paste("Calibration ", beta)),
	cex.axis=1.1, col=myPalette4(100),
	zlim=c(-12,12)
	)
mapnames <- map("world2", 
	xlim=c(min(ygrid), max(xgrid)), ylim=c(min(ygrid), max(ygrid)), 
	boundary=TRUE, interior=TRUE, add=TRUE, lwd=1.5, col="grey50"
)
contour(xgrid, ygrid, zmat, add = TRUE, levels = seq(-12,12,2), lwd = 2, labcex=1.25)


### Statistics
dev.new(width=10, height=3)
par(mfrow = c(1, 1), mar = c(4, 4, 4, 1))
zfull = rep(NaN,(nx*ny))
zfull[index1] = beta 						# add in [index1] if youre plotting India rainfall or SST
zmat = matrix(zfull,nrow=nx,ncol=ny)
image.plot(xgrid,ygrid,zmat,
	ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)),
	xlab="",ylab="", main=expression(paste("Calibration ", beta)),
	cex.axis=1.1, col=rev(myPalette3(100)),
	zlim=c(-1,1)
	)
mapnames <- map("world2", 
	xlim=c(min(ygrid), max(xgrid)), ylim=c(min(ygrid), max(ygrid)), 
	boundary=TRUE, interior=TRUE, add=TRUE, lwd=1.5, col="grey50"
)
contour(xgrid, ygrid, zmat, add = TRUE, levels = seq(-1,1,0.05), lwd = 2, labcex=1.25)


dev.new(width=10, height=3)
par(mfrow = c(1, 1), mar = c(4, 4, 4, 1))
zfull = rep(NaN,(nx*ny))
zfull[index1] = cor 						# add in [index1] if youre plotting India rainfall
zmat = matrix(zfull,nrow=nx,ncol=ny)
image.plot(xgrid,ygrid,zmat,
	ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)),
	xlab="",ylab="", main = expression(paste("Calibration ", R^{2})), 
	cex.axis=1.1, col=rev(myPalette2(100)),
	zlim=c(0,1)
	)
mapnames <- map("world2", 
	xlim=c(min(ygrid), max(xgrid)), ylim=c(min(ygrid), max(ygrid)), 
	boundary=TRUE, interior=TRUE, col="grey50", add=TRUE
)
contour(xgrid, ygrid, zmat, add = TRUE, levels = seq(-1,1,0.1), lwd = 2, labcex=1.25)


dev.new(width=10, height=4)
par(mfrow = c(2, 1), mar = c(2,2,2,2))
# BETA
zfull = rep(NaN,(nx*ny))
zfull[index1] = beta
zmat = matrix(zfull,nrow=nx,ncol=ny)
image(xgrid,ygrid,zmat,
	ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)),
	xlab="",ylab="", main=expression(paste("Zonal Wind Calibration ", beta)),
	cex.axis=1.1, yaxt="n", col=rev(myPalette3(100)),
	zlim=c(-1,1)
	)
axis(2, at = c(-10,0,10))
mapnames <- map("world2", 
	xlim=c(min(ygrid), max(xgrid)), ylim=c(min(ygrid), max(ygrid)), 
	boundary=TRUE, interior=TRUE, add=TRUE, lwd=1.5, col="grey50"
)
contour(xgrid, ygrid, zmat, add = TRUE, levels = seq(-1,1,0.1), lwd = 2, labcex=1.25)

# R^2
zfull = rep(NaN,(nx*ny))
zfull[index1] = cor
zmat = matrix(zfull,nrow=nx,ncol=ny)
image(xgrid,ygrid,zmat,
	ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)),
	xlab="",ylab="", main=expression(paste("Zonal Wind Calibration ", R^{2})),
	cex.axis=1.1, yaxt="n", col=rev(myPalette2(100)),
	zlim=c(0,1)
	)
axis(2, at = c(-10,0,10))
mapnames <- map("world2", 
	xlim=c(min(ygrid), max(xgrid)), ylim=c(min(ygrid), max(ygrid)), 
	boundary=TRUE, interior=TRUE, add=TRUE, lwd=1.5, col="grey50"
)
contour(xgrid, ygrid, zmat, add = TRUE, levels = seq(0,1,0.1), lwd = 2, labcex=1.25)

### ONLY USE FOR THE LEGENDS
dev.new(height=4,width=10)
image.plot(xgrid,ygrid,zmat,
	ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)),
	xlab="",ylab="", main=expression(paste("Zonal Wind Calibration ", beta)),
	cex.axis=1.1, yaxt="n", col=rev(myPalette3(100)),
	zlim=c(-1,1)
	)
dev.new(height=4,width=10)
image.plot(xgrid,ygrid,zmat,
	ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)),
	xlab="",ylab="", main=expression(paste("Zonal Wind Calibration ", R^{2})),
	cex.axis=1.1, yaxt="n", col=rev(myPalette2(100)),
	zlim=c(0,1)
	)


#### PLOT COMPARISON YEARS
iix = which(1949:2013 == 1988)
year = "1988-1989"
max = 6.5
min = -6.5

dev.new(width=10, height=4)
par(mfrow = c(2, 1), mar = c(2,2,2,2))
# Actual
zfull = rep(NaN,(nx*ny))
zfull[index1] = data.u[iix,]					# add in [index1] if youre plotting India rainfall or SST
zmat = matrix(zfull,nrow=nx,ncol=ny)
image(xgrid,ygrid,zmat,
	ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)),
	xlab="", ylab="", main=paste("Actual", year), 
	cex.axis=1.1, yaxt="n", col=myPalette4(100),
	zlim=c(min,max)
	)
axis(2, at = c(-10,0,10))
mapnames <- map("world2", 
	xlim=c(min(ygrid), max(xgrid)), ylim=c(min(ygrid), max(ygrid)), 
	boundary=TRUE, interior=TRUE, add=TRUE, lwd=1.5, col="grey50"
)
contour(xgrid, ygrid, zmat, add = TRUE, levels = seq(min,max,1), lwd = 2, labcex=1.25)


# Reconstructed
zfull = rep(NaN,(nx*ny))
zfull[index1] = recon.contemp.u[iix,] 						# add in [index1] if youre plotting India rainfall or SST
zmat = matrix(zfull,nrow=nx,ncol=ny)
image.plot(xgrid,ygrid,zmat,
	ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)),
	xlab="",ylab="", main=paste("Reconstructed", year),
	cex.axis=1.1, yaxt="n", col=myPalette4(100),
	zlim=c(min,max)
	)
axis(2, at = c(-10,0,10))
mapnames <- map("world2", 
	xlim=c(min(ygrid), max(xgrid)), ylim=c(min(ygrid), max(ygrid)), 
	boundary=TRUE, interior=TRUE, add=TRUE, lwd=1.5, col="grey50"
)
contour(xgrid, ygrid, zmat, add = TRUE, levels = seq(min,max,1), lwd = 2, labcex=1.25)






###################################################################################################
##### 7. ## PLOT WITHOUT VECTORS
years = seq(2,10,2)
prow = length(years) + 1
recon.zonal = c()

dev.new(height=5,width=10)
par(mfrow = c(prow, 1), mar = c(0.5, 2, 0.5, 0.5), oma=c(2,1,1,1))	

readpresent = read.table(paste("./postfiles/EqPac/CCA-PC", npc, "/", region, "_", proxy, "_", season, "_", npc, "npc_10ka_1949-2013.txt", sep=""),header=TRUE)
present = apply(readpresent,2,mean)

zfull = rep(NaN,(nx*ny))
zfull[index1] = as.numeric(unlist(present))
zmat = matrix(zfull,nrow=nx,ncol=ny)
image(xgrid,ygrid,zmat,
	ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)), axes=FALSE, 
	col=myPalette4(100), zlim=c(-10,10)
	)
	contour(xgrid, ygrid, zmat, add = TRUE, levels = seq(-10,10,1), lwd = 1.2, labcex=1.25)
mapnames <- map("world2", 
	xlim=c(min(ygrid), max(xgrid)), ylim=c(min(ygrid), max(ygrid)), 
	boundary=TRUE, interior=TRUE, col="black", add=TRUE
)
	axis(2,seq(-8,8,by=4),cex.axis=2)

for (i in 1:length(years)){
	test = read.table(paste("./postfiles/EqPac/CCA-PC", npc, "/", region, "_", proxy, "_", season, "_", npc, "npc", "_", years[i], "ka_u.txt", sep=""),header=TRUE)
	recon.zonal = rbind(recon.zonal,test)
	
	zfull = rep(NaN,(nx*ny))
	zfull[index1] = as.numeric(unlist(test))
	zmat = matrix(zfull,nrow=nx,ncol=ny)
	image(xgrid,ygrid,zmat,
		ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)), axes=FALSE, 
		col=myPalette4(100), zlim=c(-11,11)
		)
	mapnames <- map("world2", 
		xlim=c(min(ygrid), max(xgrid)), ylim=c(min(ygrid), max(ygrid)), 
		boundary=TRUE, interior=TRUE, col="grey50", add=TRUE
	)
	contour(xgrid, ygrid, zmat, add = TRUE, levels = seq(-10,10,1), lwd = 1.2, labcex=1.25)
	axis(2,seq(-8,8,by=4),cex.axis=2)

	if (i == length(years)) {
         axis(1,c(100,140,180,220,260,300),labels=c("100","140","180","-140","-100","-60"),cex.axis=2)
    }
}
