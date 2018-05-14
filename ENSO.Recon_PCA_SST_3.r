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
myPalette1 <- colorRampPalette(rev(brewer.pal(9, "RdBu")), space="Lab")
myPalette2 <- colorRampPalette(rev(brewer.pal(9, "Greens")), space="Lab")
myPalette3 <- colorRampPalette(rev(brewer.pal(9, "YlOrRd")), space="Lab")
myPalette4 <- colorRampPalette(rev(brewer.pal(9, "Greys")), space="Lab")
myPalette5 <- colorRampPalette(rev(brewer.pal(9, "GnBu")), space="Lab")
myPalette6 <- colorRampPalette(rev(brewer.pal(9, "Purples")), space="Lab")
myPalette7 <- colorRampPalette(rev(brewer.pal(9, "Blues")), space="Lab")
spectral <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")



###################################################################################################
##### *** PICK AN AREA TO DO THE RECONSTRUCTION OVER
region = "EqPac"						# EqPac or IO
proxy = "both"							# MgCa or Uk
season = "annual"						# annual, winter or summer


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

# ## Plot to check if reading the right data.
# dev.new(width=10, height=3)
# par(mfrow = c(1, 1), mar = c(4, 4, 4, 1))
# zfull = rep(NaN,(nx*ny))
# zfull[index1] = data.un[1,] 						# add in [index1] if youre plotting India rainfall or SST
# zmat = matrix(zfull,nrow=nx,ncol=ny)
# image.plot(xgrid,ygrid,zmat,
# 	ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)),
# 	xlab="",ylab="", main="",
# 	yaxt="n", col=myPalette1(100),
# 	zlim=c(20,25.01)
# 	)
# contour(xgrid, ygrid, zmat, add = TRUE, levels = seq(0,1,0.1), lwd = 2, labcex=1.25)
# mapnames <- map("world2", 
# 	xlim=c(min(ygrid), max(xgrid)), ylim=c(min(ygrid), max(ygrid)), 
# 	boundary=TRUE, interior=TRUE, col="black", add=TRUE
# )
# axis(2, at = c(-10,0,10))

# write(index1, 
# 	file=paste(paste("./postfiles/", region, "/PC3/diagnostics/", region, "_", proxy, "_", season, "_index1", sep=""),
# 	sep="\t"))

###################################################################################################
##### 8. ## PCA ON FULL CONTEMPORARY SST FIELD (AND OPTIONAL DIAGNOSTICS)
zs <- var(data)
zsvd <- svd(zs)
pcs <- t(t(zsvd$u) %*% t(data))
lambdas <- (zsvd$d/sum(zsvd$d))

# write.table(zsvd$u, 
# 	file=paste(paste("./postfiles/", region, "/PC3/diagnostics/", region, "_", proxy, "_", season, "_zvsd$u", sep=""),
# 	sep="\t"))
# write.table(pcs, 
# 	file=paste(paste("./postfiles/", region, "/PC3/diagnostics/", region, "_", proxy, "_", season, "_pcs", sep=""),
# 	sep="\t"))
# write(lambdas, 
# 	file=paste(paste("./postfiles/", region, "/PC3/diagnostics/", region, "_", proxy, "_", season, "_lambdas", sep=""),
# 	sep="\t"))

###################################################################################################
##### 3. ## READ IN PALEO DATA POINT RECONSTRUCTIONS RECONSTRUCTIONS
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

# # Calculate the number of records per 1000 years
# EP = sst.raw$EPlist
# EPres = 1:length(EP)
# for (i in 1:length(EP)){
# 	Nres = dim(EP[[i]])[1]
# 	EPres[i] = Nres/(EP[[i]][Nres,2] - EP[[i]][1,2])
# }

# WP = sst.raw$WPlist
# WPres = 1:length(WP)
# for (i in 1:length(WP)){
# 	Nres = dim(WP[[i]])[1]
# 	WPres[i] = Nres/(WP[[i]][Nres,2] - WP[[i]][1,2])	
# }


###################################################################################################
##### 4. ## SMOOTH AND PLOT PALEO RECONSTRUCTIONS
sst.smooth = gcv.recon(sst.raw$EPlist, sst.raw$WPlist, sst.raw$IOlist, 
	EPdeg=2, WPdeg=2, IOdeg=2, base=TRUE)
# pEP = plot.recon(sst.raw$EPlist,sst.smooth$EP, map=TRUE, region="EP", proxy=proxy, sst.raw$EPcoord)
# pWP = plot.recon(sst.raw$WPlist,sst.smooth$WP, map=TRUE, region="WP", proxy=proxy, sst.raw$WPcoord)
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
	
} else if (region == "IO"){
	reconpt <- rbind(sst.raw$IOcoord)
}

# recordno = c(27:39,22:26,15:21,1:14)
# recordno = recordno[-mixed3]
# recordnopt = cbind(recordno,unname(reconpt))
# write.table(recordnopt, 
# 	file=paste(paste("./postfiles/", region, "/PC3/diagnostics/", region, "_", proxy, "_", season, "_recordnopt", sep=""),
# 	sep="\t"))

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
	zsst <- data[,which(index1 == zz)]
	zclim <- climdata[which(index1 == zz)]

	print(which(index1 == zz))

	reconsst <- cbind(reconsst,zsst)			# SSTs in matrix form (160 x nreconpt)
	reconsst.a <- cbind(reconsst.a,zclim)		# SST climatology at each point (1 x nreconpt)
}

#write.table(reconsst, file="/Users/emilygill/Documents/University of Colorado/PostDoc Research/5. Bayesian/test_PC_PacSST/reconsst.txt", sep="\t")

# write(zzpt, 
# 	file=paste(paste("./postfiles/", region, "/PC3/diagnostics/", region, "_", proxy, "_", season, "_zzpt.txt", sep=""),
# 	sep="\t"))

### Optional Plot with dots for climatology
## PLOT THE SMOOTH VALUES

# plotregion = "EP"

# if (plotregion == "EP") {
# 	smooth = sst.smooth$EP
# 	yrange = c(22,30)

# 	if (proxy == "Uk") {
# 		climpts = reconsst.a[1:13]
# 	} else if (proxy == "MgCa") {
# 		climpts = reconsst.a[1:5]
# 	}

# } else {
# 	smooth = sst.smooth$WP
# 	yrange = c(25,33)

# 	if (proxy == "Uk") {
# 		climpts = reconsst.a[14:20]
# 	} else if (proxy == "MgCa") {
# 		climpts = reconsst.a[6:19]
# 	}
# }

# n = length(smooth)

# spectral <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
# cbPalette = rev(spectral(n))
# dev.new()
# par(mfrow = c(1, 1), mar = c(0.75, 0.75, 0.75, 0.75), oma = c(4, 4, 3, 1))		
# for (i in 1:(n)){
# 	xxs = data.frame(smooth[i])
# 	names(xxs) = c("X","Y")

# 	if (i == 1){
# 		plot(xxs$X, xxs$Y, type="l", lty=2, xlim=rev(c(0,10)), ylim=yrange, axes=FALSE, col=cbPalette[1], lwd=3)
# #			lines(xxs$X, xxs$Y, lty=2, col=cbPalette[1])
# 		points(0,climpts[1],col=cbPalette[1])

# 	} else {
# 		lines(xxs$X, xxs$Y, lty=2, col=cbPalette[i],lwd=3)
# #			lines(xxs$X, xxs$Y, col=cbPalette[i], lty=2)
# 		points(0,climpts[i],col=cbPalette[i])

# 	}

# 	if (i == n){
# 		axis(1, at=seq(10, 0, by = -1), las=1, cex.axis=1.5)
# 		axis(2, at=seq(yrange[1], yrange[2], by = 1), cex.axis=1.5)
# 		mtext("ka BP", side=1, line=3, cex=1)
# 		mtext(~degree~C, side=2, line=3, cex=1)
# 	}
# }


###################################################################################################
##### 6. ## SCALE ALL THE PALEO DATA (A COUPLE WAYS) and PICK YEAR TO RECONSTRUCT
if (region == "EqPac"){
	smoodat = c(sst.smooth$EP, sst.smooth$WP)
#	if (proxy == "both"){
		smoodat = smoodat[-mixed3]
#	}
} else if (region == "IO"){
	smoodat = c(sst.smooth$IO)
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

# write.table(mat, 
# 	file=paste(paste("./postfiles/", region, "/PC3/diagnostics/", region, "_", proxy, "_", season, "_mat.txt", sep=""),
# 	sep="\t"))				

# yrBP = 10

yrsBP = 0:10
for (i in 1:length(yrsBP)){
	yrBP = yrsBP[i]

	###################################################################################################
	##### 7. ## Pick the year you want to reconstruct for Holocene
	missing = which(is.na(mat[yrBP+1,]))
	if (length(missing) > 0){
		var = seq(1,nreconpt,1)[-missing]
	} else {
		var = seq(1,nreconpt,1)
	}

	###################################################################################################
	##### 9. ## PCA ON PARTIAL CONTEMPORARY SST FIELD (AND OPTIONAL DIAGNOSTICS)
	zs.p <- var(reconsst[,var])
	
	# xrecon = reconsst[,var]
	# N = dim(xrecon)[1]
	# P = dim(xrecon)[2]
	# X11 = t(t(xrecon) - apply(xrecon,2,mean))
	# X11 = unname(X11)
	# W = matrix(0,P,P)
	# # diag(W) = 1
	# # X22nw = t(X11) %*% X11 %*% W * (1/(N-1))
	# weights = c(1,1,1,1,1,1,0.6,1.2,1.2,1.2,1.2,1,1,1,1,1,1,1,1,1.2,0.5,0.7,0.6,1.2,1,0.6)[var]
	# diag(W) = 1/weights
	# zs.p = t(X11) %*% X11 %*% W * (1/(N-1))

	zsvd.p <- svd(zs.p)
	pcs.p <- t(t(zsvd.p$u) %*% t(reconsst[,var]))
	lambdas.p <- (zsvd.p$d/sum(zsvd.p$d))

	# write.table(zsvd.p$u, 
	# file=paste(paste("./postfiles/", region, "/PC3/diagnostics/", region, "_", proxy, "_", season, "_zvsd.p$u", sep=""),
	# sep="\t"))
	# write.table(pcs.p, 
	# file=paste(paste("./postfiles/", region, "/PC3/diagnostics/", region, "_", proxy, "_", season, "_pcs.p", sep=""),
	# sep="\t"))
	# write(lambdas.p, 
	# file=paste(paste("./postfiles/", region, "/PC3/diagnostics/", region, "_", proxy, "_", season, "_lambdas.p", sep=""),
	# sep="\t"))

	###################################################################################################
	##### 10. ## FIT THE MODEL
	#### Fit each retained PC of Y (full SST) with the retained PCs of X (holocene SSTs)
	xx.hol = matrix(mat[yrBP+1,var],1,length(mat[yrBP+1,var]))
	newdata = as.data.frame(xx.hol %*% zsvd.p$u)
	# write.table(newdata, 
	# file=paste(paste("./postfiles/", region, "/PC3/diagnostics/", region, "_", proxy, "_", season, "_newdata10ka", sep=""),
	# sep="\t"))

	X = as.data.frame(pcs.p[,1:length(mat[yrBP+1,var])])
	fpc = matrix(NA,nrow(pcs),ncol(pcs))
	fpc.hol = matrix(NA,1,ncol(pcs))
	fpc.hol.boot = matrix(NA,500,ncol(pcs))
	npc = 3
	npcr = 2
	for (j in 1:npc){
		Y <- as.matrix(pcs[,j])

		fit <- lm(Y ~ ., data=X[,1:npcr])
		#fit <- lm(Y ~ ., data=X)
		bmAIC <- stepAIC(fit, trace=FALSE, k=log(length(Y))) # for first 4, stepAIC has is lower than bmBIC <- stepAIC(fit, k=log(13))
		#bmAIC <- stepAIC(fit, trace=FALSE) # for first 4, stepAIC has is lower than bmBIC <- stepAIC(fit, k=log(13))

		zpredcon = predict(bmAIC, se.fit=TRUE)
		fpc[,j] <- zpredcon$fit


		zpred = predict(bmAIC, newdata, se.fit = TRUE)
		fpc.hol[1,j] <- zpred$fit

		fpc.hol.boot[,j] <- zpred$fit + rnorm(500,0,zpred$se.fit)

	}

	## Using first four and then the means for the remaining:
	pcmeans = matrix(colMeans(pcs[,-1:-npc]),nrow=1)
	fpc[,-1:-npc] = pcmeans[rep(1, nrow(pcs)),]
	fpc.hol[,-1:-npc] = pcmeans
	fpc.hol.boot[,-1:-npc] = apply(pcs[,-1:-npc],2,sample, 500, replace=TRUE)

	# Calibration measures (beta and R2)
	reconSST.sc = fpc %*% t(zsvd$u) 							# Contemporary - Yhat
	sumdiff = apply((data - reconSST.sc)^2,2,sum)
	sumYraw = apply(data^2,2,sum)
	beta = 1 - (sumdiff/sumYraw)

	reconSST.sc.hol = fpc.hol %*% t(zsvd$u)						# Paleo
	reconSST.ensemble = fpc.hol.boot %*% t(zsvd$u)				# Paleo
	ensemble.sd = apply(reconSST.ensemble, 2, sd)				# Uncertainty ensemble


	### Only using some of the PCs and then zero-ing out the rest:
	# pckeep = 1
	# fpc[,-1:-pckeep] = 0
	# fpc.hol[,-1:-pckeep] = 0
	# fpc.hol.boot[,-1:-pckeep] = 0

	# Ekeep = matrix(0,nrow=nsta,ncol=nsta)
	# Ekeep[,1:pckeep] = zsvd$u[,1:pckeep]

	#### Back transform using eigenvectors
	# reconSST.sc = fpc %*% t(Ekeep) 							# Contemporary
	# reconSST.sc.hol = fpc.hol %*% t(Ekeep)					# Paleo
	# reconSST.ensemble = fpc.hol.boot %*% t(Ekeep)			# Paleo
	# ensemble.sd = apply(reconSST.ensemble, 2, sd)			# Uncertainty ensemble

	#################################################################################################
	### 11. ## WRITE THE RECONSTRUCTION TO A TXT FILE FOR POST PROCESSING
	write.table(reconSST.sc.hol, 
		file=paste(paste("./postfiles/", region, "/PC", npc, "/", region, "_", proxy, "_", season, "_", npcr, "npcr", "_", yrBP, "ka.txt", sep=""),
		sep="\t"))
	# write.table(ensemble.sd, 
	# 	file=paste(paste("./postfiles/", region, "/PC", npc, "/", region, "_", proxy, "_", season, "_", npcr, "npcr", "_", yrBP, "ka_sd.txt", sep=""),
	# 	sep="\t"))

	# if (i == 11){
	# 	write.table(reconSST.sc, 
	# 	file=paste(paste("./postfiles/", region, "/PC", npc, "/", region, "_", proxy, "_", season, "_", npcr, "npcr", "_", yrBP, "ka_1901-2013.txt", sep=""),
	# 	sep="\t"))

	# 	write.table(data, 
 # 		file=paste(paste("./postfiles/", region, "/PC", npc, "/", region, "_", proxy, "_", season, "_", npcr, "npcr", "_actual_1901-2013.txt", sep=""),
 # 		sep="\t"))
	# }

}








# ###################################################################################################
# #### || FIGURE 5 || #### Multi-proxy reconstruction of SSTs
# validate = plot.validate2(var="sst", npc=as.character(npc),region=region, proxy=proxy, 
# 	season=season, index1, xgrid, ygrid, reconpt, zzpt, mat)


# # This one includes the scatter plots
# validate = plot.validate(npc=as.character(npc),region=region, proxy=proxy, 
# 	season=season, index1, xgrid, ygrid, reconpt, zzpt, mat)









###################################################################################################
############################# EXTRA CODES FOR VARIOUS THINGS ######################################
###################################################################################################

###################################################################################################
##### 11. ## PLOT RECONSTRUCTED HOLOCENE FIELDS USING IMAGE.PLOT
## Plot the Holocene reconstruction
dev.new(width=10, height=3)
par(mfrow = c(1, 1), mar = c(4, 4, 4, 1))
zfull = rep(NaN,(nx*ny))
zfull[index1] = reconSST.sc.hol 						# add in [index1] if youre plotting India rainfall or SST
zmat = matrix(zfull,nrow=nx,ncol=ny)
image.plot(xgrid,ygrid,zmat,
	ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)),
	xlab="",ylab="", main=paste(yrBP, "ka", sep=""),
	cex.axis=1.1, col=myPalette1(100),
	zlim=c(-2,2), horizontal=TRUE
	)
contour(xgrid, ygrid, zmat, add = TRUE, nlev = 6, lwd = 2)
mapnames <- map("world2", 
	xlim=c(min(ygrid), max(xgrid)), ylim=c(min(ygrid), max(ygrid)), 
	boundary=TRUE, interior=TRUE, col="black", add=TRUE
)


## Plot the SD map
dev.new(width=10, height=3)
par(mfrow = c(1, 1), mar = c(4, 4, 4, 1))
zfull = rep(NaN,(nx*ny))
zfull[index1] = ensemble.sd 						# add in [index1] if youre plotting India rainfall or SST
zmat = matrix(zfull,nrow=nx,ncol=ny)
image(xgrid,ygrid,zmat,
	ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)),
	xlab="",ylab="", main="",
	yaxt="n", col=rev(myPalette5(100)),
	zlim=c(0,1)
	)
contour(xgrid, ygrid, zmat, add = TRUE, levels = seq(0,1,0.1), lwd = 2, labcex=1.25)
mapnames <- map("world2", 
	xlim=c(min(ygrid), max(xgrid)), ylim=c(min(ygrid), max(ygrid)), 
	boundary=TRUE, interior=TRUE, col="black", add=TRUE
)
axis(2, at = c(-10,0,10))



## Plot the variance
data.var = apply(data,2,var)
dev.new(width=10, height=3)
par(mfrow = c(1, 1), mar = c(4, 4, 4, 1))
zfull = rep(NaN,(nx*ny))
zfull[index1] = data.var 						# add in [index1] if youre plotting India rainfall or SST
zmat = matrix(zfull,nrow=nx,ncol=ny)
image(xgrid,ygrid,zmat,
	ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)),
	xlab="",ylab="", main="",
	yaxt="n", col=rev(myPalette5(100)),
	zlim=c(0,1)
	)
contour(xgrid, ygrid, zmat, add = TRUE, levels = seq(0,1,0.1), lwd = 2, labcex=1.25)
mapnames <- map("world2", 
	xlim=c(min(ygrid), max(xgrid)), ylim=c(min(ygrid), max(ygrid)), 
	boundary=TRUE, interior=TRUE, col="black", add=TRUE
)
axis(2, at = c(-10,0,10))




#### PLOT COMPARISON YEARS
reconwind = read.table("./postfiles/EqPac/CCA-PC6/EqPac_both_annual_6npc_10ka_1949-2013.txt", header=TRUE)
actualwind = read.table("./postfiles/EqPac/CCA-PC6/EqPac_both_annual_6npc_actual_1949-2013.txt", header=TRUE)
wxgrid <- seq(100,300,by=2.5)
nwx = length(wxgrid)
wygrid <- seq(-10,10,by=2.5)
nwy = length(wygrid)

yr = 1997
iix = which(1854:2013 == yr)
iix2 = which(1949:2013 == yr)
year = "1997-1998"
max = 4
min = -4

dev.new(width=10, height=4)
par(mfrow = c(2, 1), mar = c(2,2,2,2))
# Actual
zfull = rep(NaN,(nx*ny))
zfull[index1] = data[iix,]					# add in [index1] if youre plotting India rainfall or SST
zmat = matrix(zfull,nrow=nx,ncol=ny)
image(xgrid,ygrid,zmat,
	ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)),
	xlab="",ylab="", main=paste("Actual", year),
	cex.axis=1.1, yaxt="n", col=myPalette1(100),
	zlim=c(min,max)
	)
axis(2, at = c(-10,0,10))
mapnames <- map("world2", 
	xlim=c(min(ygrid), max(xgrid)), ylim=c(min(ygrid), max(ygrid)), 
	boundary=TRUE, interior=TRUE, add=TRUE, lwd=1.5, col="grey50"
)
wfull = as.numeric(actualwind[iix2,])
wmat = matrix(wfull,nrow=nwx,ncol=nwy)
contour(wxgrid, wygrid, wmat, add=TRUE, levels = seq(-6.5,6.5,1), lwd = 2, labcex=1.25)

# Reconstructed
zfull = rep(NaN,(nx*ny))
zfull[index1] = reconSST.sc[iix,] 						# add in [index1] if youre plotting India rainfall or SST
zmat = matrix(zfull,nrow=nx,ncol=ny)
image(xgrid,ygrid,zmat,
	ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)),
	xlab="",ylab="", main=paste("Reconstructed", year),
	cex.axis=1.1, yaxt="n", col=myPalette1(100),
	zlim=c(min,max)
	)
axis(2, at = c(-10,0,10))
mapnames <- map("world2", 
	xlim=c(min(ygrid), max(xgrid)), ylim=c(min(ygrid), max(ygrid)), 
	boundary=TRUE, interior=TRUE, add=TRUE, lwd=1.5, col="grey50"
)
#contour(xgrid, ygrid, zmat, add = TRUE, levels = seq(min,max,0.5), lwd = 2, labcex=1.25)
wfull = as.numeric(reconwind[iix2,])
wmat = matrix(wfull,nrow=nwx,ncol=nwy)
contour(wxgrid, wygrid, wmat, add=TRUE, levels = seq(-6.5,6.5,1), lwd = 2, labcex=1.25)

# Average for contemporary period 1901-2013
present = apply(reconSST.sc[81:110,],2,mean)

dev.new(width=10, height=4)
par(mfrow = c(2, 1), mar = c(2,2,2,2))
# Actual
zfull = rep(NaN,(nx*ny))
zfull[index1] = present					# add in [index1] if youre plotting India rainfall or SST
zmat = matrix(zfull,nrow=nx,ncol=ny)
image(xgrid,ygrid,zmat,
	ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)),
	xlab="",ylab="",
	cex.axis=1.1, yaxt="n", col=myPalette1(100),
	zlim=c(-1.5,1.5)
	)
axis(2, at = c(-10,0,10))
mapnames <- map("world2", 
	xlim=c(min(ygrid), max(xgrid)), ylim=c(min(ygrid), max(ygrid)), 
	boundary=TRUE, interior=TRUE, add=TRUE, lwd=1.5, col="grey50"
)
contour(xgrid, ygrid, zmat, add = TRUE, levels = seq(-1,1,0.2), lwd = 2, labcex=1.25)



###################################################################################################
##### 13. ## TEST MODEL: PLOT SPATIAL R2 OF MODEL
cor <- c()
for (i in 1:nsta){
	xcor <- cor(reconSST.sc[,i],data[,i])^2
	cor <- c(cor,xcor)
}


dev.new(width=10, height=4)
par(mfrow = c(2, 1), mar = c(2,2,2,2))
# BETA
zfull = rep(NaN,(nx*ny))
zfull[index1] = beta
zmat = matrix(zfull,nrow=nx,ncol=ny)
image(xgrid,ygrid,zmat,
	ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)),
	xlab="",ylab="", main=expression(paste("SST Calibration ", beta)),
	cex.axis=1.1, yaxt="n", col=rev(myPalette3(100)),
	zlim=c(-1,1)
	)
axis(2, at = c(-10,0,10))
mapnames <- map("world2", 
	xlim=c(min(ygrid), max(xgrid)), ylim=c(min(ygrid), max(ygrid)), 
	boundary=TRUE, interior=TRUE, add=TRUE, lwd=1.5, col="grey50"
)
contour(xgrid, ygrid, zmat, add = TRUE, levels = seq(-1,1,0.05), lwd = 2, labcex=1.25)

# R^2
zfull = rep(NaN,(nx*ny))
zfull[index1] = cor
zmat = matrix(zfull,nrow=nx,ncol=ny)
image(xgrid,ygrid,zmat,
	ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)),
	xlab="",ylab="", main=expression(paste("SST Calibration ", R^{2})),
	cex.axis=1.1, yaxt="n", col=rev(myPalette2(100)),
	zlim=c(0,1)
	)
axis(2, at = c(-10,0,10))
mapnames <- map("world2", 
	xlim=c(min(ygrid), max(xgrid)), ylim=c(min(ygrid), max(ygrid)), 
	boundary=TRUE, interior=TRUE, add=TRUE, lwd=1.5, col="grey50"
)
contour(xgrid, ygrid, zmat, add = TRUE, levels = seq(0,1,0.1), lwd = 2, labcex=1.25)





####################### PLOT THE POINT MAPS SHOWING HOW WELL THE PROXY RECORDS WERE FIT
years = seq(0,10,1)	
reconSST.mat = c()
reconSST.mat.sd = c()
	
for (i in 1:length(years)){
		test = read.table(paste("./post-files/", region, "/PC", npc, "/", region, "_", proxy, "_", season, "_", npcr, "npcr", "_", i-1, "ka.txt", sep=""),header=TRUE)
		reconSST.mat = rbind(reconSST.mat,test)

		sd = read.table(paste("./post-files/", region, "/PC", npc, "/", region, "_", proxy, "_", season, "_", npcr, "npcr", "_", i-1, "ka_sd.txt", sep=""),header=TRUE)
		reconSST.mat.sd = rbind(reconSST.mat.sd,t(sd))
}

reconSST.pt = reconSST.mat[,zzpt]
reconSST.pt.sd = reconSST.mat.sd[,zzpt]
labels = c(22,23,24,25,27,28,11,12,13,14,15,16,17,18,19,20,1,2,3,4,5,6,7,8,9,10)

### E. Pacific
dev.new(width=8,height=8)
par(mfrow = c(4,3), mar = c(1, 1, 1, 1), oma = c(2, 2, 1, 1))		
for (i in 1:11){
	actual = mat[,i]
	reconsst = reconSST.pt[,i]
	sd = reconSST.pt.sd[,i]
	val = which(is.na(actual) == FALSE)
	# miny = floor(min(actual[val],reconsst[val]-sd[val]))
	# maxy = ceiling(max(actual[val],reconsst[val]+sd[val]))
	miny = -3.5
	maxy = 2

	xpoly = c(c(0:10)[val],rev(c(0:10)[val]))
	ypoly = c(actual[val]-1,rev(actual[val]+1))
	plot(0:10, actual, type="l", xlim = rev(range(0:10)), ylim = c(miny,maxy), axes=FALSE, ylab="", xlab="",col="#0072B2")
	polygon(xpoly,ypoly,col=rgb(0.337, 0.706, 0.914 ,0.25),border = "#56B4E9")

	abline(v=c(0:10),col="grey90")
	abline(h=0,col="grey90")

	points(0:10,actual,pch=20, col="#0072B2")
	points(0:10,reconsst,pch=23,bg="white",col="#D55E00")
	arrows(0:10,reconsst+sd,0:10,reconsst-sd,angle=90,code=3,length=0.05, lwd=0.5,col="#E69F00")

	text(8,1.5,paste("proxy no.",labels[i]))
	text(3,1.5,paste("RSS =",round(sum((actual[val] - reconsst[val])^2),2)))
	text(8,-3, paste("EOF1 =",round(zsvd.p$u[i,1],2)))
	text(3,-3, paste("EOF2 =",round(zsvd.p$u[i,2],2)))

	if (is.element(i,c(1,4,7,10)) == TRUE){
		axis(2, at = c(-3,-1.5,0,1.5))
	} else {
		axis(2, at = c(-3,-1.5,0,1.5), labels = c("","","",""))
	}

	if (is.element(i,c(9,10,11)) == TRUE){
		axis(1, at = c(0:10))
	} else {
		axis(1, at = c(0:10), labels = c("","","","","","","","","","",""))
	}
}


### W. Pacific
dev.new(width=8,height=10)
par(mfrow = c(5,3), mar = c(1, 1, 1, 1), oma = c(2, 2, 1, 1))		
for (i in 12:26){
	actual = mat[,i]
	reconsst = reconSST.pt[,i]
	sd = reconSST.pt.sd[,i]
	val = which(is.na(actual) == FALSE)
	#miny = floor(min(actual[val],reconsst[val]-sd[val]))
	#maxy = ceiling(max(actual[val],reconsst[val]+sd[val]))
	miny = -3.5
	maxy = 2

	xpoly = c(c(0:10)[val],rev(c(0:10)[val]))
	ypoly = c(actual[val]-1,rev(actual[val]+1))
	plot(0:10, actual, type="l", xlim = rev(range(0:10)), axes=FALSE, ylim = c(miny,maxy), ylab="", xlab="",col="#0072B2")
	polygon(xpoly,ypoly,col=rgb(0.337, 0.706, 0.914 ,0.25),border = "#56B4E9")

	abline(v=c(0:10),col="grey90")
	abline(h=0,col="grey90")

	points(0:10,actual,pch=20, col="#0072B2")
	points(0:10,reconsst,pch=23,bg="white",col="#D55E00")
	arrows(0:10,reconsst+sd,0:10,reconsst-sd,angle=90,code=3,length=0.05, lwd=0.5,col="#E69F00")

	text(8,1.5,paste("proxy no.",labels[i]))
	text(3,1.5,paste("RSS =",round(sum((actual[val] - reconsst[val])^2),2)))
	text(8,-3, paste("EOF1 =",round(zsvd.p$u[i,1],2)))
	text(3,-3, paste("EOF2 =",round(zsvd.p$u[i,2],2)))

	if (is.element(i,c(12,15,18,21,24)) == TRUE){
		axis(2, at = c(-3,-1.5,0,1.5))
	} else {
		axis(2, at = c(-3,-1.5,0,1.5), labels = c("","","",""))
	}

	if (is.element(i,c(24,25,26)) == TRUE){
		axis(1, at = c(0:10))
	} else {
		axis(1, at = c(0:10), labels = c("","","","","","","","","","",""))
	}
}



###################################################################################################
### AVERAGE VALUES OVER NINO INDICES
###################################################################################################
zfull = rep(NaN,(nx*ny))
zfull[index1] = reconSST.sc.hol 						# add in [index1] if youre plotting India rainfall or SST
zmat = matrix(zfull,nrow=nx,ncol=ny)

sdzfull = rep(NaN,(nx*ny))
sdzfull[index1] = ensemble.sd 						# add in [index1] if youre plotting India rainfall or SST
sdzmat = matrix(sdzfull,nrow=nx,ncol=ny)

lat1 = which(ygrid == -4)
lat2 = which(ygrid == 4)


# WPAC 120 - 160
lon1 = which(xgrid == 120)
lon2 = which(xgrid == 160)
wpac = mean(zmat[lon1:lon2,lat1:lat2], na.rm=TRUE)
nas = dim(which(sdzmat[lon1:lon2,lat1:lat2] == "NaN",arr.ind=TRUE))[1]
n = length(sdzmat[lon1:lon2,lat1:lat2]) - nas
wpacsd = sqrt(sum(sdzmat[lon1:lon2,lat1:lat2]^2, na.rm=TRUE))/sqrt(n)
round(wpac,2)
round(wpacsd,2)


# NINO4 160E - 150W, 5S - 5N = 160 - 210
lon1 = which(xgrid == 160)
lon2 = which(xgrid == 210)
nino4 = mean(zmat[lon1:lon2,lat1:lat2])
n = length(sdzmat[lon1:lon2,lat1:lat2])
nino4sd = sqrt(sum(sdzmat[lon1:lon2,lat1:lat2]^2, na.rm=TRUE))/sqrt(n)
round(nino4,2)
round(nino4sd,2)

# NINO34 120W - 170W, 5S - 5N = 190 - 240
lon1 = which(xgrid == 190)
lon2 = which(xgrid == 240)
nino34 = mean(zmat[lon1:lon2,lat1:lat2])
n = length(sdzmat[lon1:lon2,lat1:lat2])
nino34sd = sqrt(sum(sdzmat[lon1:lon2,lat1:lat2]^2, na.rm=TRUE))/sqrt(n)
round(nino34,2)
round(nino34sd,2)

# NINO3 90W - 150W, 5S - 5N = 210 - 270
lon1 = which(xgrid == 210)
lon2 = which(xgrid == 270)
nino3 = mean(zmat[lon1:lon2,lat1:lat2])
nino3sd = mean(sdzmat[lon1:lon2,lat1:lat2])
round(nino3,2)
round(nino3sd,2)

# NINO12 90W - 80W, 10S - 0 = 270 - 280
lat1 = which(ygrid == -10)
lat2 = which(ygrid == 0)
lon1 = which(xgrid == 270)
lon2 = which(xgrid == 280)
nino12 = mean(zmat[lon1:lon2,lat1:lat2], na.rm=TRUE)
nas = dim(which(sdzmat[lon1:lon2,lat1:lat2] == "NaN",arr.ind=TRUE))[1]
n = length(sdzmat[lon1:lon2,lat1:lat2]) - nas
nino12sd = sqrt(sum(sdzmat[lon1:lon2,lat1:lat2]^2, na.rm=TRUE))/sqrt(n)
round(nino12,2)
round(nino12sd,2)

# TNI = NINO12 - NINO4
tni = nino12 - nino4
round(tni,2)
round(sqrt(nino12sd^2 + nino4sd^2),2)

# WTNI = NINO12 - WPAC
wtni = nino12 - wpac
round(wtni,2)
round(sqrt(nino12sd^2 + wpacsd^2),2)

#### ATTEMPT 4 AT PLOTTING EQUATORIAL PACIFIC
library("PBSmapping")
library("ggplot2")
library("maps")
library("data.table")
library("dplyr")

# plot limits
xlim = c(90,330)
ylim = c(-30,30)

clip_map_polygons <- function(mapname, xlim, ylim, run_tol=300){
    require('dplyr')
    require('data.table')

    map = map_data(mapname)

    m = data.table(map_data(mapname))
    d = data.table(summarise(group_by(m,region),run=max(diff(long))))
    region_keep = d[!(run>run_tol)]$region
    map = as.data.frame(m[region %in% region_keep])
    #class(map) = 'data.frame'

    setnames(map, c('X','Y','PID','POS','region','subregion'))
    map = clipPolys(map, xlim=xlim,ylim=ylim, keepExtra=TRUE)
    map = as.data.table(as.data.frame(map))
    setnames(map, c('X','Y', 'PID'), c('lon', 'lat', 'group'))
    map
}

worldmap = clip_map_polygons('world2',xlim,ylim)

p = ggplot()+

	coord_map(xlim=xlim,ylim=ylim) +
	geom_polygon(data=worldmap, aes(lon,lat,group=group),fill='grey80',color='black',size=0.15)+
	geom_rect(aes(xmin=120, xmax=160, ymin=-4, ymax=4), linetype=3, color="black", fill="#D55E00", alpha=0.25)+
	geom_rect(aes(xmin=160, xmax=210, ymin=-4, ymax=4), linetype=3, color="black", fill="#9900CC", alpha=0.25)+
	geom_rect(aes(xmin=210, xmax=270, ymin=-4, ymax=4), linetype=3, color="black", fill="#0072B2", alpha=0.25)+
	geom_rect(aes(xmin=190, xmax=240, ymin=-4, ymax=4), alpha=0.5)+
	geom_rect(aes(xmin=270, xmax=280, ymin=-10, ymax=0), linetype=3, color="black", fill="#009E73", alpha=0.25)+

	#geom_point(data=data, aes(x=data$lon360, y=data$lat, color=region))+
	annotate("text", x = 140, y = 6.5, label = "WPAC", color="#D55E00", size=3, fontface="bold")+
	annotate("text", x = 185, y = 6.5, label = "NINO4", color="#9900CC", size=3, fontface="bold")+
	annotate("text", x = 215, y = -6.5, label = "NINO3.4", color="grey30", size=3, fontface="bold")+
	annotate("text", x = 240, y = 6.5, label = "NINO3", color="#0072B2", size=3, fontface="bold")+
	annotate("text", x = 272, y = -13.25, label = "NINO1+2", color="#009E73", size=3, fontface="bold")+
	theme_bw()+
	labs(x='',y='')
print(p)


















###################################################################################################
##### 14. ## PLOTTING IN GGPLOT!!
library("PBSmapping")
library("ggplot2")
library("maps")
library("data.table")
library("dplyr")
library("directlabels")

# plot limits
xlim = c(100,300)
ylim = c(-10,10)

clip_map_polygons <- function(mapname, xlim, ylim, run_tol=300){
    require('dplyr')
    require('data.table')

    map = map_data(mapname)

    m = data.table(map_data(mapname))
    d = data.table(summarise(group_by(m,region),run=max(diff(long))))
    region_keep = d[!(run>run_tol)]$region
    map = as.data.frame(m[region %in% region_keep])
    #class(map) = 'data.frame'

    setnames(map, c('X','Y','PID','POS','region','subregion'))
    map = clipPolys(map, xlim=xlim,ylim=ylim, keepExtra=TRUE)
    map = as.data.table(as.data.frame(map))
    setnames(map, c('X','Y', 'PID'), c('lon', 'lat', 'group'))
    map
}

worldmap = clip_map_polygons('world2',xlim,ylim)

long = rep(xgrid,ny)
lat = c()
for (i in 1:ny){
	ylat = rep(ygrid[i],nx)
	lat = c(lat,ylat)
}
SST = rep(NA,(ny*nx))
SST[index1] = cor
plot.data = data.frame(cbind(long,lat,SST))


p = ggplot()+
	geom_tile(aes(long,lat,fill=SST),data=na.omit(plot.data))+
    scale_fill_gradient2(expression(beta),low = "#FFFFCC", mid = "#FD8C3B", high = "#7F0026",limits=c(0,1)) +
    #scale_fill_gradient2(expression(R^2),low = "#F7FCF4", mid = "#73C476", high = "#00441A",limits=c(0,1)) +
    coord_map(xlim=xlim,ylim=ylim) +
    geom_polygon(data=worldmap,aes(lon,lat,group=group),fill='black',color='grey75',size=0.25) +
    theme_bw() +
    labs(x='',y='')
#p1 = p + stat_contour(aes(long,lat,z=SST), colour = "grey25", breaks = seq(0,1,0.1), data = plot.data, size = 0.25)

dev.new(width=10, height=3.5)
p
contour(xgrid, ygrid, zmat, add = TRUE, levels = seq(0,1,0.1), lwd = 2, labcex=1.5)


#Color Palettes:
#Blue - White - Red = royalblue3, white, tomato3
#Yellow - Orange - Red = #7F0026, #FD8C3B, #FFFFCC
#Greens = "#00441A" "#73C476" "#F7FCF4"


### Messing around with contour labels
# p2 <- p + 
# 	scale_color_gradient(low = "#000000", high = "#000000"))

# p2 = p + stat_contour(aes(long,lat,z=SST,colour=..level..),data=plot.data)
# direct.label((p2 + scale_color_gradient(low = "#000000", high = "#000000")))

# p = ggplot()+
# 	geom_tile(aes(long,lat,fill=SST),data=na.omit(plot.data))+
#     scale_fill_gradient2('SST',low = "royalblue3", mid = "white", high = "tomato3") +







