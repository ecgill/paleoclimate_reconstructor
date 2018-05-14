rm(list=ls())
setwd("/Users/emilygill/Documents/University of Colorado/PHD Research/6. SST Reconstruction/REDO Sep2015/")

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


###################################################################################################
##### 2. ## FULL SLAB SST DATA AND CLIMATOLOGY. SCALE FULLSLAB DATA.
results = read.contemp(region = region, season = "annual") 	   
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

#write(index, file="./postfiles/ISM/diagnostics/index-rain.txt")


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



###################################################################################################
##### 7. ## PCA ON RAINFALL (seasonal total)
# means  = apply(precip.sestot,2,mean)
# a.precip = sweep(precip.sestot,2,means)
# sc.precip = a.precip

sc.precip = scale(precip.davg)
contemp.sd = attr(sc.precip,'scaled:scale')
contemp.mu = attr(sc.precip,'scaled:center')

zs.rain <- var(sc.precip)
zsvd.rain <- svd(zs.rain)
pcs.rain <- t(t(zsvd.rain$u) %*% t(sc.precip))
lambdas.rain <- (zsvd.rain$d/sum(zsvd.rain$d))

# write.table(zsvd.rain$u, 
# 	file=paste(paste("./postfiles/ISM/diagnostics/", region, "_", proxy, "_summer_zvsd.raindavg$u", sep=""),
# 	sep="\t"))
# write.table(pcs.rain, 
# 	file=paste(paste("./postfiles/ISM/diagnostics/", region, "_", proxy, "_summer_pcs.raindavg", sep=""),
# 	sep="\t"))
# write(lambdas.rain, 
# 	file=paste(paste("./postfiles/ISM/diagnostics/", region, "_", proxy, "_summer_lambdas.raindavg", sep=""),
# 	sep="\t"))



###################################################################################################
##### 8. ## DEFINE REGIONS OF INDIA
zfull = rep(NaN,(nx*ny))
zfull[index] = contemp.mu
zmat = matrix(zfull,nrow=nx,ncol=ny)
image(xgrid,ygrid,zmat,
	ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),99), axes=FALSE, zlim=c(0,30), 
	ylab="",xlab="",col=rev(myPalette3(100))
	)
mapnames <- map("world2", 
	xlim=c(min(xgrid), max(xgrid)), ylim=c(min(ygrid), max(ygrid)), 
	boundary=TRUE, interior=TRUE, col="grey50", add=TRUE
)

library(seqinr)
colors = c("#D55E00","#E69F00",  "#0072B2", "#009E73", "black", "#9900CC")
col1 <- col2alpha("#D55E00",alpha=0.1)
col2 <- col2alpha("#E69F00",alpha=0.1)
col3 <- col2alpha("#0072B2",alpha=0.1)
col4 <- col2alpha("#009E73",alpha=0.1)
col5 <- col2alpha("black",alpha=0.1)
col6 <- col2alpha("#9900CC",alpha=0.1)
acolors = c(col1,col2,col3,col4,col5,col6)


ii = 1:(nx*ny)
imat = matrix(ii,nx,ny)

# W. Ghats: Lat 6.5 - 20.5, Lon 66.5 - 76.5
lat = c(20.5,20.5,19.5,19.5,18.5,18.5,17.5,16.5,15.5,15.5,14.5,13.5,13.5,12.5,12.5,11.5,11.5,10.5,10.5,9.5)
lon = c(72.5,73.5,72.5,73.5,72.5,73.5,73.5,73.5,73.5,74.5,74.5,74.5,75.5,74.5,75.5,75.5,76.5,75.5,76.5,76.5)
#rect(min(lon)-0.5,min(lat)-0.5,max(lon)+0.5,max(lat)+0.5,border=colors[4],lwd=3)
ista = c()
for (m in 1:length(lat)){
	points(lon[m],lat[m],pch=19,cex=0.5,col=colors[1])
	ii = imat[which(lon[m] == xgrid), which(lat[m] == ygrid)]
	ista = c(ista,which(rain$ind == ii))
}
wg.davg = apply(precip.davg[,ista],1,mean)
wg.sestot = apply(precip.sestot[,ista],1,mean)



# Indo-Gangetic Plains: Lat 24.5 - 30.5, Lon 79.5 - 89.5
lat = c(30.5,30.5,30.5,
	29.5,29.5,29.5,
	28.5,28.5,28.5,
	27.5,27.5,27.5,27.5,27.5,
	26.5,26.5,26.5,26.5,26.5,26.5,26.5,
	25.5,25.5,25.5,25.5,25.5,25.5)
lon = c(78.5,79.5,80.5,
	78.5,79.5,80.5,
	79.5,80.5,81.5,
	80.5,81.5,82.5,83.5,84.5,
	81.5,82.5,83.5,84.5,85.5,86.5,87.5,
	82.5,83.5,84.5,85.5,86.5,87.5)
#rect(min(lon)-0.5,min(lat)-0.5,max(lon)+0.5,max(lat)+0.5,border=colors[2],lwd=3)
ista = c()
for (m in 1:length(lat)){
	points(lon[m],lat[m],pch=19,cex=0.5,col=colors[2])
	ii = imat[which(lon[m] == xgrid), which(lat[m] == ygrid)]
	ista = c(ista,which(rain$ind == ii))
}
igprain.davg = apply(precip.davg[,ista],1,mean)
igprain.sestot = apply(precip.sestot[,ista],1,mean)



# Godavari: Lat 17.5 - 22.5, Lon 77.5 - 85.5
lat = c(22.5,22.5,22.5,22.5,22.5,22.5,
	21.5,21.5,21.5,21.5,21.5,21.5,
	20.5,20.5,20.5,20.5,20.5,20.5,
	19.5,19.5,19.5,19.5,19.5,19.5,
	18.5,18.5,18.5,18.5,18.5,
	17.5,17.5,17.5,
	16.5,16.5)
lon = c(78.5,79.5,80.5,81.5,82.5,83.5,
	78.5,79.5,80.5,81.5,82.5,83.5,
	78.5,79.5,80.5,81.5,82.5,83.5,
	78.5,79.5,80.5,81.5,82.5,83.5,
	78.5,79.5,80.5,81.5,82.5,
	80.5,81.5,82.5,
	81.5,82.5)
#rect(min(lon)-0.5,min(lat)-0.5,max(lon)+0.5,max(lat)+0.5,border=colors[3],lwd=3)
ista = c()
for (m in 1:length(lat)){
	points(lon[m],lat[m],pch=19,cex=0.5,col=colors[3])
	ii = imat[which(lon[m] == xgrid), which(lat[m] == ygrid)]
	ista = c(ista,which(rain$ind == ii))
}
godrain.davg = apply(precip.davg[,ista],1,mean)
godrain.sestot = apply(precip.sestot[,ista],1,mean)



# NW Rajasthan: Lat 24.5 - 30.5, Lon 66.5 - 75.5
lat = c(30.5,30.5,
	29.5,29.5,29.5,29.5,
	28.5,28.5,28.5,28.5,28.5,28.5,
	27.5,27.5,27.5,27.5,27.5,27.5,27.5,
	26.5,26.5,26.5,26.5,26.5,26.5,
	25.5,25.5,25.5,25.5,25.5,25.5,
	24.5,24.5,24.5,24.5,24.5,24.5)
lon = c(73.5,74.5,
	72.5,73.4,74.4,75.5,
	70.5,71.5,72.5,73.5,74.4,75.5,
	69.5,70.5,71.5,72.5,73.5,74.4,75.5,
	70.5,71.5,72.5,73.5,74.4,75.5,
	70.5,71.5,72.5,73.5,74.4,75.5,
	70.5,71.5,72.5,73.5,74.4,75.5)
#rect(min(lon)-0.5,min(lat)-0.5,max(lon)+0.5,max(lat)+0.5,border=colors[1],lwd=3)
ista = c()
for (m in 1:length(lat)){
	points(lon[m],lat[m],pch=19,cex=0.5,col=colors[4])
	ii = imat[which(lon[m] == xgrid), which(lat[m] == ygrid)]
	ista = c(ista,which(rain$ind == ii))
}
rajrain.davg = apply(precip.davg[,ista],1,mean)
rajrain.sestot = apply(precip.sestot[,ista],1,mean)



# NE India: Lat 17.5 - 22.5, Lon 77.5 - 85.5
lat = c(29.5,29.5,29.5,28.5,28.5,28.5,28.5,28.5,27.5,27.5,27.5,27.5,26.5)
lon = c(94.5,95.5,96.5,93.5,94.5,95.5,96.5,97.5,94.5,95.5,96.5,97.5,95.5)
#rect(min(lon)-0.5,min(lat)-0.5,max(lon)+0.5,max(lat)+0.5,border=colors[5],lwd=3)
ista = c()
for (m in 1:length(lat)){
	points(lon[m],lat[m],pch=19,cex=0.5,col=colors[5])
	ii = imat[which(lon[m] == xgrid), which(lat[m] == ygrid)]
	ista = c(ista,which(rain$ind == ii))
}
eb.davg = apply(precip.davg[,ista],1,mean)
eb.sestot = apply(precip.sestot[,ista],1,mean)



# Northern India
lat = c(36.5,36.5,36.5,36.5,36.5,36.5,
	35.5,35.5,35.5,35.5,35.5,35.5,35.5,35.5,
	34.5,34.5,34.5,34.5,34.5,34.5,34.5,
	33.5,33.5,33.5,33.5,33.5,33.5,33.5,
	32.5,32.5,32.5,32.5,32.5,32.5)
lon = c(72.5,73.5,74.5,75.5,78.6,79.5,
	73.5,74.5,75.5,76.5,77.5,78.5,79.5,80.5,
	73.5,74.5,75.5,76.5,77.5,78.5,79.5,
	73.5,74.5,75.5,76.5,77.5,78.5,79.5,
	74.5,75.5,76.5,77.5,78.5,79.5)
#rect(min(lon)-0.5,min(lat)-0.5,max(lon)+0.5,max(lat)+0.5,border=colors[6],lwd=3)
ista = c()
for (m in 1:length(lat)){
	points(lon[m],lat[m],pch=19,cex=0.5,col=colors[6])
	ii = imat[which(lon[m] == xgrid), which(lat[m] == ygrid)]
	ista = c(ista,which(rain$ind == ii))
}
north.davg = apply(precip.davg[,ista],1,mean)
north.sestot = apply(precip.sestot[,ista],1,mean) 

# Bind them together
Y.davg = unname(cbind(wg.davg, igprain.davg, godrain.davg, rajrain.davg, eb.davg, north.davg))
Y.sestot = unname(cbind(wg.sestot, igprain.sestot, godrain.sestot, rajrain.sestot, eb.sestot, north.sestot))





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

zpred.davg = matrix(NA,11,6)
zpred.davg.se = matrix(NA,11,6)
zpred.sestot = matrix(NA,11,6)
zpred.sestot.se = matrix(NA,11,6)
for (i in 1:length(yrsBP)){
	yrBP = yrsBP[i]

	missing = which(is.na(mat[yrBP+1,]))
	if (length(missing) > 0){
		var = seq(1,nreconpt,1)[-missing]
	} else {
		var = seq(1,nreconpt,1)
	}

	istart = which(1854:2013 == 1901)
	istop = which(1854:2013 == 2004)
	zreconsst = reconsst[istart:istop,var]

	zs.p <- var(zreconsst)
	zsvd.p <- svd(zs.p)
	pcs.p <- t(t(zsvd.p$u) %*% t(zreconsst))
	lambdas.p <- (zsvd.p$d/sum(zsvd.p$d))

	# write.table(zsvd.p$u, 
	# 	file=paste(paste("./postfiles/ISM/diagnostics/", region, "_", proxy, "_summer_zvsd.p$u", sep=""),
	# 	sep="\t"))
	# write.table(pcs.p, 
	# 	file=paste(paste("./postfiles/ISM/diagnostics/", region, "_", proxy, "_summer_pcs.p", sep=""),
	# 	sep="\t"))
	# write(lambdas.p, 
	# 	file=paste(paste("./postfiles/ISM/diagnostics/", region, "_", proxy, "_summer_lambdas.p", sep=""),
	# 	sep="\t"))

	###################################################################################################
	##### 9. ## PREDICT RAINFALL TIMESERIES
	## TO FIT USING THE PCS
	xx.hol = matrix(mat[yrBP+1,var],1,length(mat[yrBP+1,var]))
	newdata = as.data.frame(xx.hol %*% zsvd.p$u)

	# newdata = as.data.frame(mat %*% zsvd.p$u)
	#  X = as.data.frame(pcs.p)

	X = as.data.frame(pcs.p[,1:length(mat[yrBP+1,var])])
	npcr = 4
	for (j in 1:6){
		## Daily monsoon average
		Y1 = Y.davg[,j]
		fit = lm(Y1 ~ ., data=X[,1:npcr])
		#X = as.data.frame(X[,1])
		#fit = lm(Y1 ~ X[,1])
		#fit = lm(Y1 ~ X$V1)
		#bmAIC = stepAIC(fit, trace=FALSE, k=log(length(Y1)))
		#zz = predict(bmAIC, newdata, se.fit=TRUE)
		zz = predict(fit, newdata,se.fit=TRUE)
		zpred.davg[i,j] = zz$fit
		zpred.davg.se[i,j] = zz$se
		
		## Seasonal total
		Y2 = Y.sestot[,j]
		fit = lm(Y2 ~ ., data=X[,1:npcr])
		#fit = lm(Y2 ~ X[,1:npcr])
		#bmAIC = stepAIC(fit, trace=FALSE, k=log(length(Y2)))
		#zz = predict(bmAIC, newdata, se.fit=TRUE)
		zz = predict(fit,newdata,se.fit=TRUE)
		zpred.sestot[i,j] = zz$fit
		zpred.sestot.se[i,j] = zz$se
	}




	###################################################################################################
	##### 9. ## CCA TO PREDICT RAINFALL SPATIALLY
	###################################################################################################
	##### 7. ## PERFORM CCA to reconstruct U-WIND
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
	xx.hol = matrix(mat[yrBP+1,var],1,length(mat[yrBP+1,var]))
	newdata = xx.hol %*% zsvd.p$u[,1:ncca]

	# ndx = xx.hol %*% zsvd.p$u
	# write.table(ndx, 
	# file=paste(paste("./postfiles/ISM/diagnostics/EqPac_both_summer_newdata10ka", sep=""),
	# sep="\t"))

	### Predict past wind PC field.
	pcmeans = matrix(colMeans(pcs.rain[,-1:-ncca]),nrow=1)	
	fpc.hol = matrix(NA,1,ncol(pcs.rain))				# Set up matrix for reconstructed Holocene rain
	scaled = newdata %*% AA %*% betahat					# Scaled reconstructed Hol PCs
	unscaled = (scaled * Ysd) + Ymean					# Unscaled reconstructed Hol PCs (sd and mean of PCs have been added)
	fpc.hol[1,1:ncca] = unscaled
	fpc.hol[1,-1:-ncca] = pcmeans						# Add in means of remaning PCs
	recon.hol = fpc.hol %*% t(zsvd.rain$u)				# Scaled reconstructed Hol rain values (sd and mean may need to be added)
	recon.hol.unsc = t(t(recon.hol) * contemp.sd + contemp.mu)

	### Predict contemporary wind PC field.
	fpc = matrix(NA,nrow(pcs.rain),ncol(pcs.rain))		# Set up matrix for reconstructed contemporary rain
	scaled = X1 %*% AA %*% betahat						# Scaled reconstructed contemporary PCs
	unscaled = t(t(scaled) * Ysd + Ymean)				# Unscaled reconstructed contemporary PCs (sd and mean of PCs have been added)
	fpc[,1:ncca] = unscaled
	fpc[,-1:-ncca] = pcmeans[rep(1, nrow(pcs.rain)),]	# Add in means of remaining PCs
	recon.contemp = fpc %*% t(zsvd.rain$u)				# Scaled reconstructed contemporary rain values (sd and mean may need to be added)
	recon.contemp.unsc = t(t(recon.contemp) * contemp.sd + contemp.mu)

	perc = ((recon.hol.unsc - contemp.mu)/contemp.mu)*100
	
	# write.table(perc, 
	# 	file=paste(paste("./postfiles/ISM/EqPac_both_summer_4npc_", yrBP, "ka_rainperc.txt", sep=""),
	# 	sep="\t"))

	# if (i == 11){
	# 	write.table(recon.contemp, 
	# 		file=paste(paste("./postfiles/ISM/EqPac/CCA-PC", npc, "/", region, "_", proxy, "_", season, "_", npc, "npc", "_", yrBP, "ka_1949-2013.txt", sep=""),
	# 		sep="\t"))

	# 	write.table(data.u, 
	# 		file=paste(paste("./postfiles/EqPac/CCA-PC", npc, "/", region, "_", proxy, "_", season, "_", npc, "npc", "_actual_1949-2013.txt", sep=""),
	# 		sep="\t"))
	# }



	# ####################### CALIB STATISTICS #######################
	# cor.sc = c()
	# cor.unsc = c()

	# for (k in 1:length(index)){
	# 	xcor1 = cor(recon.contemp[,k],sc.precip[,k])^2
	# 	cor.sc = c(cor.sc, xcor1)

	# 	xcor2 = cor(recon.contemp.unsc[,k],precip.davg[,k])^2
	# 	cor.unsc = c(cor.unsc, xcor2)
	# }

	# sumdiff = apply((sc.precip - recon.contemp)^2,2,sum)
	# sumYraw = apply(sc.precip^2,2,sum)
	# beta.sc = 1 - (sumdiff/sumYraw)

	# sumdiff = apply((precip.davg - recon.contemp.unsc)^2,2,sum)
	# sumYraw = apply(precip.davg^2,2,sum)
	# beta.unsc = 1 - (sumdiff/sumYraw)


	### Statistics on just the signal data
	nkeep = ncca
	Ev = matrix(0, nrow = dim(pcs.rain)[2], ncol = dim(pcs.rain)[2])
	Ev[,1:nkeep] = zsvd.rain$u[,1:nkeep]

	rain.sig = pcs.rain %*% t(Ev)
	rain.sig.unsc = t(t(rain.sig) * contemp.sd + contemp.mu)

	# Do stats on rain.sig.unsc vs. recon.contemp.unsc

	R2_rain = c()
	for (k in 1:length(index)){
		R2_rain = c(R2_rain, cor(recon.contemp.unsc[,k], rain.sig.unsc[,k])^2)
	}

	sumdiff = apply((rain.sig.unsc - recon.contemp.unsc)^2,2,sum)
	sumYraw = apply(rain.sig.unsc^2,2,sum)
	beta_rain = 1 - (sumdiff/sumYraw)

	# write(R2_rain, file="./postfiles/ISM/diagnostics/calib_r2_unscsig_rain.txt")
	# write(beta_rain, file="./postfiles/ISM/diagnostics/calib_b_unscsig_rain.txt")

	# zY = (rain.sig * contemp.sd) + contemp.mu
	# zYhat = (recon.contemp * contemp.sd) + contemp.mu

	# sumdiff = apply((zY - zYhat)^2,2,sum)
	# sumYraw = apply(zY^2, 2, sum)
	# beta.signal = 1 - (sumdiff/sumYraw)

	# cor.signal <- c()
	# mse.rain = c()
	# mse.rain2 = c()
	# for (i in 1:length(index)){
	# 	xcor <- cor(zY[,i], zYhat[,i])^2
	# 	cor.signal <- c(cor.signal, xcor)

	# 	mse.rain = c(mse.rain, sum((zY[,i] - zYhat[,i])^2)/(ntime-1))
	# 	# mse.rain2 = c(mse.rain2, sum((sc.precip[,i] - recon.contemp[,i])^2)/(ntime-1))
	# }

	# mean.davg = apply(precip.davg,2,mean)
	# mean.davg.signal = apply(zY,2,mean)
	# mseperc = mse.rain*100/mean.davg
	# # mseperc2 = mse.rain2*100/mean.davg

}


####################################################################################
####################################################################################
### PLOT RAINFALL TIME SERIES
labels = c("W. Ghats", "Indo-Gangetic Plains", "Godavari", "Rajasthan", "N.E. India", "N. India")
means.davg = apply(Y.davg,2,mean)
means.sestot = apply(Y.sestot,2,mean)

prow=6
for (i in 1:prow){
	perc.d.max = ((zpred.davg[,i]+(zpred.davg.se[,i]*1.95))-means.davg[i])*100/means.davg[i]
	perc.d.pts = ((zpred.davg[,i]-means.davg[i])*100)/means.davg[i]
	perc.d.min = ((zpred.davg[,i]-(zpred.davg.se[,i]*1.95))-means.davg[i])*100/means.davg[i]


	dev.new(width=6, height=4.5)
	min = min(floor(perc.d.min))
	max = max(ceiling(perc.d.max))
	# min = -75
	# max = 75

	plot(1:10, perc.d.pts[2:11], type="l", xlim = c(10,1), ylim = c(min,max),
		xlab = "", ylab="",axes=FALSE,col=colors[i], cex.axis=1.25)
	title(main=labels[i], col.main=colors[i],xlab="ka BP",ylab=expression(paste(Delta,"% from present day")), mgp=c(2.5,2.5,1), cex.lab=1.25)
	polygon(c(1:10,10:1),c(perc.d.min[2:11],rev(perc.d.max[2:11])),
		lty=2,col=acolors[i],border=colors[i])
	points(1:10,perc.d.pts[2:11],col=colors[i])
	# lines(0:10, zpred.davg[,i]-(zpred.davg.se[,i]*1.95), col="grey25")
	# lines(0:10, zpred.davg[,i]+(zpred.davg.se[,i]*1.95), col="grey25")
	abline(v = 0:10, col = "grey")
	abline(h = 0, lty=2, col="grey50")
	axis(2, cex.axis=1.25)
	# axis(2, at = c(-100,-60,-20,20),cex.axis=1.25)
	axis(1, at = seq(0,10,1),cex.axis=1.25)
}



####################################################################################
####################################################################################
### PLOT SPATIAL RAINFALL

# If youre doing seasonal total
# mean.sestot = apply(precip.sestot,2,mean)
# perc = ((recon.hol.unsc - mean.sestot)/mean.sestot)*100

# If youre doing daily avg
mean.davg = apply(precip.davg,2,mean)
perc = ((recon.hol.unsc - mean.davg)/mean.davg)*100
perc.signal = ((recon.hol.unsc - mean.davg.signal)/mean.davg.signal)*100

# Plot reconstruction values
dev.new(width=5,height=5)
zfull = rep(NaN,(nx*ny))
zfull[index] = perc
zmat = matrix(zfull,nrow=nx,ncol=ny)
image(xgrid,ygrid,zmat,
	ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)),
	xlab="",ylab="", yaxt="n", xaxt="n",
	cex.axis=1.25, col=rev(myPalette1(100)), zlim=c(-60,60)
	)

axis(1,seq(70,100,10),labels=FALSE)
#axis(1,seq(70,100,10), cex.axis=1.25)
axis(2,c(10,15,20,25,30,35),cex.axis=1.25)
#axis(2,c(10,15,20,25,30,35),cex.axis=1.25, labels=FALSE)
mapnames <- map("world2", 
	xlim=c(min(ygrid), 98), ylim=c(min(ygrid), max(ygrid)), 
	boundary=TRUE, interior=TRUE, col="grey50", add=TRUE
)
contour(xgrid, ygrid, zmat, add = TRUE, levels = seq(-45,45,15), labcex=1.2, lwd = 1.2, col="black")


# Plot the climatology
dev.new(width=5, height=5)
zfull = rep(NaN,(nx*ny))
zfull[index] = mean.davg
zmat = matrix(zfull,nrow=nx,ncol=ny)
image(xgrid,ygrid,zmat,
	ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)),
	xlab="",ylab="", yaxt="n", xaxt="n", cex.axis=1.5, col=rev(myPalette3(100)), zlim=c(0,30)
	)
# image(xgrid,ygrid,zmat,ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)),
# 	xlab="",ylab="", yaxt="n", xaxt="n", cex.axis=1.5, col="firebrick4", zlim=c(100,max(meanclima)), add=TRUE
# 	)
# image(xgrid,ygrid,zmat,ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)),
# 	xlab="",ylab="", yaxt="n", xaxt="n", cex.axis=1.5, col="dodgerblue4", zlim=c(min(meanclima),-100), add=TRUE
# 	)
#axis(1,seq(30,100,10),labels=FALSE)
axis(1,seq(70,100,10), cex.axis=1.25)
axis(2,c(10,15,20,25,30,35),cex.axis=1.25)
mapnames <- map("world2", 
	xlim=c(min(ygrid), 98), ylim=c(min(ygrid), max(ygrid)), 
	boundary=TRUE, interior=TRUE, col="grey50", add=TRUE
)
contour(xgrid, ygrid, zmat, add = TRUE, levels = seq(0,30,5), lwd = 1.2, labcex=1.25, col="black")


####################################################################################
####################################################################################
### MSE Error
#recon.contemp.unsc = (recon.contemp*contemp.sd)+contemp.mu
# mean.davg = apply(precip.davg,2,mean)

# mse = c()
# sdy = c()
# sdyhat = c()
# mseperc = c()
# for (j in 1:length(index)){
# 	y = Rret[,j]
# 	yhat = recon.contemp[,j]
# 	xmse = sum((y - yhat)^2)/(length(y)-1)

# 	mse = c(mse, xmse)
# 	sdy = c(sdy, sd(y))
# 	sdyhat = c(sdyhat, sd(yhat))
# }

# mseperc = mse*100/mean.davg

dev.new(width=7,height=7)
zfull = rep(NaN,(nx*ny))
zfull[index] = mseperc
zmat = matrix(zfull,nrow=nx,ncol=ny)
image(xgrid,ygrid,zmat,
	ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)),
	xlab="",ylab="", yaxt="n", xaxt="n",
	col=rev(myPalette3(100)), zlim=c(0,1.2)
	)
axis(1,seq(30,100,10), cex.axis=2)
axis(2,cex.axis=2)
mapnames <- map("world2", 
	xlim=c(min(ygrid), 98), ylim=c(min(ygrid), max(ygrid)), 
	boundary=TRUE, interior=TRUE, col="grey50", add=TRUE
)
contour(xgrid, ygrid, zmat, add = TRUE, levels = seq(0,1.2,0.2), lwd = 1.2, labcex=1.25, col="black")



####################################################################################
####################################################################################
### Statistics
dev.new(width=7,height=7)
zfull = rep(NaN,(nx*ny))
zfull[index] = beta.unsc
zmat = matrix(zfull,nrow=nx,ncol=ny)
image(xgrid,ygrid,zmat,
	ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)),
	xlab="",ylab="", yaxt="n", xaxt="n", main=expression(paste("Rainfall Calibration ", beta)),
	cex.main=2, col=rev(myPalette4(100)), zlim=c(-1,1)
	)
axis(1,seq(30,100,10), cex.axis=2)
axis(2,cex.axis=2)
mapnames <- map("world2", 
	xlim=c(min(ygrid), 98), ylim=c(min(ygrid), max(ygrid)), 
	boundary=TRUE, interior=TRUE, col="black", add=TRUE
)
contour(xgrid, ygrid, zmat, add = TRUE, levels = seq(-1,1,0.1), lwd = 2, labcex=1.25, col="black")

dev.new(width=7,height=7)
zfull = rep(NaN,(nx*ny))
zfull[index] = cor.unsc
zmat = matrix(zfull,nrow=nx,ncol=ny)
image(xgrid,ygrid,zmat,
	ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)),
	xlab="",ylab="", yaxt="n", xaxt="n", main=expression(paste("Rainfall Calibration ", R^{2})),
	cex.main=2, col=rev(myPalette2(100)), zlim=c(0,1)
	)
axis(1,seq(30,100,10), cex.axis=2)
axis(2,cex.axis=2)
mapnames <- map("world2", 
	xlim=c(min(ygrid), 98), ylim=c(min(ygrid), max(ygrid)), 
	boundary=TRUE, interior=TRUE, col="black", add=TRUE
)
contour(xgrid, ygrid, zmat, add = TRUE, levels = seq(0,0.25,0.1), lwd = 2, labcex=1.25, col="black")



####################################################################################
####################################################################################
### Plot comparison years
iix = which(1901:2004 == 1997)
year = "1997-1998"
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
contour(xgrid, ygrid, zmat, add = TRUE, levels = 0, lwd = 2, labcex=1.25, col="black")


# Reconstructed
zfull = rep(NaN,(nx*ny))
zfull[index] = recon.contemp[iix,]*contemp.sd  						# add in [index1] if youre plotting India rainfall or SST
zmat = matrix(zfull,nrow=nx,ncol=ny)
image(xgrid,ygrid,zmat,
	ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)),
	xlab="",ylab="", main=paste("Reconstructed", year),
	cex.axis=1.1, yaxt="n", col=rev(myPalette1(100)), zlim=c(-11,11)
	)
axis(2, at = c(10,20,30))
mapnames <- map("world2", 
	xlim=c(min(ygrid), max(xgrid)), ylim=c(min(ygrid), max(ygrid)), 
	boundary=TRUE, interior=TRUE, add=TRUE, lwd=1.5, col="grey50"
)
contour(xgrid, ygrid, zmat, add = TRUE, levels = 0, lwd = 2, labcex=1.25, col="black")




