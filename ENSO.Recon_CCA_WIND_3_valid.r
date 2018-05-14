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




###################################################################################################
##### *** PICK AN AREA TO DO THE RECONSTRUCTION OVER
region = "EqPac"						# EqPac or IO
proxy = "both"							# MgCa or Uk
season = "annual"						# winter or summer

train1 = 1980
train2 = 2013
valid1 = 1949
valid2 = 1979



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


# # Remove the 1st PC from IO SSTs because its a trend.
# if (region == "IO"){
# 	Nmodes <- length(zsvd$u[1, ])
# 	Nkeep <- c(1)
# 	E <- matrix(0, nrow = Nmodes, ncol = Nmodes)
# 	E[, Nkeep] <- zsvd$u[, Nkeep]
# 	SSTannavgkeep <- pcs %*% t(E)
# 	data.detrend <- data - SSTannavgkeep

# 	zs <- var(data.detrend)
# 	zsvd <- svd(zs)
# 	pcs <- t(t(zsvd$u) %*% t(data.detrend))
# 	lambdas <- (zsvd$d/sum(zsvd$d))
# }



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
	
} else if (region == "IO"){
	reconpt <- rbind(sst.raw$IOcoord)
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
	smoodat = c(sst.smooth$EP, sst.smooth$WP)
	if (proxy == "both"){
		smoodat = smoodat[-mixed3]
	}
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
yrBP = 10

###################################################################################################
##### 7. ## Pick the year you want to reconstruct for Holocene
missing = which(is.na(mat[yrBP+1,]))
if (length(missing) > 0){
	var = seq(1,nreconpt,1)[-missing]
} else {
	var = seq(1,nreconpt,1)
}


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
index1 <- 1:(nx*ny)
data1 <- data.u[,,1]

# Put the data into a matrix of ntime x mlocation
BINdata.u = matrix(NA, nrow=ntime, ncol=nsta1)
for (i in 1:ntime){
    udata = data.u[,,i]
    data2u = udata[index1]
    BINdata.u[i,] = data2u
}

# Read in climatology and scale data.
CLIMdata.u = clim.u[index1]
data.u = sweep(BINdata.u,2,CLIMdata.u)



###################################################################################################
##### 8. ## PCA ON FULL CONTEMPORARY SST FIELD (AND OPTIONAL DIAGNOSTICS)
### Data for just the training period
aa = which(1949:2013 == train1)
bb = which(1949:2013 == train2)
xdata = data.u[aa:bb,]

## U-WIND
zs.u <- var(xdata)
zsvd.u <- svd(zs.u)
pcs.u <- t(t(zsvd.u$u) %*% t(xdata))
lambdas.u <- (zsvd.u$d/sum(zsvd.u$d))


###################################################################################################
##### 9. ## PCA ON PARTIAL CONTEMPORARY SST FIELD (AND OPTIONAL DIAGNOSTICS)
xreconsst = reconsst[aa:bb,var]
zs.p <- var(xreconsst)
zsvd.p <- svd(zs.p)
pcs.p <- t(t(zsvd.p$u) %*% t(xreconsst))
lambdas.p <- (zsvd.p$d/sum(zsvd.p$d))


###################################################################################################
##### 10. ## VALIDATE THE MODEL FOR TEH VALIDATION YEARS.
#### Fit each retained PC of Y (full SST) with the retained PCs of X (holocene SSTs)
#xx.hol = matrix(mat[yrBP+1,var],1,length(mat[yrBP+1,var]))
#newdata = as.data.frame(xx.hol %*% zsvd.p$u)

cc = which(1949:2013 == valid1)
dd = which(1949:2013 == valid2)
valid = reconsst[cc:dd,var]

npc = 6

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
xx.hol = valid
newdata = xx.hol %*% zsvd.p$u[,1:npc]

### Predict past wind PC field.
pcmeans = matrix(colMeans(pcs.u[,-1:-npc]),nrow=1)
fpc.val = matrix(NA,length(valid1:valid2),ncol(pcs.u))
scaled = newdata %*% AA %*% betahat
unscaled = (scaled * Ysd) + Ymean
fpc.val[,1:npc] = unscaled
fpc.val[,-1:-npc] = pcmeans
recon.val.u = fpc.val %*% t(zsvd.u$u)

### Predict contemporary wind PC field.
fpc = matrix(NA,nrow(pcs.u),ncol(pcs.u))
scaled = X1 %*% AA %*% betahat
unscaled = t(t(scaled) * Ysd + Ymean)
fpc[,1:npc] = unscaled
fpc[,-1:-npc] = pcmeans[rep(1, nrow(pcs.u)),]
recon.contemp.u = fpc %*% t(zsvd.u$u)
						
actualu = data.u[cc:dd,]

# Calibration measures (beta and R2)
sumdiff = apply((actualu - recon.val.u)^2,2,sum)
sumYraw = apply(actualu^2,2,sum)
beta = 1 - (sumdiff/sumYraw)

# R2 correlation
cor <- c()
cor2 = c()
for (i in 1:nsta1){
	xcor <- cor(recon.val.u[,i],actualu[,i])
	cor = c(cor,xcor)
	cor2 = c(cor2,xcor^2)
}

write(beta, 
	file=paste(paste("./postfiles/", region, "/CCA-PC6/diagnostics/", region, "_", proxy, "_", season, "_beta_v_wind.txt", sep=""),
	sep="\t"))
write(cor, 
	file=paste(paste("./postfiles/", region, "/CCA-PC6/diagnostics/", region, "_", proxy, "_", season, "_R2_v_wind.txt", sep=""),
	sep="\t"))

dev.new(width=10, height=4)
par(mfrow = c(2, 1), mar = c(2,2,2,2))
# BETA
zfull = rep(NaN,(nx*ny))
zfull[index1] = beta
zmat = matrix(zfull,nrow=nx,ncol=ny)
image(xgrid,ygrid,zmat,
	ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)),
	xlab="",ylab="", main=expression(paste("Zonal Wind Verification ", beta)),
	cex.axis=1.1, yaxt="n", col=rev(myPalette3(100)),
	zlim=c(-1,1)
	)
axis(2, at = c(-10,0,10))
mapnames <- map("world2", 
	xlim=c(min(ygrid), max(xgrid)), ylim=c(min(ygrid), max(ygrid)), 
	boundary=TRUE, interior=TRUE, add=TRUE, lwd=1.5, col="grey50"
)
contour(xgrid, ygrid, zmat, add = TRUE, levels = seq(-1,1,0.2), lwd = 2, labcex=1.25)

# R^2
zfull = rep(NaN,(nx*ny))
zfull[index1] = cor2
zmat = matrix(zfull,nrow=nx,ncol=ny)
image(xgrid,ygrid,zmat,
	ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)),
	xlab="",ylab="", main=expression(paste("Zonal Wind Verification ", R^{2})),
	cex.axis=1.1, yaxt="n", col=rev(myPalette2(100)),
	zlim=c(0,1)
	)
axis(2, at = c(-10,0,10))
mapnames <- map("world2", 
	xlim=c(min(ygrid), max(xgrid)), ylim=c(min(ygrid), max(ygrid)), 
	boundary=TRUE, interior=TRUE, add=TRUE, lwd=1.5, col="grey50"
)
contour(xgrid, ygrid, zmat, add = TRUE, levels = seq(0,1,0.1), lwd = 2, labcex=1.25)


# ### PLOT THE VALIDATION R
# dev.new(width=10, height=3)
# par(mfrow = c(1, 1), mar = c(4, 4, 4, 1))
# zfull = rep(NaN,(nx*ny))
# zfull[index1] = cor 						# add in [index1] if youre plotting India rainfall
# zmat = matrix(zfull,nrow=nx,ncol=ny)
# image.plot(xgrid,ygrid,zmat,
# 	ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)),
# 	xlab="",ylab="", main = "Validation R", 
# 	cex.axis=1.1, col=rev(myPalette1(100)),
# 	zlim=c(-1,1)
# 	)
# mapnames <- map("world2", 
# 	xlim=c(min(ygrid), max(xgrid)), ylim=c(min(ygrid), max(ygrid)), 
# 	boundary=TRUE, interior=TRUE, col="grey50", add=TRUE
# )
# contour(xgrid, ygrid, zmat, add = TRUE, levels = seq(-1,1,0.1), lwd = 2, labcex=1.25)









