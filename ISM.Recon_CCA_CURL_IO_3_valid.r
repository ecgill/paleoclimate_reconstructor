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

train1 = 1980
train2 = 2013
valid1 = 1949
valid2 = 1979

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

aa = which(1949:2013 == train1)
bb = which(1949:2013 == train2)
xdata.train = scdata.detrend.curl[aa:bb,]

zs.train <- var(xdata.train)
zsvd.train <- svd(zs.train)
pcs.train <- t(t(zsvd.train$u) %*% t(xdata.train))
lambdas.train <- (zsvd.train$d/sum(zsvd.train$d))

###################################################################################################
##### 7. ## Pick the year you want to reconstruct for Holocene
yrBP = 10
missing = which(is.na(mat[yrBP+1,]))
if (length(missing) > 0){
	var = seq(1,nreconpt,1)[-missing]
} else {
	var = seq(1,nreconpt,1)
}


###################################################################################################
##### 9. ## PCA ON PARTIAL CONTEMPORARY SST FIELD (AND OPTIONAL DIAGNOSTICS)
xreconsst = reconsst[aa:bb,var]
zs.p <- var(xreconsst)
zsvd.p <- svd(zs.p)
pcs.p <- t(t(zsvd.p$u) %*% t(xreconsst))
lambdas.p <- (zsvd.p$d/sum(zsvd.p$d))

cc = which(1949:2013 == valid1)
dd = which(1949:2013 == valid2)
valid = reconsst[cc:dd,var]

valid.curl = scdata.detrend.curl[cc:dd,]
zs.valid <- var(valid.curl)
zsvd.valid <- svd(zs.valid)
pcs.valid <- t(t(zsvd.valid$u) %*% t(valid.curl))
lambdas.valid <- (zsvd.valid$d/sum(zsvd.valid$d))


###################################################################################################
##### 7. ## PERFORM CCA to reconstruct U-WIND
npc = 8

Y = pcs.train[,1:npc]
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
pcmeans = matrix(colMeans(pcs.train[,-1:-npc]),nrow=1)
fpc.val = matrix(NA,length(valid1:valid2),ncol(pcs.train))
scaled = newdata %*% AA %*% betahat
unscaled = t(t(scaled) * Ysd + Ymean)
fpc.val[,1:npc] = unscaled
fpc.val[,-1:-npc] = pcmeans
recon.val.curl = fpc.val %*% t(zsvd.train$u)
recon.val.curl.unsc = t(t(recon.val.curl) * contemp.sd.curl + contemp.mu.curl)

# ####################### CALIB STATISTICS (ALL DATA) #######################
# actual.sc = scdata.detrend.curl[cc:dd,]
# actual.unsc = t(t(actual.sc) * contemp.sd.curl + contemp.mu.curl)

# cor.ver.sc = c()
# cor.ver.unsc = c()
# for (i in 1:nsta1){
# 	xcor1 = cor(recon.val.curl[,i], actual.sc[,i])^2
# 	cor.ver.sc = c(cor.ver.sc, xcor1)

# 	xcor2 = cor(recon.val.curl.unsc[,i], actual.unsc[,i])^2
# 	cor.ver.unsc = c(cor.ver.unsc, xcor2)
# }

# sumdiff = apply((actual.sc - recon.val.curl)^2,2,sum)
# sumYraw = apply(actual.sc^2,2,sum)
# beta.ver.sc = 1 - (sumdiff/sumYraw)

# sumdiff = apply((actual.unsc - recon.val.curl.unsc)^2,2,sum)
# sumYraw = apply(actual.unsc^2,2,sum)
# beta.ver.unsc = 1 - (sumdiff/sumYraw)


####################### CALIB STATISTICS (SIGNAL ONLY) #######################
nkeep = npc
Ev = matrix(0, nrow = dim(pcs.train)[2], ncol = dim(pcs.train)[2])
Ev[,1:nkeep] = zsvd.valid$u[,1:nkeep] 

curl.sig.sc = pcs.valid %*% t(Ev)
curl.sig.unsc = t(t(curl.sig.sc) * contemp.sd.curl + contemp.mu.curl)

# The two quantities to compare:
xcurl.unsc = sweep(curl.sig.unsc,2,CLIMdata.curl,'+')
reconcurl.unsc = sweep(recon.val.curl.unsc,2,CLIMdata.curl,'+')

R2_curl <- c()
for (i in 1:nsta1){
	R2_curl = c(R2_curl, cor(reconcurl.unsc[,i], xcurl.unsc[,i])^2)
}

sumdiff = apply((xcurl.unsc - reconcurl.unsc)^2, 2, sum)
sumYraw = apply(xcurl.unsc^2, 2, sum)
beta_curl = 1 - (sumdiff/sumYraw)

write(R2_curl, file="./postfiles/ISM/diagnostics/valid_r2_unscsig_curl.txt")
write(beta_curl, file="./postfiles/ISM/diagnostics/valid_b_unscsig_curl.txt")

####################################################################################
####################################################################################
### Statistics
dev.new(width=4, height=6)
par(mfrow = c(2, 1), mar = c(2,2,2,2))
# BETA
zfull = rep(NaN,(nx*ny))
zfull[index3] = beta_curl
zmat = matrix(zfull,nrow=nx,ncol=ny)
image(xgrid,ygrid,zmat,
	ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)),
	xlab="",ylab="", main=expression(paste("Wind Stress Curl Verification ", beta)),
	cex.axis=1.1, yaxt="n", col=rev(myPalette3(100)),
	zlim=c(-.25,1)
	)
axis(2, at = c(-5, 5, 15, 25))
world(ylim=range(min(ygrid),max(ygrid)),
	xlim=range(min(xgrid),max(xgrid)),
	add=TRUE,
	fill=TRUE,
	border="darkgrey"
)

# R^2
zfull = rep(NaN,(nx*ny))
zfull[index3] = R2_curl
zmat = matrix(zfull,nrow=nx,ncol=ny)
image(xgrid,ygrid,zmat,
	ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)),
	xlab="",ylab="", main=expression(paste("Wind Stress Curl Verification ", R^{2})),
	cex.axis=1.1, yaxt="n", col=rev(myPalette2(100)),
	zlim=c(0,1)
	)
axis(2, at = c(-5, 5, 15, 25))
world(ylim=range(min(ygrid),max(ygrid)),
	xlim=range(min(xgrid),max(xgrid)),
	add=TRUE,
	fill=TRUE,
	border="darkgrey"
)




