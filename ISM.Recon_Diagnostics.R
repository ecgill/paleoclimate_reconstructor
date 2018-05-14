rm(list=ls())

library(RColorBrewer)
library(fields)

myPalette1 <- colorRampPalette(rev(brewer.pal(9, "RdBu")), space="Lab")
myPalette2 <- colorRampPalette(rev(brewer.pal(9, "Greens")), space="Lab")
myPalette3 <- colorRampPalette(rev(brewer.pal(9, "YlOrRd")), space="Lab")
myPalette4 <- colorRampPalette(rev(brewer.pal(9, "BrBG")), space="Lab")
myPalette10 <- colorRampPalette(rev(brewer.pal(9, "Spectral")), space="Lab")

cbPalette <- c("#D55E00","#E69F00", "#56B4E9", "#009E73", 
	"#F0E442", "#0072B2", "#CC79A7","#999999",
	"#CC0000","#9900CC","#66FF66","#996633","#838B9B")

# Read in diagnostic data.
setwd("/Users/emilygill/Documents/University of Colorado/PHD Research/6. SST Reconstruction/REDO Sep2015/")
eofs.curl = read.table("./postfiles/ISM/diagnostics/EqPac_both_summer_zvsd.curl$u.txt",header=TRUE)
pcs.curl = read.table("./postfiles/ISM/diagnostics/EqPac_both_summer_pcs.curl.txt",header=TRUE)
lambdas.curl = scan("./postfiles/ISM/diagnostics/EqPac_both_summer_lambdas.curl.txt")

eofs.rain = read.table("./postfiles/ISM/diagnostics/EqPac_both_summer_zvsd.raindavg$u.txt",header=TRUE)
pcs.rain = read.table("./postfiles/ISM/diagnostics/EqPac_both_summer_pcs.raindavg.txt",header=TRUE, row.names=NULL)
lambdas.rain = scan("./postfiles/ISM/diagnostics/EqPac_both_summer_lambdas.raindavg.txt")
pcs.rain = pcs.rain[,-1]

eofs.p = read.table("./postfiles/ISM/diagnostics/EqPac_both_summer_zvsd.p$u.txt",header=TRUE)
pcs.p = read.table("./postfiles/ISM/diagnostics/EqPac_both_summer_pcs.p.txt",header=TRUE)
lambdas.p = scan("./postfiles/ISM/diagnostics/EqPac_both_summer_lambdas.p.txt")
recordnopt = read.table("./postfiles/EqPac/PC3/diagnostics/EqPac_both_annual_recordnopt",header=TRUE)

eofs.IOsst = read.table("./postfiles/ISM/diagnostics/EqPac_both_summer_zvsd.IOsst$u.txt",header=TRUE)
pcs.IOsst = read.table("./postfiles/ISM/diagnostics/EqPac_both_summer_pcs.IOsst.txt",header=TRUE)
lambdas.IOsst = scan("./postfiles/ISM/diagnostics/EqPac_both_summer_lambdas.IOsst.txt")

icurl = scan("./postfiles/ISM/diagnostics/index-curl.txt")
irain = scan("./postfiles/ISM/diagnostics/index-rain.txt")
iiosst = scan("./postfiles/ISM/diagnostics/index-iosst.txt")

curl_calib_R2 = scan("./postfiles/ISM/diagnostics/calib_r2_unscsig_curl.txt")
curl_calib_b = scan("./postfiles/ISM/diagnostics/calib_b_unscsig_curl.txt")
rain_calib_R2 = scan("./postfiles/ISM/diagnostics/calib_r2_unscsig_rain.txt")
rain_calib_b = scan("./postfiles/ISM/diagnostics/calib_b_unscsig_rain.txt")

curl_valid_R2 = scan("./postfiles/ISM/diagnostics/valid_r2_unscsig_curl.txt")
curl_valid_b = scan("./postfiles/ISM/diagnostics/valid_b_unscsig_curl.txt")
rain_valid_R2 = scan("./postfiles/ISM/diagnostics/valid_r2_unscsig_rain.txt")
rain_valid_b = scan("./postfiles/ISM/diagnostics/valid_b_unscsig_rain.txt")

# Read in reconstruction results.
years = seq(0,10,1)	
reconSST.curl = c()
reconSST.rain = c()

for (i in 1:length(years)){
	curl = read.table(paste("./postfiles/ISM/EqPac_both_summer_8npc_", i-1, "ka_curl.txt", sep=""),header=TRUE)
	rain = read.table(paste("./postfiles/ISM/EqPac_both_summer_4npc_", i-1, "ka_rainperc.txt", sep=""),header=TRUE)

	reconSST.curl = rbind(reconSST.curl,curl)
	reconSST.rain = rbind(reconSST.rain,rain)

}

xrain = seq(66.5,100.5,by=1)
nxrain = length(xrain) # nrows = 35
yrain = seq(6.5,38.5,by=1)
nyrain = length(yrain) # ncols = 33

xcurl = seq(35.625, 76.875, by=1.875)
nxcurl = length(xcurl) # nrows = 33
ycurl = seq(-6.666573, 31.42808, by=1.904732)
nycurl = length(ycurl) # ncols = 23

xiosst = seq(30,106,by=2)
nxiosst = length(xiosst)
yiosst = seq(-10,28,by=2)
nyiosst = length(yiosst)

xgrid = seq(100,300,by=2)
nx = length(xgrid)
ygrid = seq(-10,10,by=2)
ny = length(ygrid)

### FIGURE PLAN FOR MAIN PAPER
# FIGURE 1 = PROXY MAP (separate code)
# FIGURE 2 = METHODS (powerpoint)
# FIGURE 3 = RECONSTRUCTED CURL AND RAINFALL
# FIGURE 4 = GUPTA COMPARISON (separate code)
# FIGURE 5 = REGIONAL RAINFALL RECONSTRUCTIONS (separate code)

### FIGURE PLAN FOR APPENDIX
# FIG. A1 = INDIAN OCEAN FIRST PC TREND
# FIG. A2 = EIGENVALUE SPECTRA 
# FIG. A3 = CURL & RAIN PCS AND EOFS 
# FIG. A4 = CALIBRATION & VERIFICATION
# OPTIONAL = LIMITED SST PCS AND EOFS
# OPTIONAL = STANDARD ERRORS?


########################################################## FIG 2. SST proxy coverage.
wplat = c(8.8,6.635,6.3,6.514,1.25,0.32,-1.5,-4.689,-3.56,-5.003,-5.9,-7.4,-9.649,-10.592,
	9.233,8.729,5.65,6.158,6.48,1.25,-6.543)
wplon = c(121.3,113.409,125.83,126.498,146.14,159.35,100.1,117.903,119.4,133.445,103.3,115.2,118.338,125.388,
	109.383,109.869,110.65,112.213,125.83,146.14,103.833)
wpno = 1:21
wptype = c(rep("M",14), rep("U",7))
wpproxy = as.data.frame(cbind(wpno,wplat,wplon,wptype))
wpused = c("S","Y","Y","S","S","S","Y","Y","Y","Y","Y","Y","Y","Y",
	"Y","Y","S","Y","N","S","Y")

dev.new(width=14, height=4)
world(ylim=c(-11,11), xlim=c(95,160), lwd=0.25)
abline(h=c(-10,-5,0,5,10), v=seq(90,160,5), col="grey80", lwd=0.5)
world(ylim=c(-11,11), xlim=c(95,160), fill=TRUE, col="grey80", border="black", lwd=0.5, add=TRUE)
map.axes()
points(wplon, wplat,  
	bg = ifelse(wpused == "Y","#0072B2",ifelse(wpused == "N", "#E69F00", "red")),
	pch = ifelse(wptype == "M",22,21), cex = ifelse(wptype == "M",2,1.25)
	)
legend("topright",c("Mg/Ca-based SST","Uk37-based SST"), pch=c(22,21), pt.cex=c(2,1.25), bg="white")


eplat = c(7.85,0.515,0.022,-1.217,-2.51,
	8.206,7.85,4.847,2.25,1.5,0.52,-0.467,-1.217,-1.517,-2.51,-1.85,-3.383,-3.59)
eplon = c(-83.608,-92.398,-86.446,-89.683,-84.650,
	-84.122,-83.6,-77.963,-90.95,-86.485,-92.40,-82.667,-89.683,-85.817,-84.650,-82.78,-83.517,-81.18)
epno = 22:39
eptype = c(rep("M",5), rep("U",13))
epproxy = as.data.frame(cbind(epno,eplat,eplon,eptype))
epused = c("Y","Y","Y","Y","N",
	"Y","N","Y","Y","Y","N","Y","N","S","Y","Y","Y","Y")

dev.new(width=4, height=4)
world(ylim=c(-11,11), xlim=c(-95,-75), lwd=0.25)
abline(h=c(-10,-5,0,5,10), v=seq(-95,-75,5), col="grey80", lwd=0.5)
world(ylim=c(-11,11), xlim=c(-100,-75), fill=TRUE, col="grey80", border="black", lwd=0.5, add=TRUE)
map.axes()
points(eplon, eplat,  
	bg = ifelse(epused == "Y","#0072B2",ifelse(epused == "N", "#E69F00", "red")),
	pch = ifelse(eptype == "M",22,21), cex = ifelse(eptype == "M",2,1.25)
	)

########################################################## FIG 3. Reconstructed curl and rain.

prows = 5
pyrs = seq(2,10,2)
for (i in 1:prows){
	dev.new(width=6, height=5)
	zfull = rep(NaN,(nxcurl*nycurl))
	zfull[icurl] = unname(as.numeric(reconSST.curl[(pyrs[i]+1),]))
	zmat = matrix(zfull,nrow=nxcurl,ncol=nycurl)
	image(xcurl,ycurl,zmat,
		ylim=range(min(ycurl),40), xlim=range(min(xcurl),100),
		xlab="",ylab="", yaxt="n", xaxt="n", cex.axis=1.5, col=myPalette1(100),
		zlim = c(-35,35)
		)
	image(xcurl,ycurl,zmat,
		ylim=range(min(ycurl),40), xlim=range(min(xcurl),100),
		xlab="",ylab="", yaxt="n", xaxt="n", cex.axis=1.5, col="firebrick4",
		zlim = c(35.1,2500), add=TRUE
		)
	image(xcurl,ycurl,zmat,
		ylim=range(min(ycurl),40), xlim=range(min(xcurl),100),
		xlab="",ylab="", yaxt="n", xaxt="n", cex.axis=1.5, col="dodgerblue4",
		zlim = c(-375,-35.1), add=TRUE
		)
	axis(2, cex.axis=1.5)
	world(ylim=range(min(ycurl),40), xlim=range(min(xcurl),100),
		add=TRUE, fill=TRUE, col="grey80", border="grey50"
		)	
	contour(xcurl,ycurl,zmat,add=TRUE,levels = c(-30,-15,0,15,30),lwd=1.25,labcex=0.9)
	points(57.36, 18.05, pch=21, col="black", bg="#009E73", cex=1.5)
	abline(h = 0, col="grey50", lty=2)

	zfull = rep(NaN, (nxrain*nyrain))
	zfull[irain] = unname(as.numeric(reconSST.rain[(pyrs[i]+1),]))
	zmat = matrix(zfull,nrow=nxrain,ncol=nyrain)
	image(xrain,yrain,zmat,
		ylim=range(min(yrain),max(yrain)), xlim=range(min(xrain),max(xrain)),
		xlab="",ylab="", yaxt="n", xaxt="n", cex.axis=1.5, col=rev(myPalette4(100)),
		zlim = c(-62,62), add=TRUE
		)
	contour(xrain,yrain,zmat,add=TRUE,levels = seq(-60,60,15),lwd=1.25,labcex=0.9)
	world(ylim=range(min(ycurl),40), xlim=range(min(xcurl),100),
		add=TRUE, col="grey50"
		)
	if (i == prows){
		axis(1, cex.axis=1.5)
	}
}


########################################################## FIG A1. Indian Ocean First PC Trend

dev.new(width=6,height=5)
zfull = rep(NaN,(nxiosst*nyiosst))
zfull[iiosst] = unname(as.numeric(eofs.IOsst[,1]))
zmat = matrix(zfull,nrow=nxiosst,ncol=nyiosst)
image.plot(xiosst,yiosst,zmat,
	ylim=range(min(yiosst),max(yiosst)), xlim=range(min(xiosst),max(xiosst)),
	xlab="",ylab="", yaxt="n", xaxt="n", cex.axis=1.5, col=rev(myPalette1(100)),
	zlim = c(-0.07,0.07), main = "SST EOF no.1", axes=FALSE, horizontal=TRUE
	)
contour(xiosst,yiosst,zmat,add=TRUE,levels = seq(-1,1,0.01),lwd=1.25,labcex=0.9)
world(ylim=range(min(yiosst),max(yiosst)), xlim=range(min(xiosst),max(xiosst)),
	add=TRUE, col="grey50"
	)
axis(1,at=seq(30,100,20), cex.axis=1.5)
axis(2,at=seq(-10,30,10), cex.axis=1.5)

dev.new(width=6, height=4)
plot(1949:2013, scale(pcs.IOsst[,1]), type="l", main="SST PC no.1", 
	axes=FALSE, ylab="", xlab="")
axis(2, cex.axis=1.5)
axis(1, at = seq(1949,2013,5), cex.axis=1.5)

dev.new(width=4, height=6)
tot=10
plot(100 * lambdas.IOsst[1:tot], type = "l", ylab = "Variance Explained (%)", 
	main = "Indian Ocean SST Eigenvalue Spectrum", 
    bty = "n", xaxt = "n", xlab = "", ylim=c(0,100), lwd=2)
abline(v = 1:tot, h=seq(0,100,20), col = "grey")
points(100 * lambdas.IOsst[1:tot], type="p", pch=20, col="black", cex=1)
axis(1, at = seq(1,tot,1))
title(xlab = "Modes")


########################################################## FIG A2. EVS for all three
tot = 10
dev.new(width=4,height=6)
plot(100 * cumsum(lambdas.curl)[1:tot], type = "l", ylab = "Cumulative Percentage of Variance Explained (%)", 
	main = "Eigenvalue Spectra", 
    bty = "n", xaxt = "n", xlab = "", ylim=c(0,100), lwd=2, col="#B2182B")
abline(v = 1:tot, h=seq(10,100,10), col = "grey80")
lines(100 * cumsum(lambdas.rain)[1:tot], col="#00665E", lwd=2)
lines(100 * cumsum(lambdas.p)[1:tot], lwd=2)

segments(4,0,4,96.3, lty=2, col="#00665E", lwd=1.5)		# RAIN
segments(4,37.3,0,37.3, lty=2, col="#00665E", lwd=1.5)
segments(4,96.3,0,96.3, lty=2, col="#00665E", lwd=1.5)
text(9,50.5,"Precip", srt=20, col="#00665E") 

segments(8,0,8,99.1, lty=2, col="#B2182B", lwd=1.5)		# CURL
segments(8,75.5,0,75.5, lty=2, col="#B2182B", lwd=1.5)
segments(8,99.1,0,99.1, lty=2, col="#B2182B", lwd=1.5)
text(9,76,"Curl", srt=20, col="#B2182B") 

text(9,97,"Lim SST") 
text(9,50.5,"Precip", srt=20, col="#00665E") 
text(9,76,"Curl", srt=20, col="#B2182B") 
text(3.75,12,"4 Modes Retained",srt=90,col="grey40", cex=0.75)
text(7.75,12,"8 Modes Retained",srt=90,col="grey40", cex=0.75)
text(1.3,39,"37.3%", col="#00665E", cex=0.75)
text(1.3,94,"96.3%", col="#00665E", cex=0.75)
text(1.3,77.5,"75.5%", col="#B2182B", cex=0.75)
text(1.3,101,"99.1%", col="#B2182B", cex=0.75)

axis(1, at = seq(1,tot,1))
title(xlab = "Modes")



########################################################## FIG A3. EOFs and PCs for curl and rain
prow = 4
dev.new(height=10,width=3.5)
par(mfrow = c(prow, 1), mar = c(2.5, 2.5, 1, 1))
for (i in 1:prow){
	scaleEOF = eofs.curl[,i]*(lambdas.curl[i]/lambdas.curl[1])
	zfull = rep(NaN,(nxcurl*nycurl))
	zfull[icurl] = scaleEOF
	zmat = matrix(zfull,nrow=nxcurl,ncol=nycurl)
	image(xcurl,ycurl,zmat,
		ylim=range(min(ycurl),40), xlim=range(min(xcurl),100),
		xlab="",ylab="", yaxt="n", xaxt="n", cex.axis=1.5, col=myPalette1(100),
		zlim = c(-0.15,0.15), main = paste("EOF no.", i, sep="")
		)
	axis(2, cex.axis=1.5)
	world(ylim=range(min(ycurl),40), xlim=range(min(xcurl),100),
		add=TRUE, fill=TRUE, col="grey80", border="grey50"
		)	
	contour(xcurl,ycurl,zmat, add=TRUE, levels=c(-0.15,-0.1,-0.05,0,0.05,0.1,0.15), lwd=1.25, labcex=0.9)
	points(57.36, 18.05, pch=21, col="black", bg="#009E73", cex=1.5)
	abline(h = 0, col="grey50", lty=2)	

	scaleEOF = eofs.rain[,i]*(lambdas.rain[i]/lambdas.rain[1])
	zfull = rep(NaN,(nxrain*nyrain))
	zfull[irain] = scaleEOF
	zmat = matrix(zfull,nrow=nxrain,ncol=nyrain)
	image(xrain,yrain,zmat,
		ylim=range(min(yrain),max(yrain)), xlim=range(min(xrain),max(xrain)),
		xlab="",ylab="", yaxt="n", xaxt="n", cex.axis=1.5, col=rev(myPalette4(100)),
		zlim = c(-0.10,0.10), add=TRUE
		)
	contour(xrain,yrain,zmat, add=TRUE, levels=c(-0.1,-0.05,0,0.05,0.1), lwd=1.25, labcex=0.9)
	world(ylim=range(min(ycurl),40), xlim=range(min(xcurl),100),
		add=TRUE, col="grey50"
		)
	if (i == prow){
		axis(1, cex.axis=1.5)
	}
}

# PLOT OVERLYING PCs
prow=4
dev.new(height=10,width=4)
par(mfrow = c(prow, 1), mar = c(2.5, 2.5, 1, 2.5))
full = 1901:2013
for (i in 1:prow){
	curlpc = rep(NA,length(1901:2013))
	curlpc[which(full==1949):which(full==2013)] = pcs.curl[,i]

	rainpc = rep(NA,length(1901:2013))
	rainpc[which(full==1901):which(full==2004)] = pcs.rain[,i]

	plot(full, curlpc, type="l", 
		col="#B2182B", main=paste("PC no.", i, sep=""), axes=FALSE)
	axis(2, cex.axis=1.5, col="#B2182B")
	par(new=TRUE)
	plot(full, rainpc, type="l",
		col="#00665E", axes=FALSE)
	axis(4, cex.axis=1.5,col="#00665E")
	if (i == prow){
		axis(1, at = seq(1901, 2013, 10), cex.axis=1.5)
		legend(1901,10,c("Curl","Rain"), lty=c(1,1), c("#B2182B","#00665E"))
	}
}

########################################################## FIG A4. Calibration and Validation

#### CALIBRATION #####

dev.new(width=6, height=8)
par(mfrow=c(2,1), mar=c(2,2,2,2))

# Curl R2
zfull = rep(NaN,(nxcurl*nycurl))
zfull[icurl] = unname(as.numeric(curl_calib_R2))
zmat = matrix(zfull,nrow=nxcurl,ncol=nycurl)
image(xcurl,ycurl,zmat,
	ylim=range(min(ycurl),40), xlim=range(min(xcurl),100),
	xlab="",ylab="", yaxt="n", xaxt="n", cex.axis=1.5, col=rev(myPalette2(100)),
	zlim = c(0,1)
	)
world(ylim=range(min(ycurl),40), xlim=range(min(xcurl),100),
	add=TRUE, fill=TRUE, col="grey80", border="grey50"
	)
contour(xcurl,ycurl,zmat,add=TRUE,levels = seq(0,1,0.1),lwd=1.25,labcex=0.9)

# Add Rain R2
zfull = rep(NaN, (nxrain*nyrain))
zfull[irain] = unname(as.numeric(rain_calib_R2))
zmat = matrix(zfull,nrow=nxrain,ncol=nyrain)
image(xrain,yrain,zmat,
	ylim=range(min(yrain),max(yrain)), xlim=range(min(xrain),max(xrain)),
	xlab="",ylab="", yaxt="n", xaxt="n", cex.axis=1.5, col=rev(myPalette2(100)),
	zlim = c(0,1), add=TRUE
	)
contour(xrain,yrain,zmat,add=TRUE,levels = seq(0,1,0.1),lwd=1.25,labcex=0.9)
world(ylim=range(min(ycurl),40), xlim=range(min(xcurl),100),
	add=TRUE, col="grey50"
	)

###
# Curl Beta
zfull = rep(NaN,(nxcurl*nycurl))
zfull[icurl] = unname(as.numeric(curl_calib_b))
zmat = matrix(zfull,nrow=nxcurl,ncol=nycurl)
image.plot(xcurl,ycurl,zmat,
	ylim=range(min(ycurl),40), xlim=range(min(xcurl),100),
	xlab="",ylab="", yaxt="n", xaxt="n", cex.axis=1.5, col=rev(myPalette3(100)),
	zlim = c(-1,1), horizontal=TRUE
	)
world(ylim=range(min(ycurl),40), xlim=range(min(xcurl),100),
	add=TRUE, fill=TRUE, col="grey80", border="grey50"
	)
contour(xcurl,ycurl,zmat,add=TRUE,levels = seq(-1,1,0.2),lwd=1.25,labcex=0.9)

# Add Rain Beta
zfull = rep(NaN, (nxrain*nyrain))
zfull[irain] = unname(as.numeric(rain_calib_b))
zmat = matrix(zfull,nrow=nxrain,ncol=nyrain)
image(xrain,yrain,zmat,
	ylim=range(min(yrain),max(yrain)), xlim=range(min(xrain),max(xrain)),
	xlab="",ylab="", yaxt="n", xaxt="n", cex.axis=1.5, col=rev(myPalette3(100)),
	zlim = c(-1,1), add=TRUE
	)
contour(xrain,yrain,zmat,add=TRUE,levels = seq(-1,1,0.1),lwd=1.25,labcex=0.9)
world(ylim=range(min(ycurl),40), xlim=range(min(xcurl),100),
	add=TRUE, col="grey50"
	)



#### VALIDATION #####

dev.new(width=6, height=8)
par(mfrow=c(2,1), mar=c(2,2,2,2))

# Curl R2
zfull = rep(NaN,(nxcurl*nycurl))
zfull[icurl] = unname(as.numeric(curl_valid_R2))
zmat = matrix(zfull,nrow=nxcurl,ncol=nycurl)
image(xcurl,ycurl,zmat,
	ylim=range(min(ycurl),40), xlim=range(min(xcurl),100),
	xlab="",ylab="", yaxt="n", xaxt="n", cex.axis=1.5, col=rev(myPalette2(100)),
	zlim = c(0,1)
	)
world(ylim=range(min(ycurl),40), xlim=range(min(xcurl),100),
	add=TRUE, fill=TRUE, col="grey80", border="grey50"
	)
contour(xcurl,ycurl,zmat,add=TRUE,levels = seq(0,1,0.1),lwd=1.25,labcex=0.9)

# Add Rain R2
zfull = rep(NaN, (nxrain*nyrain))
zfull[irain] = unname(as.numeric(rain_valid_R2))
zmat = matrix(zfull,nrow=nxrain,ncol=nyrain)
image(xrain,yrain,zmat,
	ylim=range(min(yrain),max(yrain)), xlim=range(min(xrain),max(xrain)),
	xlab="",ylab="", yaxt="n", xaxt="n", cex.axis=1.5, col=rev(myPalette2(100)),
	zlim = c(0,1), add=TRUE
	)
contour(xrain,yrain,zmat,add=TRUE,levels = seq(0,1,0.1),lwd=1.25,labcex=0.9)
world(ylim=range(min(ycurl),40), xlim=range(min(xcurl),100),
	add=TRUE, col="grey50"
	)

###
# Curl Beta
zfull = rep(NaN,(nxcurl*nycurl))
zfull[icurl] = unname(as.numeric(curl_valid_b))
zmat = matrix(zfull,nrow=nxcurl,ncol=nycurl)
image(xcurl,ycurl,zmat,
	ylim=range(min(ycurl),40), xlim=range(min(xcurl),100),
	xlab="",ylab="", yaxt="n", xaxt="n", cex.axis=1.5, col=rev(myPalette3(100)),
	zlim = c(0,1)
	)
world(ylim=range(min(ycurl),40), xlim=range(min(xcurl),100),
	add=TRUE, fill=TRUE, col="grey80", border="grey50"
	)
contour(xcurl,ycurl,zmat,add=TRUE,levels = seq(0,1,0.1),lwd=1.25,labcex=0.9)

# Add Rain Beta
zfull = rep(NaN, (nxrain*nyrain))
zfull[irain] = unname(as.numeric(rain_valid_b))
zmat = matrix(zfull,nrow=nxrain,ncol=nyrain)
image(xrain,yrain,zmat,
	ylim=range(min(yrain),max(yrain)), xlim=range(min(xrain),max(xrain)),
	xlab="",ylab="", yaxt="n", xaxt="n", cex.axis=1.5, col=rev(myPalette3(100)),
	zlim = c(0,1), add=TRUE
	)
contour(xrain,yrain,zmat,add=TRUE,levels = seq(0,1,0.1),lwd=1.25,labcex=0.9)
world(ylim=range(min(ycurl),40), xlim=range(min(xcurl),100),
	add=TRUE, col="grey50"
	)
























########################################################## OPTIONAL. Standard Errors



########################################################## OPTIONAL. EOFs and PCs for limited SST
# prow = 4

# plotdata = as.data.frame(unname(cbind(recordnopt,eofs.p[,1:4])))
# names(plotdata) = c("No.","Lat","Lon","EOF1","EOF2","EOF3","EOF4")
# library(seqinr)
# dblue <- col2alpha("#2066AC",alpha=0.5)
# lblue <- col2alpha("#BCDAEA",alpha=0.65)
# lred <- col2alpha("#FBC8AF",alpha=0.65)
# dred <- col2alpha("#B2182B",alpha=0.5)

# dev.new(height=7,width=8)
# par(mfrow = c(prow, 1), mar = c(2.5, 2.5, 1, 1))
# for (i in 1:prow){
# 	zfull = rep(NaN,(nx*ny))
# 	zmat = matrix(zfull,nrow=nx,ncol=ny)
# 	image(xgrid,ygrid,zmat,
# 		ylim=range(-15,15), xlim=range(95,280), xlab="",ylab="", axes=FALSE,
# 		cex = 1.75, main = paste("EOF no.", i, sep=""), col=myPalette1(100),
# 		zlim=c(-0.12,0.12)
# 		)
# 	mapnames <- map("world2", 
# 	    xlim=c(min(ygrid), max(xgrid)), ylim=c(min(ygrid), max(ygrid)), 
# 	    boundary=TRUE, interior=TRUE,
# 	    col="grey50", add=TRUE
# 	)
# 	axis(2,seq(-10,10,by=5),cex.axis=1.5)

# 	# PLOT THE LIMITED FIELD EOF BUBBLES
# 	EOFlim = plotdata[,i+3]
# 	ratio = lambdas.p[1:4]/lambdas.p[1]
# 	EOFrat = abs(EOFlim)*ratio[i]
# 	radii = sqrt(EOFrat/pi)
# 	xx = rep(-1,27)
# 	xx[which(EOFlim > 0)] = 1
# 	xradii = radii * xx
# 	print(xradii)
# 	symbols(x=plotdata$Lon, y=plotdata$Lat, 
# 		circles=abs(xradii), inches = max(abs(xradii)), 
# 		bg = ifelse(xradii < -.2 ,dblue, 
# 		ifelse(xradii >= -.2 & xradii <= 0,lblue,
# 		ifelse(xradii > 0 & xradii <= 0.2, lred,dred))),
# 		add=TRUE
# 	)

# 	if (i == prow) {
#          axis(1,seq(100,280,by=20),
#          	labels=c("100","120","140","160","180","-160","-140","-120","-100","-80"),
#          	cex.axis=1.5)
#     }
# }

########################################################## OPTIONAL. Standard Errors




