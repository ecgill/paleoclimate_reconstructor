rm(list=ls())

library(RColorBrewer)

myPalette1 <- colorRampPalette(rev(brewer.pal(9, "RdBu")), space="Lab")
myPalette2 <- colorRampPalette(rev(brewer.pal(9, "Greens")), space="Lab")
myPalette3 <- colorRampPalette(rev(brewer.pal(9, "YlOrRd")), space="Lab")
myPalette4 <- colorRampPalette(rev(brewer.pal(9, "BrBG")), space="Lab")
myPalette10 <- colorRampPalette(rev(brewer.pal(9, "Spectral")), space="Lab")

cbPalette <- c("#D55E00","#E69F00", "#56B4E9", "#009E73", 
	"#F0E442", "#0072B2", "#CC79A7","#999999",
	"#CC0000","#9900CC","#66FF66","#996633","#838B9B")

### READ IN RESULTS
region = "EqPac"
proxy = "both"
season = "annual"

# Read in diagnostic data
setwd("/Users/emilygill/Documents/University of Colorado/PHD Research/6. SST Reconstruction/REDO Sep2015/")
eofs.full = read.table(paste("./postfiles/", region, "/PC3/diagnostics/", region, "_", proxy, "_", season, "_zvsd$u", sep=""),header=TRUE)
pcs.full = read.table(paste("./postfiles/", region, "/PC3/diagnostics/", region, "_", proxy, "_", season, "_pcs", sep=""),header=TRUE)
lam.full = scan(paste("./postfiles/", region, "/PC3/diagnostics/", region, "_", proxy, "_", season, "_lambdas", sep=""))
eofs.lim = read.table(paste("./postfiles/", region, "/PC3/diagnostics/", region, "_", proxy, "_", season, "_zvsd.p$u", sep=""),header=TRUE)
pcs.lim = read.table(paste("./postfiles/", region, "/PC3/diagnostics/", region, "_", proxy, "_", season, "_pcs.p", sep=""),header=TRUE)
lam.lim = scan(paste("./postfiles/", region, "/PC3/diagnostics/", region, "_", proxy, "_", season, "_lambdas.p", sep=""))
recordnopt = read.table(paste("./postfiles/", region, "/PC3/diagnostics/", region, "_", proxy, "_", season, "_recordnopt", sep=""),header=TRUE)

weofs.full = read.table(paste("./postfiles/", region, "/PC6/diagnostics/", region, "_", proxy, "_", season, "_zvsd.u$u", sep=""),header=TRUE)
wpcs.full = read.table(paste("./postfiles/", region, "/PC6/diagnostics/", region, "_", proxy, "_", season, "_pcs.u", sep=""),header=TRUE)
wlam.full = scan(paste("./postfiles/", region, "/PC6/diagnostics/", region, "_", proxy, "_", season, "_lambdas.u", sep=""))

matanom = read.table(paste("./postfiles/", region, "/PC3/diagnostics/", region, "_", proxy, "_", season, "_mat.txt", sep=""),header=TRUE)
zzpt = scan(paste("./postfiles/", region, "/PC3/diagnostics/", region, "_", proxy, "_", season, "_zzpt.txt", sep=""))

# Read in reconstructed/prediction data
i1 = scan(paste("./postfiles/", region, "/PC3/diagnostics/", region, "_", proxy, "_", season, "_index1", sep=""))
actualsst = read.table(paste("./postfiles/", region, "/PC3/", region, "_", proxy, "_", season, "_2npcr_actual_1901-2013.txt", sep=""),header=TRUE)
reconsst = read.table(paste("./postfiles/", region, "/PC3/", region, "_", proxy, "_", season, "_2npcr_10ka_1901-2013.txt", sep=""),header=TRUE)

actualwind = read.table(paste("./postfiles/", region, "/CCA-PC6/", region, "_", proxy, "_", season, "_6npc_actual_1949-2013.txt", sep=""),header=TRUE)
reconwind = read.table(paste("./postfiles/", region, "/CCA-PC6/", region, "_", proxy, "_", season, "_6npc_10ka_1949-2013.txt", sep=""),header=TRUE)

# Put reconstructed SSTs and uncertainties at each point into a matrix for the comparison/validation plots
years = seq(0,10,1)	
reconSST.mat = c()
reconSST.mat.sd = c()
	
for (i in 1:length(years)){
		test = read.table(paste("./postfiles/", region, "/PC3/", region, "_", proxy, "_", season, "_2npcr", "_", i-1, "ka.txt", sep=""),header=TRUE)
		reconSST.mat = rbind(reconSST.mat,test)

		sd = read.table(paste("./postfiles/", region, "/PC3/", region, "_", proxy, "_", season, "_2npcr", "_", i-1, "ka_sd.txt", sep=""),header=TRUE)
		reconSST.mat.sd = rbind(reconSST.mat.sd,t(sd))
}

reconSST.pt = reconSST.mat[,zzpt]
reconSST.pt.sd = reconSST.mat.sd[,zzpt]
labels = c(27,29,30,31,33,36,37,38,39,22,23,24,25,15,16,18,21,2,3,7,8,9,10,11,12,13,14) #EP Uk, EP MgCa, WP Uk, WP MgCa

xgrid = seq(100,300,by=2)
nx = length(xgrid)
ygrid = seq(-10,10,by=2)
ny = length(ygrid)

wxgrid <- seq(100,300,by=2.5)
nwx = length(wxgrid)
wygrid <- seq(-10,10,by=2.5)
nwy = length(wygrid)

###################################################################
###################################################################
############# MAIN FIGURES ########################################
###################################################################
###################################################################

## ---- FIGURE 4. EX La Nina and El Nino Years, both wind & sst ---- ######
###########################################################################
yr = 1997
iix = which(1854:2013 == yr)
iix2 = which(1949:2013 == yr)
year = "1997-1998"

dev.new(width=10, height=4)
par(mfrow = c(2, 1), mar = c(2,2,2,2))
# Actual
zfull = rep(NaN,(nx*ny))
zfull[i1] = as.numeric(actualsst[iix,])					# add in [index1] if youre plotting India rainfall or SST
zmat = matrix(zfull,nrow=nx,ncol=ny)
image(xgrid,ygrid,zmat,
	ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)),
	xlab="",ylab="", main="Actual",
	cex.axis=1.1, yaxt="n", col=myPalette1(100),
	zlim=c(-5,5)
	)
axis(2, at = c(-10,0,10))
mapnames <- map("world2", 
	xlim=c(min(ygrid), max(xgrid)), ylim=c(min(ygrid), max(ygrid)), 
	boundary=TRUE, interior=TRUE, add=TRUE, lwd=1.5, col="grey50"
)
wfull = as.numeric(actualwind[iix2,])
wmat = matrix(wfull,nrow=nwx,ncol=nwy)
contour(xgrid, ygrid, zmat, add=TRUE, levels = seq(-5,5,1), lwd=1, labcex=1.25, col="black", vfont=c("sans serif","bold"))
contour(wxgrid, wygrid, wmat, add=TRUE, levels = seq(0,6.5,1), lwd = 2, labcex=1.5, col="#2066AC", vfont=c("sans serif","bold"))
contour(wxgrid, wygrid, wmat, add=TRUE, levels = seq(-6.5,0,1), lwd = 2, labcex=1.5, col="#B2182B", vfont=c("sans serif","bold"))
contour(wxgrid, wygrid, wmat, add=TRUE, levels = 0, lwd = 2, labcex=1.5, col="black", vfont=c("sans serif","bold"))

#Reconstructed
zfull = rep(NaN,(nx*ny))
zfull[i1] = as.numeric(reconsst[iix,]) 						# add in [index1] if youre plotting India rainfall or SST
zmat = matrix(zfull,nrow=nx,ncol=ny)
image(xgrid,ygrid,zmat,
	ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)),
	xlab="",ylab="", main="Reconstructed",
	cex.axis=1.1, yaxt="n", col=myPalette1(100),
	zlim=c(-5,5)
	)
axis(2, at = c(-10,0,10))
mapnames <- map("world2", 
	xlim=c(min(ygrid), max(xgrid)), ylim=c(min(ygrid), max(ygrid)), 
	boundary=TRUE, interior=TRUE, add=TRUE, lwd=1.5, col="grey50"
)
#contour(xgrid, ygrid, zmat, add = TRUE, levels = seq(min,max,0.5), lwd = 2, labcex=1.25)
wfull = as.numeric(reconwind[iix2,])
wmat = matrix(wfull,nrow=nwx,ncol=nwy)
contour(xgrid, ygrid, zmat, add=TRUE, levels = seq(-5,5,1), lwd=1, labcex=1.25, col="black", vfont=c("sans serif","bold"))
contour(wxgrid, wygrid, wmat, add=TRUE, levels = seq(0,6.5,1), lwd = 2, labcex=1.5, col="#2066AC", vfont=c("sans serif","bold"))
contour(wxgrid, wygrid, wmat, add=TRUE, levels = seq(-6.5,0,1), lwd = 2, labcex=1.5, col="#B2182B", vfont=c("sans serif","bold"))
contour(wxgrid, wygrid, wmat, add=TRUE, levels = 0, lwd = 2, labcex=1.5, col="black", vfont=c("sans serif","bold"))

## ---- FIGURE 5. Multiproxy Reconstruction (winds & sst) ---- ############
###########################################################################
years = seq(2,10,2)
prow = 5
dev.new(height=5,width=10)
par(mfrow = c(prow, 1), mar = c(0.5, 2, 0.5, 0.5), oma=c(2,1,1,1))
for (i in 1:length(years)){
	test = read.table(paste("./postfiles/", region, "/PC3/", region, "_", proxy, "_", season, "_2npcr_", years[i], "ka.txt", sep=""),header=TRUE)
	#wtest = read.table(paste("./postfiles/", region, "/CCA-PC6/", region, "_", proxy, "_", season, "_6npc_", years[i], "ka_u.txt", sep=""),header=TRUE)
	zfull = rep(NaN,(nx*ny))
	zfull[i1] = as.numeric(unlist(test))
	zmat = matrix(zfull,nrow=nx,ncol=ny)
	image(xgrid,ygrid,zmat,
		ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)), axes=FALSE, 
		col=myPalette1(100), zlim=c(-1,1)
		)
	mapnames <- map("world2", 
		xlim=c(min(ygrid), max(xgrid)), ylim=c(min(ygrid), max(ygrid)), 
		boundary=TRUE, interior=TRUE, col="black", add=TRUE
	)
	# wfull = as.numeric(wtest)
	# wmat = matrix(wfull,nrow=nwx,ncol=nwy)
	contour(xgrid, ygrid, zmat, add=TRUE, levels = seq(-2,2,0.1), lwd=2, labcex=1.25, col="black", vfont=c("sans serif","bold"))
	# contour(wxgrid, wygrid, wmat, add=TRUE, levels = seq(-12,12,2), lwd = 2, labcex=1.25, col="#00665E")
	axis(2,seq(-8,8,by=4),cex.axis=2)	

	if (i == length(years)) {
         axis(1,c(100,140,180,220,260,300),labels=c("100","140","180","-140","-100","-60"),cex.axis=2)
    }

}


## Add in a ggplot version:


years = seq(2,10,2)
prow = 5
dev.new(height=5,width=10)
par(mfrow = c(prow, 1), mar = c(0.5, 2, 0.5, 0.5), oma=c(2,1,1,1))
for (i in 1:length(years)){
	wtest = read.table(paste("./postfiles/", region, "/CCA-PC6/", region, "_", proxy, "_", season, "_6npc_", years[i], "ka_u.txt", sep=""),header=TRUE)
	zfull = rep(NaN,(nwx*nwy))
	zfull = as.numeric(unlist(wtest))
	zmat = matrix(zfull,nrow=nwx,ncol=nwy)
	image(wxgrid,wygrid,zmat,
		ylim=range(min(wygrid),max(wygrid)), xlim=range(min(wxgrid),max(wxgrid)), axes=FALSE, 
		col=myPalette4(100), zlim=c(-8,8)
		)
	mapnames <- map("world2", 
		xlim=c(min(wygrid), max(wxgrid)), ylim=c(min(wygrid), max(wygrid)), 
		boundary=TRUE, interior=TRUE, col="black", add=TRUE
	)
	contour(wxgrid, wygrid, zmat, add=TRUE, levels = seq(-8,8,1), lwd=2, labcex=1.25, col="black", vfont=c("sans serif","bold"))
	axis(2,seq(-8,8,by=4),cex.axis=2)	

	if (i == length(years)) {
         axis(1,c(100,140,180,220,260,300),labels=c("100","140","180","-140","-100","-60"),cex.axis=2)
    }
}


## ---- FIGURE 7. E.Pacific Validation ---- ############
###########################################################################
dev.new(width=8,height=8)
par(mfrow = c(4,4), mar = c(1, 1, 1, 1), oma = c(2, 2, 1, 1))
nEP = 13

for (i in 1:nEP){
	actual = matanom[,i]
	reconsst = reconSST.pt[,i]
	sd = reconSST.pt.sd[,i]
	val = which(is.na(actual) == FALSE)
	miny = -3.5
	maxy = 3.5

	xpoly = c(c(0:10)[val],rev(c(0:10)[val]))
	ypoly = c(actual[val]-1,rev(actual[val]+1))
	plot(0:10, actual, type="l", xlim = rev(range(0:10)), ylim = c(miny,maxy), axes=FALSE, ylab="", xlab="",col="#0072B2")
	polygon(xpoly,ypoly,col=rgb(0.337, 0.706, 0.914 ,0.25),border = "#56B4E9")

	abline(v=c(0:10),col="grey75")
	abline(h=0,col="grey75")

	points(0:10,actual,pch=20, col="#0072B2")
	points(0:10,reconsst,pch=23,bg="white",col="#D55E00")
	arrows(0:10,reconsst+sd,0:10,reconsst-sd,angle=90,code=3,length=0.05, lwd=0.5,col="#E69F00")

	text(8,3.5,paste("proxy no.",labels[i]))
	text(3,3.5,paste("RSS =",round(sum((actual[val] - reconsst[val])^2),2)))
	text(8,-3, paste("EOF1 =",round(eofs.lim[i,1],2)))
	text(3,-3, paste("EOF2 =",round(eofs.lim[i,2],2)))

	if (is.element(i,c(1,5,9,13)) == TRUE){
		axis(2, at = c(-3,-1.5,0,1.5,3))
	} else {
		axis(2, at = c(-3,-1.5,0,1.5,3), labels = c("","","","",""))
	}

	if (is.element(i,c(10,11,12,13)) == TRUE){
		axis(1, at = c(0:10))
	} else {
		axis(1, at = c(0:10), labels = c("","","","","","","","","","",""))
	}
}

## ---- FIGURE 6. W.Pacific Validation ---- ############
###########################################################################
dev.new(width=8,height=10)
par(mfrow = c(5,4), mar = c(1, 1, 1, 1), oma = c(2, 2, 1, 1))		
for (i in 14:27){
	actual = matanom[,i]
	reconsst = reconSST.pt[,i]
	sd = reconSST.pt.sd[,i]
	val = which(is.na(actual) == FALSE)
	miny = -3.5
	maxy = 3.5

	xpoly = c(c(0:10)[val],rev(c(0:10)[val]))
	ypoly = c(actual[val]-1,rev(actual[val]+1))
	plot(0:10, actual, type="l", xlim = rev(range(0:10)), axes=FALSE, ylim = c(miny,maxy), ylab="", xlab="",col="#0072B2")
	polygon(xpoly,ypoly,col=rgb(0.337, 0.706, 0.914 ,0.25),border = "#56B4E9")

	abline(v=c(0:10),col="grey75")
	abline(h=0,col="grey75")

	points(0:10,actual,pch=20, col="#0072B2")
	points(0:10,reconsst,pch=23,bg="white",col="#D55E00")
	arrows(0:10,reconsst+sd,0:10,reconsst-sd,angle=90,code=3,length=0.05, lwd=0.5,col="#E69F00")

	text(8,3.5,paste("proxy no.",labels[i]))
	text(3,3.5,paste("RSS =",round(sum((actual[val] - reconsst[val])^2),2)))
	text(8,-3, paste("EOF1 =",round(eofs.lim[i,1],2)))
	text(3,-3, paste("EOF2 =",round(eofs.lim[i,2],2)))

	if (is.element((i-nEP),c(1,5,9,13,17)) == TRUE){
		axis(2, at = c(-3,-1.5,0,1.5,3))
	} else {
		axis(2, at = c(-3,-1.5,0,1.5,3), labels = c("","","","",""))
	}

	if (is.element((i-nEP),c(11,12,13,14)) == TRUE){
		axis(1, at = c(0:10))
	} else {
		axis(1, at = c(0:10), labels = c("","","","","","","","","","",""))
	}
}


## ---- FIGURE 8. Single Proxy Reconstructions ---- ############
###########################################################################
years = seq(2,10,2)
prow = 5
dev.new(height=5,width=10)
par(mfrow = c(prow, 1), mar = c(0.5, 2, 0.5, 0.5), oma=c(2,1,1,1))
for (i in 1:length(years)){
	test = read.table(paste("./postfiles/", region, "/PC3/", region, "_MgCa_", season, "_2npcr_", years[i], "ka.txt", sep=""),header=TRUE)
	zfull = rep(NaN,(nx*ny))
	zfull[i1] = as.numeric(unlist(test))
	zmat = matrix(zfull,nrow=nx,ncol=ny)
	image(xgrid,ygrid,zmat,
		ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)), axes=FALSE, 
		col=myPalette1(100), zlim=c(-1.75,1.75)
		)
	mapnames <- map("world2", 
		xlim=c(min(ygrid), max(xgrid)), ylim=c(min(ygrid), max(ygrid)), 
		boundary=TRUE, interior=TRUE, col="black", add=TRUE
	)
	contour(xgrid, ygrid, zmat, add=TRUE, levels = seq(-2,2,0.2), lwd=2, labcex=1.25, col="black", vfont=c("sans serif","bold"))
	axis(2,seq(-8,8,by=4),cex.axis=2)	

	if (i == length(years)) {
         axis(1,c(100,140,180,220,260,300),labels=c("100","140","180","-140","-100","-60"),cex.axis=2)
    }

}

years = seq(2,10,2)
prow = 5
dev.new(height=5,width=10)
par(mfrow = c(prow, 1), mar = c(0.5, 2, 0.5, 0.5), oma=c(2,1,1,1))
for (i in 1:length(years)){
	test = read.table(paste("./postfiles/", region, "/PC3/", region, "_Uk_", season, "_2npcr_", years[i], "ka.txt", sep=""),header=TRUE)
	zfull = rep(NaN,(nx*ny))
	zfull[i1] = as.numeric(unlist(test))
	zmat = matrix(zfull,nrow=nx,ncol=ny)
	image(xgrid,ygrid,zmat,
		ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)), axes=FALSE, 
		col=myPalette1(100), zlim=c(-1.75,1.75)
		)
	mapnames <- map("world2", 
		xlim=c(min(ygrid), max(xgrid)), ylim=c(min(ygrid), max(ygrid)), 
		boundary=TRUE, interior=TRUE, col="black", add=TRUE
	)
	contour(xgrid, ygrid, zmat, add=TRUE, levels = seq(-2,2,0.2), lwd=2, labcex=1.25, col="black", vfont=c("sans serif","bold"))
	axis(2,seq(-8,8,by=4),cex.axis=2)	

	if (i == length(years)) {
         axis(1,c(100,140,180,220,260,300),labels=c("100","140","180","-140","-100","-60"),cex.axis=2)
    }

}



###################################################################
###################################################################
############# APPENDIX FIGURES ####################################
###################################################################
###################################################################


## ---- FIGURE A1. SST EOFs and PCs (ENSO PAPER) ---- #####################
###########################################################################
prow = 3

plotdata = as.data.frame(unname(cbind(recordnopt,eofs.lim[,1:3])))
names(plotdata) = c("No.","Lat","Lon","EOF1","EOF2","EOF3")
library(seqinr)
dblue <- col2alpha("#2066AC",alpha=0.5)
lblue <- col2alpha("#BCDAEA",alpha=0.65)
lred <- col2alpha("#FBC8AF",alpha=0.65)
dred <- col2alpha("#B2182B",alpha=0.5)

dev.new(height=5,width=8)
par(mfrow = c(prow, 1), mar = c(2.5, 2.5, 1, 1))
for (i in 1:prow){
	# PLOT THE FULL FIELD EOF
	scaleEOF = eofs.full[,i]*(lam.full[i]/lam.full[1])
	zfull = rep(NaN,(nx*ny))
	zfull[i1] = scaleEOF 						# add in [index1] if youre plotting India rainfall
	zmat = matrix(zfull,nrow=nx,ncol=ny)
	image(xgrid,ygrid,zmat,
		ylim=range(-15,15),
		xlim=range(95,280),
		xlab="",ylab="", axes=FALSE,
		cex = 1.75,
		main = paste("EOF no.", i, sep=""), 
		col=myPalette1(100),
		zlim=c(-0.12,0.12)
		)
	contour(xgrid, ygrid, zmat, add = TRUE, levels = seq(-0.12,0.12,0.04), lwd = 1.5, labcex=1.25, vfont=c("sans serif","bold"))
	mapnames <- map("world2", 
    xlim=c(min(ygrid), max(xgrid)), 
    ylim=c(min(ygrid), max(ygrid)), 
    boundary=TRUE, ,
    interior=TRUE,
    col="grey50",
    add=TRUE
	)
	axis(2,seq(-10,10,by=5),cex.axis=1.5)

	# PLOT THE LIMITED FIELD EOF BUBBLES
	EOFlim = plotdata[,i+3]
	ratio = lam.lim[1:3]/lam.lim[1]
	EOFrat = abs(EOFlim)*ratio[i]
	radii = sqrt(EOFrat/pi)
	xx = rep(-1,27)
	xx[which(EOFlim > 0)] = 1
	xradii = radii * xx
	print(xradii)
	symbols(x=plotdata$Lon, y=plotdata$Lat, 
		circles=abs(xradii), inches = max(abs(xradii)), 
		bg = ifelse(xradii < -.2 ,dblue, 
		ifelse(xradii >= -.2 & xradii <= 0,lblue,
		ifelse(xradii > 0 & xradii <= 0.2, lred,dred))),
		add=TRUE
	)


	# sEOFlim = EOFlim*(lam.lim[i]/lam.lim[1])
	# rEOFlim = sqrt(sEOFlim/pi)
	# symbols(x=plotdata$Lon, y=plotdata$Lat, 
	# 	circles=abs(sEOFlim), inches= .4*(lam.lim[i]/lam.lim[1]), 
	# 	bg = ifelse(sEOFlim < -.2 ,dblue, 
	# 	ifelse(sEOFlim >= -.2 & sEOFlim <= 0,lblue,
	# 	ifelse(sEOFlim > 0 & sEOFlim <= 0.2, lred,dred))),
	# 	add=TRUE
	# )

	# ADD LONGITUDE AXIS TO LAST PANEL
	if (i == prow) {
         axis(1,seq(100,280,by=20),labels=c("100","120","140","160","180","-160","-140","-120","-100","-80"),cex.axis=1.5)
    }

}

# PLOT OVERLYING PCs
dev.new(height=5,width=3)
par(mfrow = c(prow, 1), mar = c(2.5, 2.5, 1, 2))
for (i in 1:prow){
	plot(1854:2013, scale(pcs.full[,i]), type="l", main=paste("PC no.", i, sep=""), axes=FALSE)
	axis(2, cex.axis=1.5)
	par(new=TRUE)
	plot(1854:2013, scale(pcs.lim[,i]), type="l" ,col="#B2182B", xaxt="n", yaxt="n", xlab="", ylab="", axes=FALSE)
	axis(4, cex.axis=1.5)
	if (i == prow) {
     	axis(1, at = seq(1854,2014,5), cex.axis=1.5)
	}
}



## ---- FIGURE A2. EVS FOR SST FULL/LIM and WINDS (ENSO PAPER) ---- #######
###########################################################################
tot = 8
par(cex.axis = 1.1, cex.lab = 1.5, cex.main = 1.5, las = 1, mar = c(5,5.5,1,1), mfrow = c(1,1), no.readonly = TRUE)
basepar = par()
par(mar = c(5,5.5,3,1))
baseparTitle = par()
dev.new(width=4,height=6)
plot(100 * lam.full[1:tot], type = "l", ylab = "Variance Explained (%)", 
	main = "Eigenvalue Spectra", 
    bty = "n", xaxt = "n", xlab = "", ylim=c(0,100), lwd=2)
abline(v = 1:tot, col = "grey")
lines(100 * lam.lim[1:tot], col="#B2182B", lwd=2)
lines(100 * wlam.full[1:tot], col="#00665E", lwd=2)

points(100 * lam.full[1:tot], type="p", pch=20, col="black", cex=1)
points(1, 100 * lam.full[1], type = "p", pch = 8, col = "black", cex=1.5)
points(3, 100 * lam.full[3], type = "p", pch = 3, col = "black", cex=1.5)
#points(4:tot, 100 * lam.full[4:tot], type="p", pch=20, col="grey50", cex=1.5)

points(100 * lam.lim[1:tot], type="p", pch=20, col="#B2182B", cex=1)
points(1, 100 * lam.lim[1], type = "p", pch = 8, col = "#B2182B", cex=1.5)
points(2, 100 * lam.lim[2], type = "p", pch = 3, col = "#B2182B", cex=1.5)
#points(3:tot, 100 * lam.lim[3:tot], type="p", pch=20, col="grey50", cex=1.5)

points(100 * wlam.full[1:tot], type="p", pch=20, col="#00665E", cex=1)
points(1, 100 * wlam.full[1], type = "p", pch = 8, col = "#00665E", cex=1.5)
points(2, 100 * wlam.full[2], type = "p", pch = 3, col = "#00665E", cex=1.5)
#points(3:tot, 100 * wlam.full[3:tot], type="p", pch=20, col="grey50", cex=1.5)

axis(1, at = seq(1,tot,1))
title(xlab = "Modes")
legend(4, 100, legend = c("Full SST Field", "Limited SST Field", "Full Wind Field"),
	col = c("black", "#B2182B", "#00665E"),
	pch = 19, cex=0.75, bg="white")


## ---- FIGURE A3. EOFs and PCs for WINDS (ENSO PAPER) ---- #######
###################################################################
prow = 6
xgrid <- seq(100,300,by=2.5)
nx <- length(xgrid)
ygrid <- seq(-10,10,by=2.5)
ny <- length(ygrid)
i1 = 1:(nx*ny)

dev.new(height=10,width=8)
par(mfrow = c(prow, 1), mar = c(2.5, 2.5, 1, 1))
for (i in 1:prow){
	scaleEOF = weofs.full[,i]*(wlam.full[i]/wlam.full[1])
	zfull = rep(NaN,(nx*ny))
	zfull[i1] = scaleEOF 						# add in [index1] if youre plotting India rainfall
	zmat = matrix(zfull,nrow=nx,ncol=ny)
	image(xgrid,ygrid,zmat,
		ylim=range(-15,15),
		xlim=range(95,280),
		xlab="",ylab="", axes=FALSE,
		cex = 1.75,
		main = paste("EOF no.", i, sep=""), 
		col=myPalette4(100), horizontal=TRUE,
		zlim=c(-0.12,0.12)
		)
	contour(xgrid, ygrid, zmat, add = TRUE, levels = seq(-0.12,0.12,0.04), lwd = 1.5, labcex=1.25, vfont=c("sans serif","bold"))
	mapnames <- map("world2", 
    xlim=c(min(ygrid), max(xgrid)), 
    ylim=c(min(ygrid), max(ygrid)), 
    boundary=TRUE, ,
    interior=TRUE,
    col="grey50",
    add=TRUE
	)
	axis(2,seq(-10,10,by=5),cex.axis=1.5)

	# ADD LONGITUDE AXIS TO LAST PANEL
	if (i == prow) {
         axis(1,seq(100,280,by=20),labels=c("100","120","140","160","180","-160","-140","-120","-100","-80"),cex.axis=1.5)
    }
}

# PLOT OVERLYING PCs
dev.new(height=10,width=3)
par(mfrow = c(prow, 1), mar = c(2.5, 2.5, 1, 1))
for (i in 1:prow){
	plot(1949:2013, wpcs.full[,i], type="l", main=paste("PC no.", i, sep=""), axes=FALSE)
	if (i == prow) {
     	axis(1, at = seq(1949,2014,5), cex.axis=1.5)
	}
}


## ---- FIGURE A4. CALIBRATION PLOTS (both sst & wind) ---- #######
###################################################################
# Calculate the R2 values
nstasst = 973
scor = c()
for (i in 1:nstasst){
	xcor = cor(actualsst[,i], reconsst[,i])^2
	scor = c(scor, xcor)
}

nstawind = 729
wcor = c()
for (i in 1:nstawind){
	xcor = cor(actualwind[,i], reconwind[,i])^2
	wcor = c(wcor, xcor)
}

# Calculate beta values
sumdiff = apply((actualsst - reconsst)^2, 2, sum)
sumYraw = apply(actualsst^2, 2, sum)
sbeta = 1 - (sumdiff/sumYraw)

sumdiff = apply((actualwind - reconwind)^2, 2, sum)
sumYraw = apply(actualwind^2, 2, sum)
wbeta = 1 - (sumdiff/sumYraw)

xgrid = seq(100,300,by=2)
nx = length(xgrid)
ygrid = seq(-10,10,by=2)
ny = length(ygrid)

wxgrid <- seq(100,300,by=2.5)
nwx = length(wxgrid)
wygrid <- seq(-10,10,by=2.5)
nwy = length(wygrid)

dev.new(width=10, height=4)
par(mfrow = c(2, 1), mar = c(2,2,2,2))
# BETA
zfull = rep(NaN,(nx*ny))
zfull[i1] = sbeta
zmat = matrix(zfull,nrow=nx,ncol=ny)
image(xgrid,ygrid,zmat,
	ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)),
	xlab="",ylab="", main=expression(paste("Calibration ", beta)),
	cex.axis=1.1, yaxt="n", col=rev(myPalette3(100)),
	zlim=c(-1,1)
	)
axis(2, at = c(-10,0,10))
mapnames <- map("world2", 
	xlim=c(min(ygrid), max(xgrid)), ylim=c(min(ygrid), max(ygrid)), 
	boundary=TRUE, interior=TRUE, add=TRUE, lwd=1.5, col="grey50"
)
wfull = as.numeric(wbeta)
wmat = matrix(wfull,nrow=nwx,ncol=nwy)
contour(xgrid, ygrid, zmat, add=TRUE, levels = (0,0.5), lwd=1, labcex=1.25, col="white", vfont=c("sans serif","bold"))
contour(wxgrid, wygrid, wmat, add=TRUE, levels = seq(-1,1,0.1), lwd = 2, labcex=1.5, col="black", vfont=c("sans serif","bold"))

# R^2
zfull = rep(NaN,(nx*ny))
zfull[i1] = scor
zmat = matrix(zfull,nrow=nx,ncol=ny)
image(xgrid,ygrid,zmat,
	ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)),
	xlab="",ylab="", main=expression(paste("Calibration ", R^{2})),
	cex.axis=1.1, yaxt="n", col=rev(myPalette2(100)),
	zlim=c(0,1)
	)
axis(2, at = c(-10,0,10))
mapnames <- map("world2", 
	xlim=c(min(ygrid), max(xgrid)), ylim=c(min(ygrid), max(ygrid)), 
	boundary=TRUE, interior=TRUE, add=TRUE, lwd=1.5, col="grey50"
)
wfull = as.numeric(wcor)
wmat = matrix(wfull,nrow=nwx,ncol=nwy)
contour(xgrid, ygrid, zmat, add=TRUE, levels = c(0,0.5), lwd=1, labcex=1.25, col="white", vfont=c("sans serif","bold"))
contour(wxgrid, wygrid, wmat, add=TRUE, levels = seq(0,1,0.1), lwd = 2, labcex=1.5, col="black", vfont=c("sans serif","bold"))


## ---- FIGURE A5. VERIFICATION PLOTS (both sst & wind) ---- ######
###################################################################
setwd("/Users/emilygill/Documents/University of Colorado/PHD Research/6. SST Reconstruction/REDO Sep2015/")
verbeta.sst = scan(paste("./postfiles/", region, "/PC3/diagnostics/", region, "_", proxy, "_", season, "_beta_v_sst.txt", sep=""))
verR2.sst = scan(paste("./postfiles/", region, "/PC3/diagnostics/", region, "_", proxy, "_", season, "_R2_v_sst.txt", sep=""))

verbeta.wind = scan(paste("./postfiles/", region, "/CCA-PC6/diagnostics/", region, "_", proxy, "_", season, "_beta_v_wind.txt", sep=""))
verR2.wind = scan(paste("./postfiles/", region, "/CCA-PC6/diagnostics/", region, "_", proxy, "_", season, "_R2_v_wind.txt", sep=""))

dev.new(width=10, height=4)
par(mfrow = c(2, 1), mar = c(2,2,2,2))
# BETA
zfull = rep(NaN,(nx*ny))
zfull[i1] = verbeta.sst
zmat = matrix(zfull,nrow=nx,ncol=ny)
image(xgrid,ygrid,zmat,
	ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)),
	xlab="",ylab="", main=expression(paste("Validation ", beta)),
	cex.axis=1.1, yaxt="n", col=rev(myPalette3(100)),
	zlim=c(-1,1)
	)
axis(2, at = c(-10,0,10))
mapnames <- map("world2", 
	xlim=c(min(ygrid), max(xgrid)), ylim=c(min(ygrid), max(ygrid)), 
	boundary=TRUE, interior=TRUE, add=TRUE, lwd=1.5, col="grey50"
)
wfull = as.numeric(verbeta.wind)
wmat = matrix(wfull,nrow=nwx,ncol=nwy)
contour(xgrid, ygrid, zmat, add=TRUE, levels = c(0,0.5), lwd=1, labcex=1.25, col="white", vfont=c("sans serif","bold"))
contour(wxgrid, wygrid, wmat, add=TRUE, levels = seq(-1,1,0.2), lwd = 2, labcex=1.5, col="black", vfont=c("sans serif","bold"))

# R^2
zfull = rep(NaN,(nx*ny))
zfull[i1] = verR2.sst
zmat = matrix(zfull,nrow=nx,ncol=ny)
image(xgrid,ygrid,zmat,
	ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)),
	xlab="",ylab="", main=expression(paste("Validation ", R^{2})),
	cex.axis=1.1, yaxt="n", col=rev(myPalette2(100)),
	zlim=c(0,1)
	)
axis(2, at = c(-10,0,10))
mapnames <- map("world2", 
	xlim=c(min(ygrid), max(xgrid)), ylim=c(min(ygrid), max(ygrid)), 
	boundary=TRUE, interior=TRUE, add=TRUE, lwd=1.5, col="grey50"
)
wfull = as.numeric(verR2.wind)
wmat = matrix(wfull,nrow=nwx,ncol=nwy)
contour(xgrid, ygrid, zmat, add=TRUE, levels = c(0,0.5), lwd=1, labcex=1.25, col="white", vfont=c("sans serif","bold"))
contour(wxgrid, wygrid, wmat, add=TRUE, levels = seq(0,1,0.2), lwd = 2, labcex=1.5, col="black", vfont=c("sans serif","bold"))



## ---- TABLE A1. NINO INDEX VALUES ---- #######
###################################################################

# lat for wpac, nino4, nino34, and nino3
lat1 = which(ygrid == -4)
lat2 = which(ygrid == 4) 

# lat for nino12
n12lat1 = which(ygrid == -10)
n12lat2 = which(ygrid == 0)

#WPAC
wplon1 = which(xgrid == 120)
wplon2 = which(xgrid == 160)

#nino4
n4lon1 = which(xgrid == 160)
n4lon2 = which(xgrid == 210)

#nino34
n34lon1 = which(xgrid == 190)
n34lon2 = which(xgrid == 240)

#nino3
n3lon1 = which(xgrid == 210)
n3lon2 = which(xgrid == 270)

#nino12
n12lon1 = which(xgrid == 270)
n12lon2 = which(xgrid == 280)

table=c()
years = seq(2,10,2)
proxies = c("both","MgCa","Uk")
for (i in 1:length(years)){
	for (j in 1:length(proxies)){
		sst = read.table(paste("./postfiles/EqPac/PC3/EqPac_", proxies[j], "_annual_2npcr_", years[i], "ka.txt", sep=""),header=TRUE)
		sst.sd = read.table(paste("./postfiles/EqPac/PC3/EqPac_", proxies[j], "_annual_2npcr_", years[i], "ka_sd.txt", sep=""),header=TRUE)

			zfull = rep(NaN,(nx*ny))
			zfull[i1] = as.numeric(sst)
			zmat = matrix(zfull, nrow=nx, ncol=ny)

			zsdfull = rep(NaN,(nx*ny))
			zsdfull[i1] = as.numeric(t(sst.sd))
			zsdmat = matrix(zsdfull, nrow=nx, ncol=ny)

			wpac = mean(zmat[wplon1:wplon2,lat1:lat2], na.rm=TRUE)
			nas = dim(which(zsdmat[wplon1:wplon2,lat1:lat2] == "NaN",arr.ind=TRUE))[1]
			n = length(zsdmat[wplon1:wplon2,lat1:lat2]) - nas
			wpacsd = sqrt(sum(zsdmat[wplon1:wplon2,lat1:lat2]^2, na.rm=TRUE))/sqrt(n)
			print(paste(proxies[j], "_wpac_", years[i], "ka = ", round(wpac,2),sep=""))
			print(paste(proxies[j], "_wpac_sd_", years[i], "ka = ", round(wpacsd,2),sep=""))

			nino4 = mean(zmat[n4lon1:n4lon2,lat1:lat2])
			n = length(zsdmat[n4lon1:n4lon2,lat1:lat2])
			nino4sd = sqrt(sum(zsdmat[n4lon1:n4lon2,lat1:lat2]^2, na.rm=TRUE))/sqrt(n)
			print(paste(proxies[j], "_nino4_", years[i], "ka = ", round(nino4,2),sep=""))
			print(paste(proxies[j], "_nino4_sd_", years[i], "ka = ", round(nino4sd,2),sep=""))
		
			nino34 = mean(zmat[n34lon1:n34lon2,lat1:lat2])
			n = length(zsdmat[n34lon1:n34lon2,lat1:lat2])
			nino34sd = sqrt(sum(zsdmat[n34lon1:n34lon2,lat1:lat2]^2, na.rm=TRUE))/sqrt(n)
			print(paste(proxies[j], "_nino34_", years[i], "ka = ", round(nino34,2),sep=""))
			print(paste(proxies[j], "_nino34_sd_", years[i], "ka = ", round(nino34sd,2),sep=""))

			nino3 = mean(zmat[n3lon1:n3lon2,lat1:lat2])
			nino3sd = mean(zsdmat[n3lon1:n3lon2,lat1:lat2])
			print(paste(proxies[j], "_nino3_", years[i], "ka = ", round(nino3,2),sep=""))
			print(paste(proxies[j], "_nino3_sd_", years[i], "ka = ", round(nino3sd,2),sep=""))

			nino12 = mean(zmat[n12lon1:n12lon2,n12lat1:n12lat2], na.rm=TRUE)
			nas = dim(which(zsdmat[n12lon1:n12lon2,n12lat1:n12lat2] == "NaN",arr.ind=TRUE))[1]
			n = length(zsdmat[n12lon1:n12lon2,n12lat1:n12lat2]) - nas
			nino12sd = sqrt(sum(zsdmat[n12lon1:n12lon2,n12lat1:n12lat2]^2, na.rm=TRUE))/sqrt(n)
			print(paste(proxies[j], "_nino12_", years[i], "ka = ", round(nino12,2),sep=""))
			print(paste(proxies[j], "_nino12_sd_", years[i], "ka = ", round(nino12sd,2),sep=""))

			tni = nino12 - nino4
			tnisd = sqrt(nino12sd^2 + nino4sd^2)
			print(paste(proxies[j], "_tni_", years[i], "ka = ", round(tni,2),sep=""))
			print(paste(proxies[j], "_tni_sd_", years[i], "ka = ", round(tnisd,2),sep=""))

			wtni = nino12 - wpac
			wtnisd = sqrt(nino12sd^2 + wpacsd^2)
			print(paste(proxies[j], "_wtni_", years[i], "ka = ", round(wtni,2),sep=""))
			print(paste(proxies[j], "_wtni_sd_", years[i], "ka = ", round(wtnisd,2),sep=""))

			results = cbind(wpac, wpacsd, nino4, nino4sd, nino34, nino34sd, nino3, nino3sd, nino12, nino12sd, tni, tnisd, wtni, wtnisd)
			table = rbind(table, results)
	}
}

	# write.table(table, 
	# 	file="./postfiles/EqPac/PC3/diagnostics/ninotable.txt",
	# 	sep="\t")
 


































