read.contemp <- function(region, season){

		ntime <- 160		# number of years
		nclim <- 1
	if (region == "EqPac"){
		xgrid <- seq(100,300,by=2)
		nx <- length(xgrid) # nrows = 101
		ygrid <- seq(-10,10,by=2)
		ny <- length(ygrid) # ncols = 11

		nsta <- nx*ny
		if (season == "annual"){
			data <- readBin("datafiles/1854-2014_MAYtoAPRavg_ERSSST.sst.r4",what="numeric",n=(nx*ny*ntime),size=4,endian="swap")	
			clim <- readBin("datafiles/1981-2010_MAYtoAPRavg_ERSSST.sst-CLIM.r4",what="numeric",n=(nx*ny*nclim),size=4,endian="swap")
		} else if (season == "winter"){
			data <- readBin("datafiles/1854-2014_DECtoFEBavg_ERSSST.sst.r4",what="numeric",n=(nx*ny*ntime),size=4,endian="swap")
			clim <- readBin("datafiles/1981-2010_DECtoFEBavg_ERSSST.sst-CLIM.r4",what="numeric",n=(nx*ny*nclim),size=4,endian="swap")
		} else if (season == "summer"){
			data <- readBin("datafiles/1854-2013_JUNtoAUGavg_ERSSST.sst.r4",what="numeric",n=(nx*ny*ntime),size=4,endian="swap")
			clim <- readBin("datafiles/1981-2010_JUNtoAUGavg_ERSSST.sst-CLIM.r4",what="numeric",n=(nx*ny*nclim),size=4,endian="swap")
		}		
	} else if (region == "IO") {
		xgrid <- seq(30,106,by=2)
		nx <- length(xgrid) # nrows = 39
		ygrid <- seq(-10,28,by=2)
		ny <- length(ygrid) # ncols = 20

		nsta <- nx*ny
		if (season == "annual"){
			data <- readBin("1854-2014_MAYtoAPRavg_ERSSST.sst.IO.r4",what="numeric",n=(nx*ny*ntime),size=4,endian="swap")	
			clim <- readBin("1981-2010_MAYtoAPRavg_ERSSST.sst.IO-CLIM.r4",what="numeric",n=(nx*ny*nclim),size=4,endian="swap")
		} else if (season == "winter"){
			data <- readBin("1854-2014_NOVtoAPRavg_ERSSST.sst.IO.r4",what="numeric",n=(nx*ny*ntime),size=4,endian="swap")
			clim <- readBin("1981-2010_NOVtoAPRavg_ERSSST.sst.IO-CLIM.r4",what="numeric",n=(nx*ny*nclim),size=4,endian="swap")
		} else if (season == "summer"){
			data <- readBin("datafiles/1854-2014_JUNtoSEPavg_ERSSST.sst.IO.r4",what="numeric",n=(nx*ny*ntime),size=4,endian="swap")
			clim <- readBin("datafiles/1981-2010_JUNtoSEPavg_ERSSST.sst.IO-CLIM.r4",what="numeric",n=(nx*ny*nclim),size=4,endian="swap")
		}	

	} 
	else if (region == "fullslab"){
		xgrid <- seq(30,300,by=2)
		nx <- length(xgrid) # nrows = 39
		ygrid <- seq(-10,28,by=2)
		ny <- length(ygrid) # ncols = 20

		nsta <- nx*ny
		
		if (season == "annual"){
			data <- readBin("1854-2014_MAYtoAPRavg_ERSSST.sst.fullslab.r4",what="numeric",n=(nx*ny*ntime),size=4,endian="swap")
			clim <- readBin("1981-2010_MAYtoAPRavg_ERSSST.sst.fullslab-CLIM.r4",what="numeric",n=(nx*ny*nclim),size=4,endian="swap")	
	
		} else if (season == "summer"){
			data <- readBin("1854-2014_JUNtoSEPavg_ERSSST.sst.30-300.10S-28N.r4",what="numeric",n=(nx*ny*ntime),size=4,endian="swap")
			clim <- readBin("1981-2010_JUNtoSEPavg_ERSSST.sst.30-300.10S-28N-CLIM.r4",what="numeric",n=(nx*ny*nclim),size=4,endian="swap")	
		}
		

	}

	data[is.nan(data)] <- -999
	BINresults <- BINxygrid(data,xgrid,ygrid,ntime,type="special", undefined=-999, condition="greater")
	results = list(xgrid = BINresults$xgrid, ygrid = BINresults$ygrid, grid = BINresults$grid,
		index1 = BINresults$ind, nsta = BINresults$nsta, data = BINresults$data, clim=clim)

	return(results)

}



read.altdata <- function(var,region){
	if (var == "winds"){
		if (region == "EqPac"){	
			xgrid <- seq(100,300,by=2)
			nx <- length(xgrid) # nrows = 101
			ygrid <- seq(-10,10,by=2)
			ny <- length(ygrid) # ncols = 11

			ntime <- 66		# number of years
			nsta <- nx*ny
			udata <- readBin("1949-2015_JUNtoSEPavg_NCEPNCAR.u850.Pac.r4",
				what="numeric",n=(nx*ny*ntime),size=4,endian="swap")
			vdata <- readBin("1949-2015_JUNtoSEPavg_NCEPNCAR.v850.Pac.r4",
				what="numeric",n=(nx*ny*ntime),size=4,endian="swap")
		} else if (region == "IO"){
			xgrid <- seq(30,105,by=2.5)
			nx <- length(xgrid) # nrows = 31
			ygrid <- seq(-10,30,by=2.5)
			ny <- length(ygrid) # ncols = 17

			ntime <- 66		# number of years
			nsta <- nx*ny
			udata <- readBin("1949-2015_JUNtoSEPavg_NCEPNCAR.u850.IO.r4",
				what="numeric",n=(nx*ny*ntime),size=4,endian="swap")
			vdata <- readBin("1949-2015_JUNtoSEPavg_NCEPNCAR.v850.IO.r4",
				what="numeric",n=(nx*ny*ntime),size=4,endian="swap")
		}

		data[is.nan(data)] <- -999
		uwind <- BINxygrid(udata,xgrid,ygrid,ntime,type="special", undefined=-999, condition="greater")
		vwind <- BINxygrid(vdata,xgrid,ygrid,ntime,type="special", undefined=-999, condition="greater")
		results = list(uwind, vwind)
		return(results)
	}

}


gcv.recon <- function(EP, WP, IO, EPdeg, WPdeg, IOdeg, base){
	nEP = length(EP)			# number of east Pacific records
	nWP = length(WP)			# number of west Pacific records
	nIO = length(IO)			# number of Indian Ocean records
	# ntot = nEP + nWP + nIO

	EP.smooth = list()
	for (i in 1:nEP){
		N = dim(EP[[i]])[1]
		zfit = locfit(EP[[i]]$SST2 ~ EP[[i]]$age, deg=EPdeg, kern="bisq", scale=TRUE)

		if (round(EP[[i]]$age[N]) > 10){
			maxyr = 10
		} else {
			maxyr = round(EP[[i]]$age[N])
		}
		minyr = floor(EP[[i]]$age[1])

		if (base == "TRUE"){
			aa = 10
			bb = 0
		}
		yestbase = predict(zfit,seq(bb,aa,1))

		xsmo = rep(NaN,11)
		xsmo[(minyr+1):(maxyr+1)] = yestbase[(minyr+1):(maxyr+1)]

		xEP.smooth = matrix(cbind(seq(bb,aa,1), yestbase, xsmo), 11, 3)

		#xEP.smooth = matrix(cbind(seq(minyr,maxyr,1),yestbase),length(yestbase),2)
		EP.smooth[[i]] = xEP.smooth
	}


	WP.smooth = list()
	for (i in 1:nWP){
		N = dim(WP[[i]])[1]
		zfit = locfit(WP[[i]]$SST2 ~ WP[[i]]$age, deg=WPdeg, kern="bisq", scale=TRUE)

		if (round(WP[[i]]$age[N]) > 10){
			maxyr = 10
		} else {
			maxyr = round(WP[[i]]$age[N])
		}
		minyr = floor(WP[[i]]$age[1])

		if (base == "TRUE"){
			aa = 10
			bb = 0
		}
		yestbase = predict(zfit,seq(bb,aa,1))

		xsmo = rep(NaN,11)
		xsmo[(minyr+1):(maxyr+1)] = yestbase[(minyr+1):(maxyr+1)]

		xWP.smooth = matrix(cbind(seq(bb,aa,1), yestbase, xsmo), 11, 3)

		#xWP.smooth = matrix(cbind(seq(minyr,maxyr,1),yest),length(yest),2)
		WP.smooth[[i]] = xWP.smooth
	}

	IO.smooth = list()
	for (i in 1:nIO){
		N = dim(IO[[i]])[1]
		zfit = locfit(IO[[i]]$SST1 ~ IO[[i]]$age, deg=IOdeg, kern="bisq", scale=TRUE)

		if (round(IO[[i]]$age[N]) > 10){
			maxyr = 10
		} else {
			maxyr = round(IO[[i]]$age[N])
		}
		minyr = floor(IO[[i]]$age[1])

		if (base == "TRUE"){
			aa = 10
			bb = 0
		}
		yestbase = predict(zfit,seq(bb,aa,1))

		xsmo = rep(NaN,11)
		xsmo[(minyr+1):(maxyr+1)] = yestbase[(minyr+1):(maxyr+1)]

		xIO.smooth = matrix(cbind(seq(bb,aa,1), yestbase, xsmo), 11, 3)

		#xIO.smooth = matrix(cbind(seq(minyr,maxyr,1),yest),length(yest),2)
		IO.smooth[[i]] = xIO.smooth
	}

	results = list(EP = EP.smooth, WP = WP.smooth, IO = IO.smooth)
	return(results)

}


plot.recon <- function(raw, smooth, map, region, proxy, coord){
	n = length(raw)

	if (region == "EP"){
		rawsst = raw$SST2
		yrange = c(22,30)
	} else if (region == "WP"){
		rawsst = raw$SST2
		yrange = c(25,33)
	}

	# xmin = c()
	# xmax = c()
	# ymin = c()
	# ymax = c()
	# for (i in 1:n){
	# 	xmin = c(xmin,min(raw[[i]]$age))
	# 	xmax = c(xmax,max(raw[[i]]$age))

	# 	if (region == "EP"){
	# 		ymin = c(ymin,min(raw[[i]]$SST3))
	# 		ymax = c(ymax,max(raw[[i]]$SST3))
	# 	} else if (region == "WP"){
	# 		ymin = c(ymin,min(raw[[i]]$SST3))
	# 		ymax = c(ymax,max(raw[[i]]$SST3))
	# 	}
	# }
	# xrange = c(sort(xmin)[1], sort(xmax)[n]) 
	# xrange[1] = 0
	# yrange = c(sort(ymin)[1], sort(ymax)[n]) 
	# yrange = c(floor(yrange[1]), ceiling(yrange[2]))

	par(mfrow = c(1, 1), mar = c(0.75, 0.75, 0.75, 0.75), oma = c(4, 4, 3, 1))		
	# cbPalette <- c("#D55E00","#E69F00", "#56B4E9", "#009E73", 
	# 	"#F0E442", "#0072B2", "#CC79A7","#999999",
	# 	"#CC0000","#9900CC","#66FF66","#996633","#838B9B")

	# cbPalette <- c("#FF0000FF","#FF6D00FF","#FFDB00FF","#B6FF00FF",
	# 	"#49FF00FF", "#00FF24FF","#00FF92FF","#00FFFFFF","#0092FFFF",
	# 	"#0024FFFF", "#4900FFFF", "#B600FFFF", "#FF00DBFF", "#FF006DFF")

	# cbPalette <- rev(c("#5E4EA1", "#437BB6", "#54A6B0", "#7DCAA4", "#AFDEA3", 
	# 	"#DDF199", "#F5FBB0", "#FFF3AA", "#FED884", "#FDB164", "#F7824B", 
	# 	"#E55849", "#C8334C", "#9E0041"))

	spectral <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
	cbPalette = rev(spectral(n))

	
	## PLOT THE RAW VALUES
	for (i in 1:(n)){
		if (i == 1){
			if (region == "EP"){
				plot(raw[[i]]$age, raw[[i]]$SST2, type="l", xlim=rev(c(0,10)), ylim=yrange, axes=FALSE, col=cbPalette[1], lwd=3)
			} else if (region == "WP"){
				plot(raw[[i]]$age, raw[[i]]$SST2, type="l", xlim=rev(c(0,10)), ylim=yrange, axes=FALSE, col=cbPalette[1], lwd=3)
			}
			# lines(xxs$X, xxs$Y, lty=2, col=cbPalette[1])
			#text(0,head(xxs$Y, 1)+0.25,i, col=cbPalette[1],)
			#mtext(i, side=3, line=1, cex=1, at=head(xxs$Y, 1))
		} else {
			if (region == "EP"){
				lines(raw[[i]]$age, raw[[i]]$SST2, col=cbPalette[i],lwd=3)
			} else if (region == "WP"){
				lines(raw[[i]]$age, raw[[i]]$SST2, col=cbPalette[i],lwd=3)
			}
			# lines(xxs$X, xxs$Y, col=cbPalette[i], lty=2)
			#text(0,head(xxs$Y, 1)+0.25, i, col=cbPalette[i])
		}

		if (i == n){
			axis(1, at=seq(10, 0, by = -1), las=1, cex.axis=1.5)
			axis(2, at=seq(yrange[1], yrange[2], by = 1), cex.axis=1.5)
			mtext("ka BP", side=1, line=3, cex=1)
			mtext(~degree~C, side=2, line=3, cex=1)
		}
	}

	## PLOT THE SMOOTH VALUES
	dev.new()
	par(mfrow = c(1, 1), mar = c(0.75, 0.75, 0.75, 0.75), oma = c(4, 4, 3, 1))		

	for (i in 1:(n)){
		xxs = data.frame(smooth[i])
		names(xxs) = c("X","Y")

		if (i == 1){
			plot(xxs$X, xxs$Y, type="l", lty=2, xlim=rev(c(0,10)), ylim=yrange, axes=FALSE, col=cbPalette[1], lwd=3)
#			lines(xxs$X, xxs$Y, lty=2, col=cbPalette[1])

		} else {
			lines(xxs$X, xxs$Y, lty=2, col=cbPalette[i],lwd=3)
#			lines(xxs$X, xxs$Y, col=cbPalette[i], lty=2)

		}

		if (i == n){
			axis(1, at=seq(10, 0, by = -1), las=1, cex.axis=1.5)
			axis(2, at=seq(yrange[1], yrange[2], by = 1), cex.axis=1.5)
			mtext("ka BP", side=1, line=3, cex=1)
			mtext(~degree~C, side=2, line=3, cex=1)
		}
	}


	if (map == TRUE){
		dev.new()
		par(mfrow = c(1, 1), mar = c(0.75, 0.75, 0.75, 0.75), oma = c(4, 4, 3, 1))		
		
		# Plot boundaries for East Pacific and change coordinates to East (negative)
		if (region == "EP"){
			xgrid = seq(-100,-70,2)
			nx = length(xgrid)
			ygrid = seq(-12,12,2)
			ny = length(ygrid)	

			coord[,2] = -1 * (360 - coord[,2]) 
		}

		# Plot boundaries for West Pacific
		if (region == "WP"){
			xgrid = seq(100,160,2)
			nx = length(xgrid)
			ygrid = seq(-12,12,2)
			ny = length(ygrid)	
		}

		# Plot boundaries for East Pacific and change coordinates to East (negative)
		if (region == "IO"){
			dev.new(width=10, height=7)
			par(mfrow = c(1, 1), mar = c(0.75, 0.75, 0.75, 0.75), oma = c(4, 4, 3, 1))		

			xgrid = seq(30,110,10)
			nx = length(xgrid)
			ygrid = seq(-15,30,5)
			ny = length(ygrid)	 
		}

		# Add countries and reconstruction points
		image(xgrid,ygrid,matrix(rep(NaN,(nx*ny)),nx,ny),
			xlab="",ylab="", cex.axis=1.25)
		
		if (region == "IO"){
			abline(h=c(-10), lty=2, col="grey75")
		} else {
			abline(h=c(-10,10), lty=2, col="grey75")		
		}

		
		world(ylim=range(min(ygrid),max(ygrid)),
			xlim=range(min(xgrid),max(xgrid)),
			add=TRUE,
			fill=TRUE,
			border="darkgrey")
		

		if (proxy == "MgCa"){
			if (region == "EP"){
				labels = c(22:26)  #FOR EPAC Mg/Ca
				pos = c(1,2,2,2,2)  #FOR EPAC Mg/Ca
			} else {
				labels = c(1:14) #FOR WPAC Mg/Ca 
				pos = c(2,2,1,1,2,2,1,1,1,1,1,1,1,1)  #FOR WPAC Mg/Ca				
			}
		} else {
			if (region == "EP"){
				labels = c(27:39) #FOR EPAC Uk
				pos = c(2,1,2,1,1,1,1,1,1,1,1,1,1) # FOR EPAC Uk							
			} else {
				labels = c(15:21) #FOR WPAC Uk
				pos = c(4,4,2,1,4,4,1) #FOR WPAC UK
			}
		}

		#labels = c(30:34) #FOR IO Mg/Ca
		#pos = c(1,2,2,1,1) #FOR IO Mg/Ca

		#labels = c(35:43)  #FOR IO Uk
		#pos = c(2,3,1,2,2,1,1,4,2)

		for (i in 1:n){
			points(coord[i,2], coord[i,1], pch=20, cex=7, col=cbPalette[i])
			# if (i == 2){
			# 	col = cbPalette[2]
			# } else {
			# 	col = "black"
			# }
			# text(coord[i,2], coord[i,1], labels=labels[i], pos=pos[i], col=col, cex=1.5, lwd=1.5)
			text(coord[i,2], coord[i,1], labels=labels[i], col="black", cex=1.75, lwd=1.5)
			#text(coord[i,2], coord[i,1], labels=labels[i], pos=pos[i], col="grey50", cex=2, lwd=1.5)
		}
		
	} 

}


scale.recon <- function(smoodat, sites, clim){
	sst.mat = matrix(NaN,length(seq(0:10)),(sites))
	sst.mat.a = sst.mat 		# Initialize matrix for anomalies (using climatology)
	sst.mat.base = sst.mat 		# Initialize matrix for using base 0 ka as zero
	sst.mat.baselim = sst.mat 		# Initialize matrix for using base 0 ka as zero

	for (i in 1:sites){
		xmu = clim[i]

		xi = smoodat[[i]][,1]
		#paleomu = mean(smoodat[[i]][,2], na.rm=TRUE)
		#paleosig = sd(smoodat[[i]][,2], na.rm=TRUE)
		sst.mat[xi+1,i] = smoodat[[i]][,3]
		sst.mat.a[xi+1,i] =  (smoodat[[i]][,3] - xmu )
		sst.mat.base[xi+1,i] = smoodat[[i]][,2] - smoodat[[i]][1,2]
		sst.mat.baselim[xi+1,i] =  smoodat[[i]][,3] - smoodat[[i]][1,2]

	}

	results = list(unscaled = sst.mat, anom = sst.mat.a, base = sst.mat.base, baselim=sst.mat.baselim)
	return(results)

}



plot.PCA.test <- function(){

	## Make a function that will plot years in comparison and spit out the spatial R2 of the model
}


monthly.imd <- function(start,stop,type){
	source("BINxygrid.R")

	xgrid <- seq(66.5,100.5,by=1)
	nx <- length(xgrid) # nrows = 101
	ygrid <- seq(6.5,38.5,by=1)
	ny <- length(ygrid) # ncols = 11

	ntime <- 37986		# number of years
	nsta <- nx*ny
	data <- readBin("1901-2004_daily_IMDrain.r4",what="numeric",n=(nx*ny*ntime),size=4,endian="swap")
	results <- BINxygrid(data,xgrid,ygrid,ntime,type="special", undefined=-1, condition="greater")
	data.un <- results$data
	leap <- seq(425,1461*104,1461)
	data.un <- data.un[-leap,]
	grid <- results$grid
	index1 <- results$ind
	nsta <- results$nsta

	starts = seq(start,ntime,365)
	stops = seq(stop,ntime,365)

	xdata = matrix(NA,ntime/365,nsta)

	if (type == "sestotal"){
		for (i in 1:(ntime/365)){
			xdata[i,] = colSums(data.un[starts[i]:stops[i],])
		}
	} else if (type == "sesavg"){
		for (i in 1:(ntime/365)){
			xdata[i,] = colMeans(data.un[starts[i]:stops[i],])	
		}
	}

	results = list(data = xdata, nsta = results$nsta, index1 = results$ind, grid = results$grid)
	return(results)
}


plot.validate <- function(npc, region, proxy, season, index1, xgrid, ygrid, reconpt, zzpt, mat){
	nx = length(xgrid)
	ny = length(ygrid)
	nreconpt = dim(reconpt)[1]
	years = seq(0,10,1)
	reconSST.mat = c()
	reconSST.mat.sd = c()

	prow = 11
	dev.new(length=5,width=8)
	par(mfrow = c(prow, 2), mar = c(0.5, 2, 0.5, 0.5))	
	for (i in 1:length(years)){
			test = read.table(paste("./postfiles/", region, "/PC", npc, "/", region, "_", proxy, "_", season, "_", npcr, "npcr", "_", i-1, "ka.txt", sep=""),header=TRUE)
			reconSST.mat = rbind(reconSST.mat,test)

			sd = read.table(paste("./postfiles/", region, "/PC", npc, "/", region, "_", proxy, "_", season, "_", npcr, "npcr", "_", i-1, "ka_sd.txt", sep=""),header=TRUE)
			reconSST.mat.sd = rbind(reconSST.mat.sd,t(sd))
			
			zfull = rep(NaN,(nx*ny))
			zfull[index1] = as.numeric(unlist(test))
			zmat = matrix(zfull,nrow=nx,ncol=ny)
			image.plot(xgrid,ygrid,zmat,
				ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)), axes=FALSE, 
				col=myPalette1(100), zlim=c(-2,2)
				)
			contour(xgrid, ygrid, zmat, add = TRUE, nlev = 6, lwd = 1.2)
			mapnames <- map("world2", 
				xlim=c(min(ygrid), max(xgrid)), ylim=c(min(ygrid), max(ygrid)), 
				boundary=TRUE, interior=TRUE, col="black", add=TRUE
			)
			text(294,-6,paste(i-1,"ka",sep=""),col="red")

			zfull = rep(NaN,(nx*ny))
			zfull[index1] = as.numeric(unlist(sd))
			zmat = matrix(zfull,nrow=nx,ncol=ny)
			image.plot(xgrid,ygrid,zmat,
				ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)), axes=FALSE,
				col=myPalette1(100), zlim=c(-2,2)
				)
			contour(xgrid, ygrid, zmat, add = TRUE, nlev = 6, lwd = 1.2)
			mapnames <- map("world2", 
				xlim=c(min(ygrid), max(xgrid)), ylim=c(min(ygrid), max(ygrid)), 
				boundary=TRUE, interior=TRUE, col="black", add=TRUE
			)
	}

	reconSST.pt = reconSST.mat[,zzpt]
	reconSST.pt.sd = reconSST.mat.sd[,zzpt]

	prow = ceiling(nreconpt/5)
	dev.new(width=12,height=8)
	par(mfrow = c(prow,5), mar = c(0.75, 0.75, 0.75, 0.75), oma = c(2, 2, 1, 1))		
	for (i in 1:nreconpt){
		actual = mat[,i]
		reconsst = reconSST.pt[,i]
		sd = reconSST.pt.sd[,i]
		min = min(actual,reconsst-sd)
		max = max(actual,reconsst+sd)
		val = which(is.na(actual) == FALSE)
		beta = 1-((sum((actual[val]-reconsst[val])^2))/(sum(actual[val]^2)))
		beta = round(beta,2)
		R2 = round(cor(actual[val], reconsst[val]),2)

		plot(0:10, actual, type="l", xlim = rev(range(0:10)), ylim = c(-2,2), axes=FALSE)
		abline(v=c(0:10),col="grey90")
		abline(h=0,col="grey90")
		points(0:10,actual,pch=20)
		points(0:10,reconsst,pch=23,bg="white")
		arrows(0:10,reconsst+sd,0:10,reconsst-sd,angle=90,code=3,length=0.05, lwd=0.5,col="grey50")
		text(8,1,round(beta,2))
		text(5,1,round(R2,2))
		text(2,1,round(sum((actual[val] - reconsst[val])^2),2))
		if (reconpt[i,2] < 180){
			text(9,3,labels = paste("(",reconpt[i,1],",",reconpt[i,2],")",sep=""))	
		} else {
			text(9,3,labels = paste("(",reconpt[i,1],",",reconpt[i,2]-360,")",sep=""))
		}
		
		if (is.element(i,c(1,6,11,16,21,26)) == TRUE){
			axis(2, at = c(-4:4))
		}

		if (is.element(i,seq(nreconpt-4,nreconpt,1)) ==TRUE) {
			axis(1, at = c(0:10))
		}
	}



}








plot.validate2 <- function(var, npc, region, proxy, season, index1, xgrid, ygrid, reconpt, zzpt, mat){
	nx = length(xgrid)
	ny = length(ygrid)
	nreconpt = dim(reconpt)[1]
	years = seq(2,10,2)
	reconSST.mat = c()
	reconSST.mat.sd = c()

	# prow = length(years)+1
	prow = 5
	dev.new(height=5,width=10)
	par(mfrow = c(prow, 1), mar = c(0.5, 2, 0.5, 0.5), oma=c(2,1,1,1))

	# readpresent = read.table(paste("./postfiles/", region, "/PC", npc, "/", region, "_", proxy, "_", season, "_", npcr, "npcr", "_10ka_1901-2013.txt", sep=""),header=TRUE)
	# present = apply(readpresent,2,mean)

	# zfull = rep(NaN,(nx*ny))
	# zfull[index1] = as.numeric(unlist(present))
	# zmat = matrix(zfull,nrow=nx,ncol=ny)
	# image(xgrid,ygrid,zmat,
	# 	ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)), axes=FALSE, 
	# 	col=myPalette7(100), zlim=c(-1.5,0)
	# 	)
	# contour(xgrid, ygrid, zmat, add = TRUE, levels = seq(-3,3,by=0.2), lwd = 1.2, labcex=1.25)
	# mapnames <- map("world2", 
	# 	xlim=c(min(ygrid), max(xgrid)), ylim=c(min(ygrid), max(ygrid)), 
	# 	boundary=TRUE, interior=TRUE, col="black", add=TRUE
	# )
	# 	axis(2,seq(-8,8,by=4),cex.axis=2)

	for (i in 1:length(years)){
			
			if (var == "wind"){
				test = read.table(paste("./postfiles/", region, "/CCA-PC", npc, "/", region, "_", proxy, "_", season, "_6npc_", years[i], "ka_u.txt", sep=""),header=TRUE)	
				colpalette = myPalette4(100)
				zlimit = c(-12,12)
				contlev = seq(-12,12,by=2)
			} else if (var == "sst"){
				test = read.table(paste("./postfiles/", region, "/PC", npc, "/", region, "_", proxy, "_", season, "_", npcr, "npcr", "_", years[i], "ka.txt", sep=""),header=TRUE)
				colpalette = myPalette1(100)
				zlimit = c(-2.5,2.5)
				contlev = seq(3,3,by=0.2)
			}
			
			reconSST.mat = rbind(reconSST.mat,test)

			# sd = read.table(paste("./postfiles/", region, "/PC", npc, "/", region, "_", proxy, "_", season, "_", npcr, "npcr", "_", years[i], "ka_sd.txt", sep=""),header=TRUE)
			# reconSST.mat.sd = rbind(reconSST.mat.sd,t(sd))
			
			zfull = rep(NaN,(nx*ny))
			zfull[index1] = as.numeric(unlist(test))
			zmat = matrix(zfull,nrow=nx,ncol=ny)
			image(xgrid,ygrid,zmat,
				ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)), axes=FALSE, 
				col=colpalette, zlim=zlimit
				)
			contour(xgrid, ygrid, zmat, add = TRUE, levels = contlev, lwd = 1.2, labcex=1.25)
			mapnames <- map("world2", 
				xlim=c(min(ygrid), max(xgrid)), ylim=c(min(ygrid), max(ygrid)), 
				boundary=TRUE, interior=TRUE, col="black", add=TRUE
			)
				axis(2,seq(-8,8,by=4),cex.axis=2)

			if (i == length(years)) {
		         axis(1,c(100,140,180,220,260,300),labels=c("100","140","180","-140","-100","-60"),cex.axis=2)
		    }
			#text(294,-6,paste(i-1,"ka",sep=""),col="red")

			# zfull = rep(NaN,(nx*ny))
			# zfull[index1] = as.numeric(unlist(sd))
			# zmat = matrix(zfull,nrow=nx,ncol=ny)
			# image.plot(xgrid,ygrid,zmat,
			# 	ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)), axes=FALSE,
			# 	col=myPalette1(100), zlim=c(-2,2)
			# 	)
			# contour(xgrid, ygrid, zmat, add = TRUE, nlev = 6, lwd = 1.2)
			# mapnames <- map("world2", 
			# 	xlim=c(min(ygrid), max(xgrid)), ylim=c(min(ygrid), max(ygrid)), 
			# 	boundary=TRUE, interior=TRUE, col="black", add=TRUE
			# )
	}

	# reconSST.pt = reconSST.mat[,zzpt]
	# reconSST.pt.sd = reconSST.mat.sd[,zzpt]

	# prow = round(nreconpt/3)
	# dev.new()
	# par(mfrow = c(prow,3), mar = c(0.75, 0.75, 0.75, 0.75), oma = c(2, 2, 1, 1))		
	# for (i in 1:nreconpt){
	# 	actual = mat[,i]
	# 	reconsst = reconSST.pt[,i]
	# 	sd = reconSST.pt.sd[,i]
	# 	min = min(actual,reconsst-sd)
	# 	max = max(actual,reconsst+sd)
	# 	plot(0:10, actual, type="l", xlim = rev(range(0:10)), ylim = c(-2,2), axes=FALSE)
	# 	abline(v=c(0:10),col="grey90")
	# 	abline(h=0,col="grey90")
	# 	points(0:10,actual,pch=20)
	# 	points(0:10,reconsst,pch=23,bg="white")
	# 	arrows(0:10,reconsst+sd,0:10,reconsst-sd,angle=90,code=3,length=0.05, lwd=0.5,col="grey50")

	# 	if (reconpt[i,2] < 180){
	# 		text(9,3,labels = paste("(",reconpt[i,1],",",reconpt[i,2],")",sep=""))	
	# 	} else {
	# 		text(9,3,labels = paste("(",reconpt[i,1],",",reconpt[i,2]-360,")",sep=""))
	# 	}
		
	# 	if (is.element(i,c(1,4,7,10,13)) == TRUE){
	# 		axis(2, at = c(-4:4))
	# 	}

	# 	if (is.element(i,seq(nreconpt-2,nreconpt,1)) ==TRUE) {
	# 		axis(1, at = c(0:10))
	# 	}
	# }



}


