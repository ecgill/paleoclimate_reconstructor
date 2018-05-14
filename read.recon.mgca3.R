read.recon.mgca3 <- function(plot){
	
	###################################################################################################
	## Read in the East Pacific SST Reconstructions
	dir <- "/Users/emilygill/Documents/University of Colorado/PHD Research/6. SST Reconstruction/cores/MgCa/EP/"

	# *** This one is the same location as ODP-1242, but shows less of a cooling trend. Find out whats different.
	# EP1 <- c(8,276)			#Benway_ME0005A-43JC <- c(7.9,-83.6) = 276.4
	# sEP1 <- data.frame(matrix(scan(paste(dir,"ME0005A-43JC_SST-age_MgCa.txt",sep="")),ncol=4,byrow=TRUE))
	# names(sEP1) <- c("depth","age","var","SST1")

	EP1 <- c(8,278)			#Benway_ODP-1242 <- c(7.9,-83.6) = 276.4
	sEP1 <- data.frame(matrix(scan(paste(dir,"ODP-1242_SST-age_MgCa.txt",sep="")),ncol=4,byrow=TRUE))
	names(sEP1) <- c("depth","age","var","SST1")

	EP2 <- c(0,268)			#Lea_TR163-22 <- c(0.5,-92.4) = 267.6
	sEP2 <- data.frame(matrix(scan(paste(dir,"TR163-22_SST-age_MgCa.txt",sep="")),ncol=4,byrow=TRUE))
	names(sEP2) <- c("depth","age","var","SST1")

	# EP3 <- c(2,268)			#Lea_TR163-19 <- c(2.3,-90.95) = 269.1
	# sEP3 <- data.frame(matrix(scan(paste(dir,"TR163-22_SST-age_MgCa.txt",sep="")),ncol=4,byrow=TRUE))
	# names(sEP3) <- c("depth","age","var","SST1")

	EP3 <- c(0,274) 		#Pena_ODP-1240 <- c(0,-86.4) = 273.6
	sEP3 <- data.frame(matrix(scan(paste(dir,"ODP-1240_SST-age_MgCa.txt",sep="")),ncol=4,byrow=TRUE))
	names(sEP3) <- c("depth","age","var","SST1")

	EP4 <- c(-2,270) 		#Koutavas2006_V21.30 <- c(-1.2,-89.7) = 270.3
	sEP4 <- data.frame(matrix(scan(paste(dir,"V21-30_SST-age_MgCa.txt",sep="")),ncol=4,byrow=TRUE))
	names(sEP4) <- c("depth","age","var","SST1")

	EP5 <- c(-2,276)		#Koutavas2006_V19.28 <- c(-2.5,-84.7) = 275.3
	sEP5 <- data.frame(matrix(scan(paste(dir,"V19-28_SST-age_MgCa.txt",sep="")),ncol=4,byrow=TRUE))
	names(sEP5) <- c("depth","age","var","SST1")

	cleanEP = c(2,0,0,0,0)									# 0 = reductive; 1 = oxidative; 2 = other
	depthsEP = c(1364,2830,2921,617,2720)					# water depth (m)
	depthsEP = depthsEP/1000								# water depth (km)
	dco3EP = c(1.23, -6.23, -7.74, 5.74, -3.57)

	coord.EP <- rbind(EP1,EP2,EP3,EP4,EP5)
	EP <- list(sEP1,sEP2,sEP3,sEP4,sEP5)
	
	for (i in 1:length(EP)){
		Dek = EP[[i]][,3]
		SST2 = (log(Dek/0.33)/0.09) - (0.042*dco3EP[i])

		if (cleanEP[i] == 1){
			Dekc = EP[[i]][,3] - (EP[[i]][,3] * 0.15)
			SST3 = (log(Dekc/0.33)/0.09) - (0.042*dco3EP[i])
		} else {
			Dekc = EP[[i]][,3]
			SST3 = (log(Dekc/0.33)/0.09) - (0.042*dco3EP[i])			
		}

		EP[[i]] = cbind(EP[[i]],SST2,SST3)

	}


	for (i in 1:length(EP)){
		if (cleanEP[i] == 1){							# Oxidative
			Dek = EP[[i]][,3] - (EP[[i]][,3] * 0.15)
			Ana = EP[[i]][,3]
		} else if (cleanEP[i] == 0) {					# Reductive
			Dek = EP[[i]][,3]
			Ana = EP[[i]][,3]/0.85
		} else {										# Flow-thru
			Dek = EP[[i]][,3]
			Ana = EP[[i]][,3]
		}

		SST4 = (log(Dek/0.38)/0.09) + 0.61*depthsEP[i] + 1.6 				# Dekens with cleaning (ox)
		SST5 = log(Ana/0.38)/0.09											# Anand with cleaning (red)
		SST6 = (log(EP[[i]][,3]/0.38)/0.09) + 0.61*depthsEP[i] + 1.6 		# Dekens et al. w.o "oxidative" corrections
		SST7 = log(EP[[i]][,3]/0.38)/0.09									# Anand et al. w.o "reductive" corrections

		EP[[i]] = cbind(EP[[i]],SST4,SST5,SST6,SST7)
	}

	# depthsEP = c(1364,2830,2921,617,2720)					# water depth (m)
	# depthsEP = depthsEP/1000								# water depth (km)

	# for (i in 1:length(EP)){
	# 	SST2 = (log(EP[[i]][,3]/0.38)/0.09) + 0.61*depthsEP[i] 		# Dekens et al. w.o "oxidative" corrections
	# 	SST3 = log(EP[[i]][,3]/0.38)/0.09									# Anand et al. w.o "reductive" corrections
	# 	EP[[i]] = cbind(EP[[i]], SST2, SST3)
	# }

	names(EP) <- c("1","2","3","4","5")


	###################################################################################################
	## Read in the West Pacific SST Reconstructions	
	dir <- "/Users/emilygill/Documents/University of Colorado/PHD Research/6. SST Reconstruction/cores/MgCa/WP/"	

	WP1 <- c(8,122)			#Rosenthal_MD97.2141 <- c(8.8,121)
	sWP1 <- data.frame(matrix(scan(paste(dir,"MD97-2141_SST-age_MgCa.txt",sep="")),ncol=4,byrow=TRUE))
	names(sWP1) <- c("depth","age","var","SST1")

	WP2 <- c(6,114)			#Steinke_MD01.2390 <- c(6.6,113)
	sWP2 <- data.frame(matrix(scan(paste(dir,"MD01-2390_SST-age_MgCa.txt",sep="")),ncol=4,byrow=TRUE))
	names(sWP2) <- c("depth","age","var","SST1")

	WP3 <- c(6,126)			#Stott_MD98.2181 <- c(6.3,126)
	sWP3 <- data.frame(matrix(scan(paste(dir,"MD98-2181_SST-age_MgCa.txt",sep="")),ncol=4,byrow=TRUE))
	names(sWP3) <- c("depth","age","var","SST1")

	WP4 <- c(6,128)					#Bolliet_MD06-3067 <- c(6.514, 126.5) ** Had to move this slightly east because of Stott record at 6,126
	sWP4 <- data.frame(matrix(scan(paste(dir,"MD06-3067_bolliet2011_MgCa.txt",sep="")),ncol=4,byrow=TRUE))
	names(sWP4) <- c("depth","age","var","SST1")

	WP5 <- c(2,146)			#deGaridelThron2007_MD97-2138 <- c(1.25,146.1)
	sWP5 <- data.frame(matrix(scan(paste(dir,"MD97-2138_degaridel-thoron2007_MgCa.txt",sep="")),ncol=4,byrow=TRUE))
	names(sWP5) <- c("depth","age","var","SST1")


	# *** This only has four points so if you don't need it, don't use it.
	# WP6 <- c(6,172)			#Dyaz2013_ODP-871 <- c(5.6,172)
	# sWP6 <- data.frame(matrix(scan(paste(dir,"ODP-871_dyaz2013_MgCa.txt",sep="")),ncol=2,byrow=TRUE))
	# names(sWP6) <- c("depth","age","var","SST1")

	WP6 <- c(0,160)			#Lea2000_ODP-806b <- c(0.3,159.3)
	sWP6 <- data.frame(matrix(scan(paste(dir,"ODP-806b_lea2000_MgCa.txt",sep="")),ncol=4,byrow=TRUE))
	names(sWP6) <- c("depth","age","var","SST1")

	WP7 <- c(-2,100)			#Mohtadi2010_GeoB10029-4 <- c(-1.5,100.1)
	sWP7 <- data.frame(matrix(scan(paste(dir,"GeoB10029-4_MgCa.txt",sep="")),ncol=4,byrow=TRUE))
	names(sWP7) <- c("depth","age","var","SST1")	

	WP8 <- c(-4,118)		#Visser_MD98.2162 <- c(-4.6,118)
	sWP8 <- data.frame(matrix(scan(paste(dir,"MD98-2162_SST-age_MgCa.txt",sep="")),ncol=4,byrow=TRUE))
	names(sWP8) <- c("depth","age","var","SST1")

	WP9 <- c(-4,120)			#Linsley2010_70GGC <- c(-3.56,119.4)
	sWP9 <- data.frame(matrix(scan(paste(dir,"70GGC_MgCa.txt",sep="")),ncol=4,byrow=TRUE))
	names(sWP9) <- c("depth","age","var","SST1")

	WP10 <- c(-4,132)		#Stott_MD98.2176 <- c(-5,133)
	sWP10 <- data.frame(matrix(scan(paste(dir,"MD98-2176_SST-age_MgCa.txt",sep="")),ncol=4,byrow=TRUE))
	names(sWP10) <- c("depth","age","var","SST1")

	WP11 <- c(-6,102)			#Mohtadi2010_GeoB10038-4 <- c(-5.9,103.3)
	sWP11 <- data.frame(matrix(scan(paste(dir,"GeoB10038-4_MgCa.txt",sep="")),ncol=4,byrow=TRUE))
	names(sWP11) <- c("depth","age","var","SST1")

	WP12 <- c(-8,116)			#Linsley2010_13GGC <- c(-7.4,115.2)
	sWP12 <- data.frame(matrix(scan(paste(dir,"13GGC_MgCa.txt",sep="")),ncol=4,byrow=TRUE))
	names(sWP12) <- c("depth","age","var","SST1")	

	WP13 <- c(-10,118)			#Levi_MD98-2165 <- c(-9.6,118.3)
	sWP13 <- data.frame(matrix(scan(paste(dir,"MD98-2165_MgCa.txt",sep="")),ncol=4,byrow=TRUE))
	names(sWP13) <- c("depth","age","var","SST1")

	WP14 <- c(-10,126)			#Stott_MD98-2170 <- c(-10.6,125.3)
	sWP14 <- data.frame(matrix(scan(paste(dir,"MD98-2170_MgCa.txt",sep="")),ncol=4,byrow=TRUE))
	names(sWP14) <- c("depth","age","var","SST1")

	cleanWP = c(1,1,0,0,1,0,1,1,0,0,1,0,1,0)				# 0 = reductive; 1 = oxidative; 2 = other
	depthsWP = c(3633,1545,2114,1575,1960,2520,
		964,1855,482,2382,1819,594,2100,832)				# water depth (m)
	depthsWP = depthsWP/1000								# water depth (km)
	dco3WP =  c(0.73, 12.33, 2.26, 7.16, 11.48, 3.47, 19.33, 5.52, 23.31, 4.22, 11.72, 23.30, 12.98, 16.75)

	coord.WP <- rbind(WP1,WP2,WP3,WP4,WP5,WP6,WP7,WP8,WP9,WP10,WP11,WP12,WP13,WP14)
	WP <- list(sWP1,sWP2,sWP3,sWP4,sWP5,sWP6,sWP7,sWP8,sWP9,sWP10,sWP11,sWP12,sWP13,sWP14)
	

	for (i in 1:length(WP)){
		Dek = WP[[i]][,3]
		SST2 = (log(Dek/0.33)/0.09) - (0.042*dco3WP[i])

		if (cleanWP[i] == 1){
			Dekc = WP[[i]][,3] - (WP[[i]][,3] * 0.15)
			SST3 = (log(Dekc/0.33)/0.09) - (0.042*dco3WP[i])
		} else {
			Dekc = WP[[i]][,3]
			SST3 = (log(Dekc/0.33)/0.09) - (0.042*dco3WP[i])			
		}

		WP[[i]] = cbind(WP[[i]],SST2,SST3)

	}
	
	for (i in 1:length(WP)){

		if (cleanWP[i] == 1){							# Oxidative
			Dek = WP[[i]][,3] - (WP[[i]][,3] * 0.15)
			Ana = WP[[i]][,3]
		} else if (cleanWP[i] == 0) {					# Reductive
			Dek = WP[[i]][,3]
			Ana = WP[[i]][,3]/0.85
		} else {										# Flow-thru
			Dek = WP[[i]][,3]
			Ana = WP[[i]][,3]
		}

		SST4 = (log(Dek/0.38)/0.09) + 0.61*depthsWP[i] + 1.6 				# Dekens et al. with cleaning corrections
		SST5 = log(Ana/0.38)/0.09											# Anand et al. with cleaning corrections
		SST6 = (log(WP[[i]][,3]/0.38)/0.09) + 0.61*depthsWP[i] + 1.6 		# Dekens et al. w.o "oxidative" corrections
		SST7 = log(WP[[i]][,3]/0.38)/0.09									# Anand et al. w.o "reductive" corrections

		WP[[i]] = cbind(WP[[i]],SST4,SST5,SST6,SST7)
	}

		# depthsWP = c(3633,1545,2114,1575,1960,2520,
		# 	964,1855,482,2382,1819,594,2100,832)				# water depth (m)
		# depthsWP = depthsWP/1000								# water depth (km)

		# for (i in 1:length(WP)){
		# 	SST2 = (log(WP[[i]][,3]/0.38)/0.09) + 0.61*depthsWP[i] 		# Dekens et al. w.o "oxidative" corrections
		# 	SST3 = log(WP[[i]][,3]/0.38)/0.09									# Anand et al. w.o "reductive" corrections
		# 	WP[[i]] = cbind(WP[[i]], SST2, SST3)
		# }

	names(WP) <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14")	


	###################################################################################################
	## Read in the Indian Ocean SST Reconstructions
	dir <- "/Users/emilygill/Documents/University of Colorado/PHD Research/6. SST Reconstruction/cores/MgCa/IO/"

	IO1 <- c(14,72)				# Govil and Naidu 2010 AAS9/21 <- c(14.5,72)
	sIO1 <- data.frame(matrix(scan(paste(dir,"AAS9.21_MgCa.txt",sep="")),ncol=4,byrow=TRUE))
	names(sIO1) <- c("depth","age","var","SST1")								

	IO2 <- c(10,74)			#Saraswat2014_SK237-GC04 <- c(10.9,74.9)
	sIO2 <- data.frame(matrix(scan(paste(dir,"SK237-GC04_saraswat2013_MgCa.txt",sep="")),ncol=4,byrow=TRUE))
	names(sIO2) <- c("depth","age","var","SST1")

	IO3 <- c(2,78)				#Saraswat2005_SK157-GC04 <- c(-2.4,78)
	sIO3 <- data.frame(matrix(scan(paste(dir,"SK157-GC04_MgCa.txt",sep="")),ncol=4,byrow=TRUE))
	names(sIO3) <- c("depth","age","var","SST1")						

	# IO4 <- c(-10,50)			# Kiefer2006_WIND 28K  <- c(-10.1,51)
	# sIO4 <- data.frame(matrix(scan(paste(dir,"WIND28K_MgCa.txt",sep="")),ncol=4,byrow=TRUE))
	# names(sIO4) <- c("depth","age","var","SST1")								



	coord.IO <- rbind(IO1,IO2,IO3)
	IO <- list(sIO1,sIO2,sIO3)
	for (i in 1:length(IO)){
		SST2 = (log(IO[[i]][,3]/0.38))/0.09  			# Anand et al. 2003/Dekens et al. 2002
		SST3 = ((log(IO[[i]][,3]/0.38))/0.09) 			# Anand et al. 2003/Dekens et al. 2002
		if (i == 2 | i == 3){
			SST2 = IO[[i]][,4]
			SST3 = IO[[i]][,4]
		}
		IO[[i]] = cbind(IO[[i]],SST2,SST3)
	}
	names(IO) <- c("1","2","3")	


	###################################################################################################
	## Return data
	results = list(EPlist=EP, WPlist=WP, IOlist=IO, 
		EPcoord=coord.EP, WPcoord=coord.WP, IOcoord=coord.IO)
	return(results)



	###################################################################################################
	## Optional plot to show all the reconstruction points
	if (plot == TRUE){
		par(mfrow = c(1, 1), mar = c(0.75, 0.75, 0.75, 0.75), oma = c(4, 4, 3, 1))		
		allcoord = rbind(coord.IO, coord.WP, coord.EP)

		xgrid <- seq(30,300,by=2)
		nx <- length(xgrid) # nrows = 101
		ygrid <- seq(-10,28,by=2)
		ny <- length(ygrid) # ncols = 11

		zfull = rep(NaN,(nx*ny))
		zmat = matrix(zfull,nrow=nx,ncol=ny)
		image(xgrid,ygrid,zmat,
				ylim=range(min(ygrid),max(ygrid)), xlim=range(min(xgrid),max(xgrid)),
				xlab="",ylab="", cex.axis=1.5,
				)
		mapnames <- map("world2", xlim=c(min(ygrid), max(xgrid)), ylim=c(min(ygrid), max(ygrid)), 
		    boundary=TRUE, interior=TRUE, col="black", add=TRUE
			)
		points(allcoord[,2],allcoord[,1],col="red")

	}

}