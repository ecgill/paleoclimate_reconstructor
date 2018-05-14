read.recon.uk3 <- function(plot){
	
	###################################################################################################
	## Read in the East Pacific SST Reconstructions
	dir <- "/Users/emilygill/Documents/University of Colorado/PHD Research/6. SST Reconstruction/cores/Uk37/EP/"

	EP1 <- c(8,276)			#Leduc_MD02.2529 <- c(8.2,-84.1) = 275.9
	sEP1 <- data.frame(matrix(scan(paste(dir,"MD02-2529_SST-age_alk.txt",sep="")),ncol=4,byrow=TRUE))
	names(sEP1) <- c("depth","age","var","SST1")

	EP2 <- c(8,278)			#Dubois_ME0005A-43JC <- c(7.85,-83.6) = 276.4
	sEP2 <- data.frame(matrix(scan(paste(dir,"ME0005A-43JC_dubois.txt",sep="")),ncol=4,byrow=TRUE))
	names(sEP2) <- c("depth","age","var","SST1")	

	EP3 <- c(4,282)			#Pahnke_KNR176_JPC32 <- c(4.8,-78) = 282
	sEP3 <- data.frame(matrix(scan(paste(dir,"KNR176-JPC32_SST-age_alk.txt",sep="")),ncol=4,byrow=TRUE))
	names(sEP3) <- c("depth","age","var","SST1")

	EP4 <- c(2,270)			#Dubois_TR163-19 <- c(2.25, -90.95) = 269.05
	sEP4 <- data.frame(matrix(scan(paste(dir,"TR163-19P_dubois.txt",sep="")),ncol=4,byrow=TRUE))
	names(sEP4) <- c("depth","age","var","SST1") 

	EP5 <- c(2,274) 		#Kienst_ME0005A.24JC <- c(1.5,-86) = 274
	sEP5 <- data.frame(matrix(scan(paste(dir,"ME0005A-24JC_SST-age_alk.txt",sep="")),ncol=4,byrow=TRUE))
	names(sEP5) <- c("depth","age","var","SST1")

	EP6 <- c(0,268) 		#Dubois_TR162-22 <- c(0.52,-92.40) = 267
	sEP6 <- data.frame(matrix(scan(paste(dir,"TR163-22P_dubois.txt",sep="")),ncol=4,byrow=TRUE))
	names(sEP6) <- c("depth","age","var","SST1")	

	EP7 <- c(0,278)			#KoutavasSachs_V19.27 <- c(-0.5,-82.67) = 277.33
	sEP7 <- data.frame(matrix(scan(paste(dir,"V19-27_SST-age_alk.txt",sep="")),ncol=4,byrow=TRUE))
	names(sEP7) <- c("depth","age","var","SST1")

	EP8 <- c(-2,270) 		#KoutavasSachs_V21.30 <- c(-1.2,-90) = 270
	sEP8 <- data.frame(matrix(scan(paste(dir,"V21-30_SST-age_alk.txt",sep="")),ncol=4,byrow=TRUE))
	names(sEP8) <- c("depth","age","var","SST1")

	EP9 <- c(-2,274)		#KoutavasSachs_RC11.238 <- c(-1.5,-86) = 274
	sEP9 <- data.frame(matrix(scan(paste(dir,"RC11-238_SST-age_alk.txt",sep="")),ncol=4,byrow=TRUE))
	names(sEP9) <- c("depth","age","var","SST1")

	EP10 <- c(-2,276)		#KoutavasSachs_V19.28 <- c(-2.5,-85) = 275.35
	sEP10 <- data.frame(matrix(scan(paste(dir,"V19-28_SST-age_alk.txt",sep="")),ncol=4,byrow=TRUE))
	names(sEP10) <- c("depth","age","var","SST1")

	EP11 <- c(-2,278)		#Dubois_ME0005A-27JC <- c(-1.85,-82.78) = 277.2
	sEP11 <- data.frame(matrix(scan(paste(dir,"ME0005A-27JC_dubois.txt",sep="")),ncol=4,byrow=TRUE))
	names(sEP11) <- c("depth","age","var","SST1")		

	# EP12 <- c(-4,276)		#Dubois_TR163-31 <- c(-3.62, -83.97) = 276.03
	# sEP12 <- data.frame(matrix(scan(paste(dir,"TR163-31P_dubois.txt",sep="")),ncol=4,byrow=TRUE))
	# names(sEP12) <- c("depth","age","var","SST1")

	EP12 <- c(-4,276)		#KoutavasSachs_V19.30 <- c(-3.4,-84) = 276
	sEP12 <- data.frame(matrix(scan(paste(dir,"V19-30_SST-age_alk.txt",sep="")),ncol=4,byrow=TRUE))
	names(sEP12) <- c("depth","age","var","SST1")

	EP13 <- c(-4,278)		#Bova_CDH26 <- c(-3.59,-81.18) = 278.82
	sEP13 <- data.frame(matrix(scan(paste(dir,"CDH26_bova.txt",sep="")),ncol=4,byrow=TRUE))
	names(sEP13) <- c("depth","age","var","SST1")

	# EP9 <- c(0,274)		#Prahl_Y69-71P <- c(0,-86.4) = 273.6
	# sEP9 <- data.frame(matrix(scan(paste(dir,"Y69-71P_SST-age_alk.txt",sep="")),ncol=4,byrow=TRUE))
	# names(sEP9) <- c("depth","age","var","SST1")

	coord.EP <- rbind(EP1,EP2,EP3,EP4,EP5,EP6,EP7,EP8,EP9,EP10,EP11,EP12,EP13)
	EP <- list(sEP1,sEP2,sEP3,sEP4,sEP5,sEP6,sEP7,sEP8,sEP9,sEP10,sEP11,sEP12,sEP13)

	for (i in 1:length(EP)){
		SST2 = (EP[[i]][,3] - 0.044) / 0.033  			# Muller et al. 1998
		SST3 = (EP[[i]][,3] - 0.039) / 0.034 			# Prahl et al. 1988
		EP[[i]] = cbind(EP[[i]],SST2,SST3)
	}

	names(EP) <- c("1","2","3","4","5","6","7","8","9","10","11","12","13")

	# coord.EP <- rbind(EP2,EP3,EP4,EP5,EP6,EP7,EP8)
	# EP <- list(sEP2,sEP3,sEP4,sEP5,sEP6,sEP7,sEP8)
	# names(EP) <- c("1","2","3","4","5","6","7")

	# coord.EP <- rbind(EP1,EP3,EP4,EP5,EP6,EP7,EP9)
	# EP <- list(sEP1,sEP3,sEP4,sEP5,sEP6,sEP7,sEP9)
	# names(EP) <- c("1","2","3","4","5","6","7")


	###################################################################################################
	## Read in the West Pacific SST Reconstructions	
	dir <- "/Users/emilygill/Documents/University of Colorado/PHD Research/6. SST Reconstruction/cores/Uk37/WP/"

	WP1 <- c(10,110)		#Kienast_GIK18252-3 c(9.2,109.4)
	sWP1 <- data.frame(matrix(scan(paste(dir,"GIK18252-3_SST-age_alk.txt",sep="")),ncol=4,byrow=TRUE))
	names(sWP1) <- c("depth","age","var","SST1")

	WP2 <- c(8,110)		#Zhao_MD97-2151 c(8.7,109.7)
	sWP2 <- data.frame(matrix(scan(paste(dir,"MD97-2151_SST-age_alk.txt",sep="")),ncol=4,byrow=TRUE))
	names(sWP2) <- c("depth","age","var","SST1")

	WP3 <- c(6,110)		#Kienast_GIK18287-3 c(5.7,110.7)
	sWP3 <- data.frame(matrix(scan(paste(dir,"GIK18287-3_SST-age_alk.txt",sep="")),ncol=4,byrow=TRUE))
	names(sWP3) <- c("depth","age","var","SST1")
	
	WP4 <- c(6,112)		#Pelejero_GIK17964-1 c(6.15,112.2)
	sWP4 <- data.frame(matrix(scan(paste(dir,"GIK17964-1_SST-age_alk.txt",sep="")),ncol=4,byrow=TRUE))
	names(sWP4) <- c("depth","age","var","SST1")

	WP5 <- c(6,126)		#Fraser_MD06-3075 c(6.48, 125.83) *** had to move to 8N to not coincide with MgCa Stott2181
	sWP5 <- data.frame(matrix(scan(paste(dir,"MD06-3075_fraser2014_alk.txt",sep="")),ncol=4,byrow=TRUE))
	names(sWP5) <- c("depth","age","var","SST1")

	WP6 <- c(2,146)			#deGaridelThron2007_MD97-2138 <- c(1.25,146.1)
	sWP6 <- data.frame(matrix(scan(paste(dir,"MD97-2138_degaridel-thoron2007_alk.txt",sep="")),ncol=4,byrow=TRUE))
	names(sWP6) <- c("depth","age","var","SST1")

	WP7 <- c(-6,104)		#Luckge_SO139-74KL c(-6.5,103.8)
	sWP7 <- data.frame(matrix(scan(paste(dir,"SO139-74KL_SST-age_alk.txt",sep="")),ncol=4,byrow=TRUE))
	names(sWP7) <- c("depth","age","var","SST1")

	coord.WP <- rbind(WP1,WP2,WP3,WP4,WP5,WP6,WP7)
	WP <- list(sWP1,sWP2,sWP3,sWP4,sWP5,sWP6,sWP7)

	for (i in 1:length(WP)){
		SST2 = (WP[[i]][,3] - 0.044) / 0.033  			# Muller et al. 1998
		SST3 = (WP[[i]][,3] - 0.039) / 0.034 			# Prahl et al. 1988
		WP[[i]] = cbind(WP[[i]], SST2, SST3)
	}

	names(WP) <- c("1","2","3","4","5","6","7")	




	###################################################################################################
	## Read in the Indian Ocean SST Reconstructions
	dir <- "/Users/emilygill/Documents/University of Colorado/PHD Research/6. SST Reconstruction/cores/Uk37/IO/"
	
	# IO1 <- c(26,36)			#Arz2006_GeoB-5836-2 <- c(26.21,35.35)
	# sIO1 <- data.frame(matrix(scan(paste(dir,"GeoB-5836-2_arz2006_alk.txt",sep="")),ncol=4,byrow=TRUE))
	# names(sIO1) <- c("depth","age","var","SST1")

	IO1 <- c(24,64)			#Schulz2002_SO90-93KL <- c(23.6,64.2)
	sIO1 <- data.frame(matrix(scan(paste(dir,"SO90-93KL_SST-age_alk.txt",sep="")),ncol=4,byrow=TRUE))
	names(sIO1) <- c("depth","age","var","SST1")

	IO2 <- c(24,66)			#vonRad1999&Doose-Rolnski2001_SO90-39KG <- c(24.8,65.9)
	sIO2 <- data.frame(matrix(scan(paste(dir,"SO90-39KG_SST-age_alk.txt",sep="")),ncol=4,byrow=TRUE))
	names(sIO2) <- c("depth","age","var","SST1")

	IO3 <- c(22,66)			#SchulteMueller2001_SO90-136KL <- c(23.1,66.5) *ok to make it 22?
	sIO3 <- data.frame(matrix(scan(paste(dir,"SO90-136KL_SST-age_alk.txt",sep="")),ncol=4,byrow=TRUE))
	names(sIO3) <- c("depth","age","var","SST1")

	IO4 <- c(20,90)			#Kudrass2001_SO93-126KL <- c(19.9,90)
	sIO4 <- data.frame(matrix(scan(paste(dir,"SO93-126KL_SST-age_alk.txt",sep="")),ncol=4,byrow=TRUE))
	names(sIO4) <- c("depth","age","var","SST1")

	# IO5 <- c(16,60)			#Herbert2010_ODP-722 <- c(16.6,59.8)
	# sIO5 <- data.frame(matrix(scan(paste(dir,"ODP-722_herbert2010_alk.txt",sep="")),ncol=4,byrow=TRUE))
	# names(sIO5) <- c("depth","age","var","SST1")

	IO5 <- c(14,54)		#Rostek1997_TY93929/P <- c(13.7,53.3)
	sIO5 <- data.frame(matrix(scan(paste(dir,"TY93929_P_SST-age_alk.txt",sep="")),ncol=4,byrow=TRUE))
	names(sIO5) <-  c("depth","age","var","SST1")

	IO6 <- c(12,52)		#Kim2004_TY93-905 <- c(11.1,51.9)
	sIO6 <- data.frame(matrix(scan(paste(dir,"TY93-905_SST-age_alk.txt",sep="")),ncol=4,byrow=TRUE))
	names(sIO6) <- c("depth","age","var","SST1")

	IO7 <- c(10,76)			#Sonzogni1998_MD77-194 <- c(10.3,75.1)
	sIO7 <- data.frame(matrix(scan(paste(dir,"MD77-194_SST-age_alk.txt",sep="")),ncol=4,byrow=TRUE))
	names(sIO7) <- c("depth","age","var","SST1")

	IO8 <- c(0,46)			#Bard1997_MD85-668 <- c(0,46)
	sIO8 <- data.frame(matrix(scan(paste(dir,"MD85-668_bard1997_alk.txt",sep="")),ncol=4,byrow=TRUE))
	names(sIO8) <- c("depth","age","var","SST1")

	IO9 <- c(-6,104)		#Luckge2009_ SO139-74KL<- c(-6.5,103.8)
	sIO9 <- data.frame(matrix(scan(paste(dir,"SO139-74KL_SST-age_alk.txt",sep="")),ncol=4,byrow=TRUE))
	names(sIO9) <- c("depth","age","var","SST1")

	coord.IO <- rbind(IO1,IO2,IO3,IO4,IO5,IO6,IO7,IO8,IO9)
	IO <- list(sIO1,sIO2,sIO3,sIO4,sIO5,sIO6,sIO7,sIO8,sIO9)

	for (i in 1:length(IO)){
		SST2 = (IO[[i]][,3] - 0.044) / 0.033  			# Muller et al. 1998
		SST3 = (IO[[i]][,3] - 0.039) / 0.034 			# Prahl et al. 1988
		IO[[i]] = cbind(IO[[i]],SST2,SST3)
	}

	names(IO) <- c("1","2","3","4","5","6","7","8","9")

	

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