library(RColorBrewer)
library(fields)

# East Pacific
dev.new()
par(mfrow = c(1, 1), mar = c(0.75, 0.75, 0.75, 0.75), oma = c(4, 4, 3, 1))		
xgrid = seq(-100,-70,2)
nx = length(xgrid)
ygrid = seq(-12,12,2)
ny = length(ygrid)	
image(xgrid,ygrid,matrix(rep(NaN,(nx*ny)),nx,ny),
	xlab="",ylab="", cex.axis=1.25)

world(ylim=range(min(ygrid),max(ygrid)),
	xlim=range(min(xgrid),max(xgrid)),
	add=TRUE,
	fill=TRUE,
	border="darkgrey")

abline(v=seq(-99,-71,by=2), col="grey25")
abline(h=seq(-11,11,by=2), col="grey25")
spectral <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")

# For alkenones
cbPalette = rev(spectral(14))
points(-84.12, 8.21, col=cbPalette[1])	#27	
points(-83.6, 7.85, col=cbPalette[2])	#28
points(-77.96, 4.85, col=cbPalette[3])	#29
points(-90.95, 2.25, col=cbPalette[4])	#30
points(-86.485, 1.5, col=cbPalette[5])	#31
points(-92.40, 0.52, col=cbPalette[6])	#32
points(-82.67, -0.47, col=cbPalette[7])	#33
points(-89.68, -1.22, col=cbPalette[8])	#34
points(-85.82, -1.52, col=cbPalette[9])	#35
points(-82.78, -1.85, col=cbPalette[10])	#36
points(-84.65, -2.51, col=cbPalette[11])	#37
points(-83.52, -3.38, col=cbPalette[12])	#38 v19-30 k&s
points(-83.97, -3.62, col=cbPalette[13])	#   tr163-31 dubois - this was omitted bc #38 has higher res
points(-81.18,-3.59, col=cbPalette[14])		#39 bova


# For MgCa
cbPalette = rev(spectral(5))
points(-83.61, 7.86, col=cbPalette[1])	#22
points(-92.40, 0.52, col=cbPalette[2])	#23
points(-86.45, 0.02, col=cbPalette[3])	#24
points(-89.68, -1.22, col=cbPalette[4])	#25
points(-84.65, -2.51, col=cbPalette[5])	#26








# West Pacific
dev.new()
par(mfrow = c(1, 1), mar = c(0.75, 0.75, 0.75, 0.75), oma = c(4, 4, 3, 1))		
xgrid = seq(100,160,2)
nx = length(xgrid)
ygrid = seq(-12,12,2)
ny = length(ygrid)	
image(xgrid,ygrid,matrix(rep(NaN,(nx*ny)),nx,ny),
	xlab="",ylab="", cex.axis=1.25)

world(ylim=range(min(ygrid),max(ygrid)),
	xlim=range(min(xgrid),max(xgrid)),
	add=TRUE,
	fill=TRUE,
	border="darkgrey")

abline(v=seq(99,161,by=2), col="grey25")
abline(h=seq(-11,11,by=2), col="grey25")
spectral <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")

# For alkenones
cbPalette = rev(spectral(7))
points(109.38, 9.23, col=cbPalette[1])	#15
points(109.87, 8.72, col=cbPalette[2])	#16
points(110.65, 5.65, col=cbPalette[3])	#17
points(125.83, 6.48, col=cbPalette[4])	#18
points(112.21, 6.16, col=cbPalette[5])	#19
points(146.14, 1.25, col=cbPalette[6])	#20
points(103.83, -6.54, col=cbPalette[7])	#21


# For MgCa
cbPalette = rev(spectral(14))
points(121.30, 8.80, col=cbPalette[1])	#1
points(113.41, 6.64, col=cbPalette[2])	#2
points(125.83, 6.30, col=cbPalette[3])	#3
points(126.50, 6.51, col=cbPalette[4])	#4
points(146.14, 1.25, col=cbPalette[5])	#5
points(159.35, 0.32, col=cbPalette[6])	#6
points(100.1, -1.5, col=cbPalette[7])	#7
points(117.90, -4.69, col=cbPalette[8])	#8
points(119.4, -3.56, col=cbPalette[9])	#9
points(133.45, -5.0, col=cbPalette[10])	#10
points(103.3, -5.9, col=cbPalette[11])	#11
points(115.2, -7.4, col=cbPalette[12])	#12
points(118.34, -9.65, col=cbPalette[13])	#13
points(125.39, -10.59, col=cbPalette[14])	#14













