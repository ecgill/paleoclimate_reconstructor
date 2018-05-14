BINxygrid <- function(data, xgrid, ygrid, ntime, type, undefined, condition) {
    ny <- length(ygrid)
    nx <- length(xgrid)
    ngrid <- nx*ny
    xygrid <- matrix(0, nrow=(nx*ny), ncol=2)
    
    # Create the xygrid of lat and lon for each gridcell.
    i=0
    for (iy in 1:ny){
        for (ix in 1:nx){
            i=i+1
            xygrid[i,1] = ygrid[iy]
            xygrid[i,2] = xgrid[ix] 
        }
    }

    # Option for data that is defined in every single gridcell.
    if (type == "ngrid"){
        nsta = nx * ny
        data <- array(data=data, dim=c(nx,ny,ntime))
        BINdata <- matrix(NA, nrow=ntime, ncol=ngrid)
    
        for (i in 1:ntime){
            data1 <- data[,,i]
            BINdata[i,] <- data1
        }
        rm("data")

        # Return the results.
        results <- list(data = BINdata, grid = xygrid, nsta=nsta, xgrid=xgrid, ygrid=ygrid)
        return(results)
    }
    

    # Option for data that undefined in some gridcells (i.e. rainfall data just over India).
    if (type == "special"){
        data <- array(data=data, dim=c(nx,ny,ntime))
        
        index <- 1:(nx*ny)
        data1 <- data[,,1]


        if (condition == "greater"){
            index1 <- index[data1 > undefined]    
        }
        
        if (condition == "less"){
            index1 <- index[data1 < undefined]   
        }

        nsta <- length(index1)
        BINdata <- matrix(NA, nrow=ntime, ncol=nsta)
        for (i in 1:ntime){
            data1 <- data[,,i]
            data2 <- data1[index1]
            BINdata[i,] <- data2
        }
    
        # Return the results.
        results <- list(data = BINdata, grid = xygrid, ind=index1, nsta=nsta, xgrid=xgrid, ygrid=ygrid)
        return(results)    
    }


}