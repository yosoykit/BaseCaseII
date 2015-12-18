genHexGrid <- function(ll = c(0, 0), ur = c(1, 1), totalcrypts)
{   
    ### Search for optimal dx s.t. dx*dy is within epsilon close of totalcrypts ###
    crypt_num <- function(x_temp){
        dy <- sqrt(3) * x_temp / 2

        x <- seq(ll[1], ur[1] - x_temp / 2,x_temp)
        y <- seq(ll[2], ur[2], dy)

        if(length(y) %% 2 != 0) y <- y[-length(y)] 
        nx <- length(x)
        ny <- length(y)
        return(abs(nx*ny-totalcrypts))
    }
    dx <- optim(.14,crypt_num,method="Brent",lower=.12,upper=.15)
    dy <- sqrt(3) * dx$par / 2

    x <- seq(ll[1], ur[1] - dx$par / 2, dx$par)
    y <- seq(ll[2], ur[2], dy)

    if(length(y) %% 2 != 0) y <- y[-length(y)] 
    
    nx <- length(x)
    ny <- length(y)
    y <- rep(y, each = nx)
    x <- rep(c(x, x + dx$par / 2), length.out = length(y))

    x <- x + (ur[1] - max(x)) / 2
    y <- y + (ur[2] - max(y)) / 2
    list(x = x, y = y, nx = nx, ny = ny)
}

genPolyList <- function(hexGrid)
{
        dx <- hexGrid$x[2] - hexGrid$x[1]
        dy <- dx / sqrt(3)

        x.offset <- c(-dx / 2, 0, dx / 2, dx / 2, 0, -dx / 2)
        y.offset <- c(dy / 2, dy, dy / 2, -dy / 2, -dy, -dy / 2)

        f <- function(i) list(x = hexGrid$x[i] + x.offset,
                              y = hexGrid$y[i] + y.offset)

        lapply(1:length(hexGrid$x), f)
}

#xy <- genHexGrid(0.1)
#plot(xy, asp = 1, axes = F, pch = 19, xlab = NA, ylab = NA, cex = 0.5, col = 'mistyrose')
# lapply(genPolyList(xy), polygon, xpd = NA, ypd = NA)

