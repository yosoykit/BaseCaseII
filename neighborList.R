neighborList = function (xy) {

# neighbor list ... with quasi-periodic bc
N = length(xy$x); nx = xy$nx; ny = xy$ny
nbr = as.list(N)

for (i in (nx+1):(N-nx)) {
  r = ceiling(i/nx) #current row
  if(r%%2==0) {  # if row is even
    nbr[[i]] = c(i-1, i+1, i+nx, i+nx+1, i-nx, i-nx+1)}
  else {  # if row is odd
    nbr[[i]] = c(i-1, i+1, i+nx, i+nx-1, i-nx, i-nx-1)}
}

# bc in first row
nbr[[1]] = c(nx, 2, 1+nx, 2*nx, N-nx+1, N)
for(i in 2:nx) {nbr[[i]] = c(i-1, i+1, i+nx, i+nx-1, N-nx+i, N-nx+i-1)}

# bc in last row
nbr[[N]] = c(N-1, N-nx+1, 1, nx, N-nx, N-2*nx+1) 
for(i in (N-nx+1):(N-1)) {nbr[[i]] = c(i-1, i+1, i-N+nx, i-N+nx+1, i-nx, i-nx+1)}  

return(nbr = nbr)
}
