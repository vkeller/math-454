
#load gpu dynamic library
dyn.load("avg_wrap.so")

# Create a matrix
N <- 1000L
matrix.in  <- abs( rnorm(N * N) )
matrix.out  <- rep(0, N*N);


# Call C function
Radius <- 3
.Call("avg_wrap", matrix.in,N,N,Radius,matrix.out)
            
#convert  to matrix
matrix.out <- matrix( matrix.out, nrow=N, byrow=TRUE)
