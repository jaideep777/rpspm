
createSolver <- function(size, xbirth, xmax, method){
  x = x = seq(xb, xm, length.out=J+1) # grid edges
  X = x[-1]-diff(x)/2                 # grid centers
  list(J=size, xb=xbirth, xm=xmax, method=method, 
       x=x,
       X=X
       )
}

S = createSolver(25, 0, 1, "FMU")

initialize = function(f){
  
}