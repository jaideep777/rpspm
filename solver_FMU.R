rm(list=ls())

require(deSolve)

#### Load model ####

# source("codes/pspm_package/Rscripts/wave_model.R")
source("codes/pspm_package/Rscripts/test_model.R")


#### FMU core ####
phi <- function(rr){
  sapply(X = rr, FUN = function(r){max(max(0.0,min(2*r,1.0)),min(r,2.0))})
}

calc_dU_FMU <- function(U, X, t, E){
  
  birthArray = birthRate(X,t,E)
  growthArray = growthRate(x,t,E)
  mortalityArray = mortalityRate(X,t,E)
  
  birthFlux = sum(birthArray*U*h)
  
  u = numeric(J+1)

  u[1] = birthFlux/(growthArray[1]+1e-12)
  u[2] = 2*U[1]-u[1] # Calc with trapezoildal rule. But for g(x)<0, This could be calc with upwind scheme
  
  for (i in 3:(J-1)){  
    if(growthArray[i] >=0){ # needs U @ i-2, i-1, i
      rMinus = (((U[i]-U[i-1])/(x[i]-x[i-1]))/((U[i-1]-U[i-2]+1e-12)/(x[i-1]-x[i-2]))) 
      u[i] = U[i-1] + phi(rMinus)*(U[i-1]-U[i-2])*(x[i]-x[i-1])/(x[i+1]-x[i-1])
    }
    else{ # needs U @ i-1, i, i+1
      rPlus  = ((U[i]-U[i-1])/(x[i]-x[i-1]))/((U[i+1]-U[i]+1e-12)/(x[i+1]-x[i]));
      u[i] = U[i] - phi(rPlus)*(U[i+1]-U[i])*(x[i+1]-x[i])/(x[i+2]-x[i]);
    }
    # u[i] = 2*U[i-1] - u[i-1]
  }
  
  u[J] = 2*U[J-1] - u[J-1] # For g(x)>0, This could be calc with upwind scheme
  u[J+1] = 2*U[J] - u[J] # Calc with trapezoidal rule 

  dU = numeric(J)
  for (i in 1:J){
    # dU[i] = -mortalityArray[i]*U[i] - (growthArray[i+1]*u[i+1] - growthArray[i]*u[i])/h[i]
    dU[i] = -mortalityArray[i]*U[i] - (growthArray[i+1]*u[i+1] - growthArray[i]*u[i])/h[i]
  }
  # cat("du=", dU, "\n")
  
  dU
}

func <- function(t, U, par){
  E = calcEnv(x, U)
  dU = calc_dU_FMU(U,X,t,E)  
  
  list(dU=dU)
}

#### Solve PSPM ####

J = 25                          # Number of grid cells
xb = 0                          # Minimum size (size at birth)
xm = 1                          # Maximum size (size at maturity)
x = seq(xb, xm, length.out=J+1) # grid edges
X = x[-1]-diff(x)/2             # grid centers
h = diff(x)                     # grid widths

# if logarithmiz grid desired:
#x = c(xb, exp(seq(log(xb+0.01), log(xm), length.out=J))) # seq(0,1, length.out=J+1)

U = calcIC(X)

Y = lsoda(U, seq(0,7.5,length.out=20), func, par=0, rtol=1e-4, atol=1e-4)

matplot(x=X, y=t(Y[,-1]), type = "l", lty=1, col=colorRampPalette(colors = c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)), alpha=T)(20), ylim=c(0,2), xlim=c(0,1), xlab="Size", ylab="Density", main="FMU")
points(x=X, y=Y[1,-1], col="grey", type="l", lwd=2)
points(x=x, y=Ueq(x), col="green3", type="l", lwd=3)
points(x=X, y=Y[nrow(Y),-1], col="blue", pch=20, cex=1)
abline(v=X, col="grey")





