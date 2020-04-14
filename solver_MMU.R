rm(list=ls())

require(deSolve)
require(limSolve)

#### Load model ####

# source("codes/pspm_package/Rscripts/wave_model.R")
source("codes/pspm_package/Rscripts/test_model.R")


#### FMU core ####
phi <- function(rr){
  sapply(X = rr, FUN = function(r){max(max(0.0,min(2*r,1.0)),min(r,2.0))})
}


calc_dU_MMU <- function(U, x, uprev, t, E){
  X = x[-1]-diff(x)/2             # grid centers
  h = diff(x)
  
  birthArray = birthRate(X,t,E)
  growthArray = growthRate(x,t,E)
  mortalityArray = mortalityRate(X,t,E)
  
  birthFlux = sum(birthArray*U*h)
  
  dx = calc_dx_MMU(x,uprev)

  u = numeric(J+1)

  # JAI: Is is correct to use updated dx for calculation of u?
  u[1] = birthFlux/(growthArray[1]+1e-12)
  u[2] = 2*U[1]-u[1] # Calc with trapezoildal rule. But for g(x)<0, This could be calc with upwind scheme
  
  for (i in 3:(J-1)){  
    if(growthArray[i]-dx[i] >=0){ # needs U @ i-2, i-1, i
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
    dU[i] = -mortalityArray[i]*U[i] - ((growthArray[i+1]-dx[i+1])*u[i+1] - (growthArray[i]-dx[i])*u[i])/h[i] - (dx[i+1]-dx[i])*U[i]/h[i]
  }
  # cat("du=", dU, "\n")

  list(u=u, rates=c(dx,dU))
}


monitor = function(i, x, u, regularizingFactor){
  dx = x[i+1] - x[i] 
  du = (u[i+1] - u[i])/dx
  dxu = (u[i+1]*x[i+1] - u[i]*x[i])/dx
  z1 = sqrt(1 + (1/regularizingFactor)*du^2)
  z2 = sqrt(1 + (1/regularizingFactor)*dxu^2)
  # cat("i/z" = i,z,"\n")
  (z1+z2)/2
  # z1
}


calc_dx_MMU = function(x,u){
  factor = sum(diff(u)^2/diff(x))/(xm-xb)
  a = numeric(J+1)
  b = numeric(J+1)
  c = numeric(J+1)
  d = numeric(J+1)
  b[1] = 1
  b[J+1] = 1 
  
  rho1 = sapply(X = 1:J, FUN = monitor, x=x,u=u,regularizingFactor=factor)
  pm=4
  gamma=2
  rho = rho1
  for(i in (1):(J)){
    p = min(i-1, pm, J-i)
    weights = (gamma/(1+gamma))^(abs(-p:p))
    rho[i] = sqrt(sum(rho1[(i-p):(i+p)]^2*weights)/sum(weights))
  }
  
  for (i in 2:(J)){
    rhoMinus = (rho[i-1]+rho[i])/2
    if (i < J){
      rhoPlus = (rho[i+1]+rho[i])/2
    }else{
      rhoPlus = rho[i]
    }
    rightSide = -(1/tau)*(rhoPlus*(x[i+1]-x[i]) - rhoMinus*(x[i]-x[i-1]));
    a[i] = rhoMinus;
    b[i] = -(rhoMinus + rhoPlus);
    c[i] = rhoPlus;
    d[i] = rightSide;
  }
  # cat(c,"\n")
  # dx = rep(0,J+1) 
  dx = as.numeric(Solve.tridiag(a[-1], b, c[-(J+1)], d))
  dx
}

func <- function(t, y, par){
  x = y[1:(J+1)]
  U = y[(J+2):length(y)]

  E = calcEnv(x, U)

  dxdU = calc_dU_MMU(U,x,uprev,t,E)  
  
  uprev = dxdU$u
  list(dU=dxdU$rates)
}

#### Solve PSPM ####

J = 25                          # Number of grid cells
xb = 0                          # Minimum size (size at birth)
xm = 1                          # Maximum size (size at maturity)
x = seq(xb, xm, length.out=J+1) # grid edges
X = x[-1]-diff(x)/2             # grid centers
# h = diff(x)                     # grid widths
tau = 1

mids = function(x){
  x[-1]-diff(x)/2
}

# if logarithmiz grid desired:
#x = c(xb, exp(seq(log(xb+0.01), log(xm), length.out=J))) # seq(0,1, length.out=J+1)

U = calcIC(X)
uprev = calcIC(x)

Y = lsoda(c(x,U), seq(0,10,length.out=20), func, par=0, rtol=1e-4, atol=1e-4)

plot(x=mids(Y[1,2:(J+2)]), y=Y[1,(J+3):(2*J+2)], col="red", type="l",ylim=c(0,2), xlim=c(0,1), xlab="Size", ylab="Density", main="MMU")
for(i in 2:20){
  points(x=mids(Y[i,2:(J+2)]), y=Y[i,(J+3):(2*J+2)], col=colorRampPalette(c(rgb(1,0,0,0.5),rgb(0,0,1,0.5)),alpha=T)(20)[i], type="l")
}
points(x=seq(0,1,.02), y=Ueq(seq(0,1,.02)), col="green3", type="l", lwd=3)
points(x=mids(Y[1,2:(J+2)]), y=Y[1,(J+3):(2*J+2)], col="red", pch=20)
points(x=mids(Y[20,2:(J+2)]), y=Y[20,(J+3):(2*J+2)], col="blue", pch=20)
abline(v=mids(Y[20,2:(J+2)]), col="grey")

# matplot(x=X, y=t(Y[,-1]), type = "l", lty=1, col=colorRampPalette(colors = c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)), alpha=T)(20))
# points(x=X, y=Y[1,-1], col="grey", type="l", lwd=2)
# points(x=X, y=Y[nrow(Y),-1], col="blue", pch=20, cex=1)





