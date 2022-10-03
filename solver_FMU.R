rm(list=ls())

require(deSolve)

#### Load model ####

# source("codes/pspm_package/Rscripts/wave_model.R")


#### FMU core ####
phi <- function(rr){
  sapply(X = rr, FUN = function(r){max(max(0.0,min(2*r,1.0)),min(r,2.0))})
}

calc_dU_FMU <- function(U, x, t, E){
  X = x[-1]-diff(x)/2             # grid centers
  h = diff(x)

  birthArray = birthRate(X,t,E)
  growthArray = growthRate(x,t,E)
  mortalityArray = mortalityRate(X,t,E)
  
  birthFlux = sum(birthArray*U*h)
  
  u = numeric(J+1)

  u[1] = birthFlux/(growthArray[1]+1e-12)
  # cat("u0 = ", u[1], birthFlux, "\n")
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
  # for (i in 3:(J-1)){
  #   if(growthArray[i] >=0){ # needs U @ i-2, i-1, i
  #     rMinus = (((U[i]-U[i-1])/(X[i]-X[i-1]))/((U[i-1]-U[i-2]+1e-12)/(X[i-1]-X[i-2])))
  #     u[i] = U[i-1] + phi(rMinus)*(U[i-1]-U[i-2])*(X[i]-X[i-1])/(X[i+1]-X[i-1])
  #   }
  #   else{ # needs U @ i-1, i, i+1
  #     rPlus  = ((U[i]-U[i-1])/(X[i]-X[i-1]))/((U[i+1]-U[i]+1e-12)/(X[i+1]-X[i]));
  #     u[i] = U[i] - phi(rPlus)*(U[i+1]-U[i])*(X[i+1]-X[i])/(X[i+2]-X[i]);
  #   }
  #   # u[i] = 2*U[i-1] - u[i-1]
  # }
  
    
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

#### Simple version #####
calc_dU_FMU <- function(U, x, t, E, B){
  X = x[-1]-diff(x)/2             # grid centers
  h = diff(x)
  
  birthArray = birthRate(X,t,E)
  growthArray = growthRate(X,t,E)
  mortalityArray = mortalityRate(X,t,E)
  
  birthFlux = sum(birthArray*U*h)
  dU0 =birthFlux/(growthArray[1]+1e-12)
  
  dU = numeric(J)
  dU[1]=dU0/dt -mortalityArray[1]*U[1]
  for (i in 2:(J)){
    dU[i] = -mortalityArray[i]*U[i] - (growthArray[i]*U[i] - growthArray[i-1]*U[i-1])/h[i-1]
  }
  
  # dU = numeric(J)
  # for (i in 1:(J-1)){
  #   dU[i] = -mortalityArray[i]*U[i] - (growthArray[i+1]*U[i+1] - growthArray[i]*U[i])/h[i]
  #   # dU[i] = -mortalityArray[i]*U[i] - (growthArray[i]*U[i] - growthArray[i-1]*U[i-1])/h[i-1]
  # }
  # dU[1]= dU[1] + dU0/dt
  # dU[J]= -mortalityArray[J]*U[J]
  
  # cat("du=", dU, "\n")
  
  dU
}


#### Lindh18 version #####
calc_dU_FMU <- function(U, x, t, E, B){
  X = x[-1]-diff(x)/2             # grid centers
  h = diff(x)
 
  birthArray = birthRate(X,t,E)
  growthArray = growthRate(X,t,E)
  mortalityArray = mortalityRate(X,t,E)
  
  birthFlux = sum(birthArray*U*h)
  
  Unew = numeric(length(U))
  Unew[1] = U[1]
}

func <- function(t, U, par){
  E = calcEnv(x, U)
  dU = calc_dU_FMU(U,x,t,E)  
  
  list(dU=dU)
}





# #### Solve PSPM ####
# source("test_model.R")
# 
# J = 100                          # Number of grid cells
# xb = 0                        # Minimum size (size at birth)
# xm = 1                         # Maximum size (size at maturity)
# # if logarithmiz grid desired:
# # x = exp(seq(log(xb), log(xm), length.out=J+1)) # seq(0,1, length.out=J+1)
# x = seq(xb, xm, length.out=J+1) # grid edges
# X = x[-1]-diff(x)/2             # grid centers
# h = diff(x)                     # grid widths
# 
# 
# U = calcIC(X)
# 
# # calcEnv(x,U)
# #
#  Y = lsoda(U, seq(0,8,length.out=20), func, par=0, rtol=1e-5, atol=1e-5)
# 
# matplot(x=X, y=t(Y[,-1]), type = "l", lty=1, col=colorRampPalette(colors = c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)), alpha=T)(20), xlab="Size", ylab="Density", main="FMU")
# points(x=X, y=Y[1,-1], col="grey", type="l", lwd=2)
# points(x=x, y=Ueq(x), col="green3", type="l", lwd=3)
# points(x=X, y=Y[nrow(Y),-1], col="blue", pch=20, cex=1)
# abline(v=X, col="grey")
# 
# bo = numeric(nrow(Y))
# boo = bo
# for (i in 1:length(bo)){
#   U = Y[i,-1]
#   bo[i] = sum(birthRate(X,Y[i,1], calcEnv(x, U))*U*h)
#   Uo = (1-X)^2/(1+X)^4 + (1-X)/(1+X)^3*exp(-0.225*Y[i,1]^2)
#   boo[i] = sum(birthRate(X,Y[i,1], calcEnv(x, Uo))*Uo*h)
# }
# 
# plot(boo~Y[,1], type="l")
# points(bo~Y[,1], col="red")
# 


#### Solve PSPM ####
source("codes/rpspm/RED_model_v2.R")

J = 50                          # Number of grid cells
xb = 1                        # Minimum size (size at birth)
xm = 1e6                         # Maximum size (size at maturity)
# if logarithmiz grid desired:
x = exp(seq(log(xb), log(xm), length.out=J)) # seq(0,1, length.out=J+1)
# x = seq(xb, xm, length.out=J+1) # grid edges
X = x[-1]-diff(x)/2             # grid centers
h = diff(x)                     # grid widths

U = calcIC(X)
Y = lsoda(U, seq(0,500,length.out=20), func, par=0, rtol=1e-5, atol=1e-5)

U = calcIC(X)
Y = NULL
y = c(U)
dt = .5
times = seq(0,5000,by=dt)
for (i in 1:length(times)){
  y = y + func(times[i],y,0)$dU*dt
  Y = rbind(Y, c(times[i],y))
}


matplot(x=X, y=t(Y[,-1]), type = "l", lty=1, col=colorRampPalette(colors = c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)), alpha=T)(nrow(Y)), xlab="Size", ylab="Density", main="FMU", log="xy")
points(x=X, y=Y[1,-1], col="grey", type="l", lwd=2)
points(x=x, y=Ueq(x), col="green3", type="l", lwd=3)
points(x=X, y=Y[nrow(Y),-1], col="blue", pch=20, cex=1)
abline(v=X, col="grey")


bo = numeric(nrow(Y))
boo = bo
for (i in 1:length(bo)){
  U = Y[i,-1]
  bo[i] = sum(birthRate(X,Y[i,1], calcEnv(x,U))*U*h)
  # Uo = (1-X)^2/(1+X)^4 + (1-X)/(1+X)^3*exp(-0.225*Y[i,1]^2)
  # boo[i] = sum(birthRate(X,Y[i,1], calcEnv(x, Uo))*Uo*h)
}

# plot(boo~Y[,1], type="l")
plot(bo~Y[,1], type="o")
