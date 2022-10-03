rm(list=ls())

require(deSolve)

#### Load model ####

# source("codes/pspm_package/Rscripts/wave_model.R")
source("~/codes/rpspm/test_model.R")


#### CM core ####

integrate_trapezium <- function(x,y){
  idx = 2:length(x)
  return (as.double( (x[idx] - x[idx-1]) %*% (y[idx] + y[idx-1])) / 2)
}


calc_dU_CM <- function(u, x, t, E){
  dx = growthRate(x,t,E)
  grad_dx = 0.001
  growthGrad = (growthRate(x+grad_dx,t,E) - growthRate(x,t,E))/grad_dx
  du = -mortalityRate(x,t,E)*u - growthGrad*u
  c(dx, du)
}


addCohort <- function(u,x,t,E, B){
  
  # # predicted u0 as a function of u(x) including u0
  # f <- function(unew,x,t){
  #   Enew = calcEnv_CM(c(xb,x),c(unew,u))
  #   b = birthRate(c(xb,x),t,Enew)
  #   birthFlux = integrate_trapezium(c(xb,x), b*c(unew,u))
  #   u0 = birthFlux/(growthRate(xb,t,Enew)+1e-12)
  #   u0
  # }
  # 
  # # Iterate to solve for f(u0)=u0
  # u0 = 10
  # err=100
  # while(err > 1e-4){
  #   cat(".")
  #   u1 = f(u0,x,t)
  #   err = abs(u1-u0)
  #   u0 = u1
  # }
  # cat("\n")
  
  # remove most crowded cohort
  J = length(x)
  Dx = rep(NA,J)
  for (j in 2:(J-1)){
    Dx[j] = x[j+1]-x[j-1]
  }
  # cat("Dx = ", Dx, "\n")
  id = which(Dx == min(Dx, na.rm=T))
  x = x[-id]
  u = u[-id]
  
  # # Remove almost dead cohorts  
  # id = which(u<1e-6)
  # if (length(id) > 0){
  #   x = x[-id]
  #   u = u[-id]
  # }
  
  insert_lower = T
  # return the result
  if (insert_lower){
    c(xb,x,B,u)
  }
  else{
    c(x,xm,u,B)
  }
}


# y = c(x,u)
func_CM <- function(t, y, par){
  n = length(y)/2
  x = y[1:n]
  u = y[(n+1):(2*n)]
  # cat(x, "\n")
  
  E = calcEnv_CM(x,u)
  
  dxdu = calc_dU_CM(u,x,t,E)
  list(dxdu=dxdu)
}

#### Solve PSPM ####

J = 25
xb = 0
xm = 1

schedule = seq(0,8,.05)  # New cohort introduction times

x = seq(0,1, length.out=J+1)
x2 = x
u0 = c(x, calcIC(x))
B = 2
plot(x=1, y=NA, xlim=c(0,1), ylim=c(0,2), xlab="Size", ylab="Density", main="CM")
# plot(x=1, y=NA, xlim=c(0,100), ylim=c(0,1))
for (i in 1:(length(schedule)-1)){
  ut = lsoda(u0, c(schedule[i],schedule[i+1]), func_CM, par=0, rtol=1e-6, atol=1e-6)
  y = ut[2,-1]
  
  # y = u0 + unlist(func_CM(schedule[i+1], u0, 0))*(schedule[i+1]-schedule[i])
  
  n = length(y)/2
  x1 = y[1:n]
  u1 = y[(n+1):(2*n)]
  # if (i %% 200 == 0){
  #   points(x1~rep(schedule[i+1], length(x1)), pch=20, cex=u1)
  # }

  Enew = calcEnv_CM(c(xb,x1), c(B,u1))
  b = birthRate(c(xb,x1),schedule[i+1],Enew)
  birthFlux = integrate_trapezium(c(xb,x1), b*c(B,u1))
  cat("u0 = ", birthFlux/(growthRate(xb,schedule[i+1],Enew)+1e-12), "\n")
    
  cat(schedule[i+1],"\n")
  u0 = addCohort(u1,x1,schedule[i+1],calcEnv_CM(x1,u1), B)
  
  
  if (i %% 10 == 0){
    n = length(u0)/2
    x = u0[1:n]
    u = u0[(n+1):(2*n)]
    points(u~x, type="l", col=colorRampPalette(c(rgb(1,0,0,0.5),rgb(0,0,1,0.5)),alpha=T)(length(schedule))[i])
  }
}

points(x=x2, y=calcIC(x2), col="grey", type="l", lwd=1)
points(x=x2, y=Ueq(x2), col="green3", type="l", lwd=3)
points(u~x, col="blue", pch=20)


