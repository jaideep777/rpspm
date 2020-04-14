rm(list=ls())

require(deSolve)

source("codes/pspm_package/Rscripts/wave_model.R")
# source("codes/pspm_package/Rscripts/test_model.R")


#### EBT core ####

calc_dU_EBT <- function(Nint,xint,N0,pi0, t, E){

  dx = growthRate(xint,t,E)
  dN = -mortalityRate(xint,t,E)*Nint

  grad_dx = 0.001
  growthGrad = (growthRate(xb+grad_dx,t,E) - growthRate(xb,t,E))/grad_dx
  mortGrad = (mortalityRate(xb+grad_dx,t,E) - mortalityRate(xb,t,E))/grad_dx
  
  x0 = xb + pi0/(N0+1e-12)
  
  dN0  = -mortalityRate(xb,t,E)*N0 - mortGrad*pi0 + sum(birthRate(c(x0,xint),t,E)*c(N0,Nint))
  dpi0 = growthRate(xb,t,E)*N0 + growthGrad*pi0 - mortalityRate(xb,t,E)*pi0

  # cat("rates", dpi0, dx, dN0, dN, "\n")
  c(dpi0, dx, dN0, dN)
}


addCohort_EBT <- function(Nint,xint,N0,pi0){
  
  # remove almost dead cohorts
  id = which(Nint<1e-10)
  if (length(id) > 0){
    # cat("id = ", id, Nint[id], "\n")
    xint = xint[-id]
    Nint = Nint[-id]
  }
  
  if (N0 > 0){
    N = c(0, N0, Nint)
    X = c(xb, xb+pi0/N0, xint)
  }else{
    N = c(N0, Nint)
    X = c(pi0, xint)
  }

  # return the result
  c(X,N)
  
}



# y = c(x,u)
func_EBT <- function(t, y, par){
  n = length(y)/2
  pi0 = y[1]
  xint = y[2:n]
  N0 = y[n+1]
  Nint = y[(n+2):(2*n)]
  # cat(x, "\n")
  
  E = calcEnv_EBT(pi0, xint, N0, Nint) 

  dy = calc_dU_EBT(Nint, xint, N0, pi0, t,E)
  list(dy=dy)
}

J = 25
xb = 0
xm = 1


#schedule = c(0, exp(seq(log(0.01), log(20), length.out=500)))
schedule = seq(0,10,0.02)

#x = c(xb, exp(seq(log(xb+0.01), log(xm), length.out=J))) # seq(0,1, length.out=J+1)
x = seq(xb, xm, length.out=J) # seq(0,1, length.out=J+1)
h = diff(x)
X = x[-1]-diff(x)/2

y0 = c(0, X, 0, calcIC(X)*h)
#     ^ pi0 xint N0 Nint
plot(x=1, y=NA, xlim=c(0,1), ylim=c(0,2), xlab="Size", ylab="Density", main="EBT")
# plot(x=1, y=NA, xlim=c(0,20), ylim=c(0,1))
for (i in 1:(length(schedule)-1)){
  # cat("y0 = ", y0, "\n")
  yt = lsoda(y0, c(schedule[i],schedule[i+1]), func_EBT, par=0, rtol=1e-6, atol=1e-6, hmax=1)
  # cat("y = ", y, "\n")
  
  y = yt[2,-1]
  n = length(y)/2
  pi0 = y[1]
  xint = y[2:n]
  N0 = y[n+1]
  Nint = y[(n+2):(2*n)]
  x0 = xb + pi0/(N0+1e-12)
  
  id = as.integer(cut(c(x0,xint), breaks = x))
  x1 = tapply(c(x0,xint), id, FUN = mean)
  u1 = tapply(c(N0,Nint), id, FUN = sum)
  h1 = diff(c(xb, x1[-1]-diff(x1)/2, xm))
  # if (i %% 5 == 0){
  #   points(c(x0,xint)~rep(schedule[i+1], length(xint)+1), pch=20, cex=c(N0,Nint))
  # }
  cat(schedule[i+1],"\n")
  # cat(pi0, N0, "\n")
  
  y0 = addCohort_EBT(Nint, xint, N0, pi0)
  
  if ((i-1) %% 20 == 0){
    # cat("X:", c(x0, xint), "\n")
    # cat("N:", c(N0, Nint), "\n")
    points(y=u1/h1, x=x1, type="l", col=colorRampPalette(c(rgb(1,0,0,0.5),rgb(0,0,1,0.5)),alpha=T)(length(schedule))[i])
  }
}
points(x=x, y=Ueq(x), col="green3", type="l", lwd=3)
points(x=x, y=calcIC(x), col="grey", type="l", lwd=1)
points(y=u1/h1, x=x1, col="blue", pch=20)


