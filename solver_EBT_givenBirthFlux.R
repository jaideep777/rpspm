rm(list=ls())

require(deSolve)

setwd("~/codes/rpspm")
# source("codes/pspm_package/Rscripts/wave_model.R")
source("test_model.R")
# source("RED_model_v2.R")


#### EBT core ####

calc_dU_EBT <- function(Nint,xint,N0,pi0, t, E, B){

  dx = growthRate(xint,t,E)
  dN = -mortalityRate(xint,t,E)*Nint

  x0 = xb + pi0/(N0+1e-12)

  grad_dx = 0.001 #max(x0-xb, 0.0001)
  growthGrad = (growthRate(xb+grad_dx,t,E) - growthRate(xb,t,E))/grad_dx
  mortGrad = (mortalityRate(xb+grad_dx,t,E) - mortalityRate(xb,t,E))/grad_dx
  
  # u0 = sum(birthRate(c(x0,xint),t,E)*c(N0,Nint))/growthRate(xb,t,E)
  
  birthFlux = B*growthRate(xb,t,E)
  # cat(t, " ", B, " -> ", u0, "\n")
  
  dN0  = -mortalityRate(xb,t,E)*N0 - mortGrad*pi0 + birthFlux
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
  # cat(E, "\n")

  dy = calc_dU_EBT(Nint, xint, N0, pi0, t,E, par)
  list(dy=dy)
}

J = 25
xb = 0
xm = 1


#schedule = c(0, exp(seq(log(0.01), log(20), length.out=500)))
schedule = seq(0,8,.05)

x = seq(0,1, length.out=J+1)
# x = exp(seq(log(xb), log(xm), length.out=J)) # seq(0,1, length.out=J+1)
h = diff(x)
X = x[-1]-diff(x)/2

y0 = c(0, X, 0, calcIC(X)*h)
#     ^ pi0 xint N0 Nint

plot(x=1, y=NA, xlim=c(0,1), ylim=c(0,2), xlab="Size", ylab="Density", main="EBT")
# plot(x=1, y=NA, xlim=c(1,1e4), ylim=c(1e-14,1e4), xlab="Size", ylab="Density", main="EBT", log="xy")
# plot(x=1, y=NA, xlim=c(0,20), ylim=c(0,1))
for (i in 1:(length(schedule)-1)){
  # cat("y0 = ", y0, "\n")
  yt = lsoda(y0, c(schedule[i],schedule[i+1]), func_EBT, par=2, rtol=1e-6, atol=1e-6, hmax=1)
  # cat("y = ", y, "\n")
  
  y = yt[2,-1]
  n = length(y)/2
  pi0 = y[1]
  xint = y[2:n]
  N0 = y[n+1]
  Nint = y[(n+2):(2*n)]
  x0 = xb + pi0/(N0+1e-12)
  
  id = as.integer(cut(c(x0,xint), breaks = x))
  id# x1 = tapply(c(x0,xint), id, FUN = mean)
  x1_1 = tapply(c(x0,xint)*c(N0,Nint), id, FUN = sum)
  x1_2 = tapply(c(N0,Nint), id, FUN = sum)
  x1 = x1_1/x1_2
  u1 = tapply(c(N0,Nint), id, FUN = sum)
  h1 = diff(c(xb, x1[-1]-diff(x1)/2, xm))
  # if (i %% 5 == 0){
  #   points(c(x0,xint)~rep(schedule[ i+1], length(xint)+1), pch=20, cex=c(N0,Nint))
  # }
  # cat(schedule[i+1],"\n")
  # cat(pi0, N0, "\n")
  
  y0 = addCohort_EBT(Nint, xint, N0, pi0)
    
  E = calcEnv_EBT(pi0, xint, N0, Nint)
  u0 = sum(birthRate(c(x0,xint),schedule[i+1],E)*c(N0,Nint))/growthRate(xb,schedule[i+1],E)
  cat(schedule[i+1], " 2 ->", u0, "\n")

  if ((i-1) %% 20 == 0){
    # cat("X:", c(x0, xint), "\n")
    # cat("N:", c(N0, Nint), "\n")
    points(y=u1/h1, x=x1, type="l", col=colorRampPalette(c(rgb(1,0,0,0.5),rgb(0,0,1,0.5)),alpha=T)(length(schedule))[i])
  }
  
}
points(x=x, y=Ueq(x), col="green3", type="l", lwd=3)
points(x=x, y=calcIC(x), col="grey", type="l", lwd=1)
points(y=u1/h1, x=x1, col="blue", pch=20)
abline(v=c(xb, x1[-1]-diff(x1)/2, xm), col="grey")

# 
# 
# # #### C++ with cohortsToDensity ####
file = "~/codes/pspm_package/ebt_testmodel.txt"
n_col <- max(count.fields(file, sep = "\t"))
YY <- readLines(file)
plot(x=1, y=NA, xlim=c(0,1), ylim=c(0,2), xlab="Size", ylab="Density", main="EBT")
for (i in 1:(length(schedule)-1)){
  yt = as.numeric(unlist(strsplit(YY[i], split = "\t")))

  y = yt[-c(1,2)]
  n = length(y)/2

  if ((i-1) %% 20 == 0){
    points(x=y[1:n], y=y[(n+1):(2*n)], type="l", col=colorRampPalette(c(rgb(1,0,0,0.5),rgb(0,0,1,0.5)),alpha=T)(length(schedule))[i])
  }
}
points(x=x, y=Ueq(x), col="green3", type="l", lwd=3)
points(x=x, y=calcIC(x), col="grey", type="l", lwd=1)
points(x=y[1:n], y=y[(n+1):(2*n)], col="blue", pch=20)
abline(v=c(xb, x1[-1]-diff(x1)/2, xm), col="grey")

# 
# dat = read.delim("/home/jaideep/codes/pspm_package/ebt_testmodel.txt", header=F)
# dat=dat[,-ncol(dat)]
# dat[dat==0] = NA
# matplot(x=X, y=t(as.matrix(dat[,c(-1,-2)])), type="l", lty=1, col=colorRampPalette(c(rgb(1,0,0,0.5),rgb(0,0,1,0.5)),alpha=T)(nrow(dat)))
# matplot(x=X, y=t(as.matrix(dat[160,c(-1,-2)])), pch=20, col=colorRampPalette(c(rgb(1,0,0,0.5),rgb(0,0,1,0.5)),alpha=T)(nrow(dat)), add=T)
