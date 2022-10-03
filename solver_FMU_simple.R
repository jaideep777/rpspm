# rm(list=ls())

require(deSolve)

#### Load model ####

# source("codes/pspm_package/Rscripts/wave_model.R")


#### FMU core ####
step_FMU <- function(U, x, t, E, dt){
  cat("step_FMU: ", t, "-->", t+dt, "(", dt, ")", "\n")
  X = x[-1]-diff(x)/2             # grid centers
  h = diff(x)

  birthArray = birthRate(X,t,E)
  growthArray = growthRate(X,t,E)
  mortalityArray = mortalityRate(X,t,E)
  
  birthFlux = sum(birthArray*U*h)
  
  Unew = U*0
  
  B1 = 1 + dt/h[1]*growthArray[1] + dt*mortalityArray[1]
  C1 = U[1]+ dt/h[1]*birthFlux
  Unew[1] = C1/B1
  
  for (w in 2:J){
    Aw = -growthArray[w-1]*dt/h[w]
    Bw = 1 + dt/h[w]*growthArray[w] + dt*mortalityArray[w]
    Cw = U[w]
    
    Unew[w] = (Cw - Aw*Unew[w-1])/Bw
  }
  
  Unew
}

step_to <- function(t0, tf, U0, stepsize){
  U = U0
  t = t0
  while(t<tf){
    dt_now = min(stepsize, tf-t)
    E = calcEnv(x, U)
    U = step_FMU(U,x,t,E, dt_now)
    t = t+dt_now
  }
  U
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
# tf_vec = seq(0,8,length.out=20)
# 
# y = matrix(nrow = length(tf_vec), ncol=1+length(U))
# y[1,2:ncol(y)]=U
# y[1,1] = tf_vec[1]
# 
# y0 = y
# bo = tf_vec*0
# boo = tf_vec*0
# for (i in 2:length(tf_vec)){
#   cat("stepping ", tf_vec[i-1], "---->" ,tf_vec[i], "\n")
#   y[i,1] = tf_vec[i]
#   U = step_to(tf_vec[i-1], tf_vec[i], U, stepsize = .01)
#   y[i,2:ncol(y)] = U
#   Uo = (1-X)^2/(1+X)^4 + (1-X)/(1+X)^3*exp(-0.225*tf_vec[i]^2)
#   y0[i,2:ncol(y0)] = Uo
#   bo[i] = sum(birthRate(X,tf_vec[i], calcEnv(x, U))*U*h)
#   boo[i] = sum(birthRate(X,tf_vec[i], calcEnv(x, Uo))*Uo*h)
#   
# }
# 
# 
# Y=y
# matplot(x=X, y=t(Y[,-1]), type = "l", lty=1, col=colorRampPalette(colors = c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)), alpha=T)(20), xlab="Size", ylab="Density", main="FMU")
# matlines(x=X, y=t(y0[,-1]), type="p", pch=20, lty=1, col=colorRampPalette(colors = c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)), alpha=T)(20), xlab="Size", ylab="Density", main="FMU")
# points(x=X, y=Y[1,-1], col="grey", type="l", lwd=2)
# points(x=x, y=Ueq(x), col="green3", type="l", lwd=3)
# points(x=X, y=Y[nrow(Y),-1], col="blue", pch=20, cex=1)
# abline(v=X, col="grey")
# 
# plot(boo~tf_vec, type="l")
# points(bo~tf_vec, col="red")




#### Solve PSPM ####
source("RED_model_v2.R")

J = 100                          # Number of grid cells
xb = 1                        # Minimum size (size at birth)
xm = 1e6                         # Maximum size (size at maturity)
# if logarithmiz grid desired:
x = exp(seq(log(xb), log(xm), length.out=J+1)) # seq(0,1, length.out=J+1)
# x = seq(xb, xm, length.out=J+1) # grid edges
X = x[-1]-diff(x)/2             # grid centers
h = diff(x)                     # grid widths


U = calcIC(X)


tf_vec = seq(0,5000,length.out=20)

y = matrix(nrow = length(tf_vec), ncol=1+length(U))
y[1,2:ncol(y)]=U
y[1,1] = tf_vec[1]

for (i in 2:length(tf_vec)){
  cat("stepping ", tf_vec[i-1], "---->" ,tf_vec[i], "\n")
  y[i,1] = tf_vec[i]
  U = step_to(tf_vec[i-1], tf_vec[i], U, stepsize = 1)
  y[i,2:ncol(y)] = U
}


Y = y
matplot(x=X, y=t(Y[,-1]), type = "l", lty=1, col=colorRampPalette(colors = c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)), alpha=T)(20), xlab="Size", ylab="Density", main="FMU", log="xy")
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
points(bo~Y[,1], type="o", col="orange")


