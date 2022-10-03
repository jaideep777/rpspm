rm(list=ls())
source("~/codes/rpspm/test_model.R")


Ut <- function(x,t){
  (1-x)^2/(1+x)^4 + (1-x)/(1+x)^3 * exp(-0.225*t^2)
}
xhires = seq(0,1,length.out=100)

schedule = seq(0.05,8, 0.05)

par(mfrow = c(2,2), mar=c(4,4,1,1), oma=c(1,1,1,1))


#### C++ FMU ####


Y = read.delim("~/codes/pspm_package/fmu_testmodel.txt", header=F)
Y=Y[,-length(Y)]

J = ncol(Y)-2                  # Number of grid cells
xb = 0                         # Minimum size (size at birth)
xm = 1                         # Maximum size (size at maturity)
# if logarithmiz grid desired:
# x = exp(seq(log(xb), log(xm), length.out=J+1)) # seq(0,1, length.out=J+1)
x = seq(xb, xm, length.out=J+1) # grid edges
X = x[-1]-diff(x)/2             # grid centers
h = diff(x)                     # grid widths

B_fmu = Y[,2]

matplot(x=X, y=t(Y[seq(1,nrow(Y),10),c(-1,-2)]), ylim=c(0,2), cex=.5, pch=20, col=colorRampPalette(colors = c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)), alpha=T)(20), xlab="Size", ylab="Density", main="FMU C++")
for (i in 1:nrow(Y)){
  if (i %% 10 == 1) points(Ut(xhires,Y[i,1])~xhires, type = "l", lty=1, col=colorRampPalette(colors = c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)), alpha=T)(20)) 
}
points(x=xhires, y=Ueq(xhires), col=scales::alpha("green3", .8), type="l", lwd=2)
points(x=X, y=Y[1,c(-1,-2)], col="grey", type="l", lwd=2)
points(x=X, y=Y[nrow(Y),c(-1,-2)], col="blue", pch=20, cex=1)
# abline(v=X, col="grey")


#### CM C++ ####

YY <- readLines("~/codes/pspm_package/cm_testmodel.txt")
B_cm = numeric(length(YY))
plot(x=1, y=NA, xlim=c(0,1), ylim=c(0,2), xlab="Size", ylab="Density", main="CM C++")
# plot(x=1, y=NA, xlim=c(0,100), ylim=c(0,1))
for (i in 1:(length(YY))){
  yt = as.numeric(unlist(strsplit(YY[i], split = "\t")))

  t = yt[1]
  B_cm[i] = yt[2]
  y = yt[c(-1,-2)]
  n = length(y)/2
  cat("n = ", n, "\n")
  x1 = y[1:n]
  u1 = exp(y[(n+1):(2*n)])
  # if (i %% 200 == 0){
  #   points(x1~rep(schedule[i+1], length(x1)), pch=20, cex=u1)
  # }
  cat(schedule[i+1],"\n")
  # u0 = addCohort(u1,x1,schedule[i+1],calcEnv_CM(x1,u1))
  
  if (i %% 10 == 0){
    # n = length(u0)/2
    # x = u0[1:n]
    # u = u0[(n+1):(2*n)]
    points(u1~x1, pch=20, cex=.5, col=colorRampPalette(c(rgb(1,0,0,0.5),rgb(0,0,1,0.5)),alpha=T)(length(schedule))[i])
    points(Ut(xhires,t)~xhires, type="l", col=colorRampPalette(c(rgb(1,0,0,0.5),rgb(0,0,1,0.5)),alpha=T)(length(schedule))[i])
  }
}
# points(x=x1, y=calcIC(x1), col="grey", type="l", lwd=1)
points(x=xhires, y=Ueq(xhires), col=scales::alpha("green3", .8), type="l", lwd=2)
points(u1~x1, col="blue", pch=20)




#### C++ EBT ####

file = "~/codes/pspm_package/ebt_testmodel.txt"
YY <- readLines(file)
B_ebt = numeric(length(YY))
plot(x=1, y=NA, xlim=c(0,1), ylim=c(0,2), xlab="Size", ylab="Density", main="EBT C++")
for (i in 1:length(YY)){
  yt = as.numeric(unlist(strsplit(YY[i], split = "\t")))
  
  B_ebt[i] = yt[2]
  y = yt[c(-1,-2)]
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
  
  if ((i-1) %% 10 == 0){
    # cat("X:", c(x0, xint), "\n")
    # cat("N:", c(N0, Nint), "\n")
    points(y=u1/h1, x=x1, pch=20, cex=.5, type="o", lwd=.5, col=colorRampPalette(c(rgb(1,0,0,0.5),rgb(0,0,1,0.5)),alpha=T)(length(schedule))[i])
    points(y=Ut(xhires, yt[1]), x=xhires, type="l", col=colorRampPalette(c(rgb(1,0,0,0.5),rgb(0,0,1,0.5)),alpha=T)(length(schedule))[i])
  }
}
points(x=xhires, y=Ueq(xhires), col=scales::alpha("green3", .8), type="l", lwd=2)
points(x=x, y=calcIC(x), col="grey", type="l", lwd=1)
points(y=u1/h1, x=x1, col="blue", pch=20)
# abline(v=c(xb, x1[-1]-diff(x1)/2, xm), col="grey")



matplot(x=seq(0.05, 8, 0.05), y=t(rbind(B_fmu, B_ebt, B_cm)), type="l", lty=1, col=colorRampPalette(colors = c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)), alpha=T)(3), xlab="time", ylab="number", main="Newborns")


