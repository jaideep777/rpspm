#### ~~~~~ Model description ~~~~~~~~~~~ ####

calcIC <- function(x){
  u = numeric(length(x))
  u[x<.2]=1
  u[x>=.2]=0
  # u[x>.8]=1
  # u[x<=.8]=0
  u
}

growthRate = function(x,t,E){
  rep(.1, length(x))
}

mortalityRate = function(x,t,E){
  rep(0, length(x))
}

birthRate = function(x,t,E){
  rep(0, length(x))
}

calcEnv <- function(x,U){
  rep(0, length(U))
}

calcEnv_CM <- function(x,u){
  rep(0, length(u))
}

calcEnv_EBT <- function(pi0,xint,N0,Nint){
  rep(0, length(xint)+1)
}

Ueq = function(x){
 x*0
}

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
