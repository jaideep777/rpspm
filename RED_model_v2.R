require(pracma)

#### ~~~~~ Model description ~~~~~~~~~~~ ####

calcIC <- function(x){
  IC = 100/x^4
  
  IC
}

growthRate = function(x,t,E){
  g = 0.0838*x^0.7134
  g
}

mortalityRate = function(x,t,E){
  m = 0.035
  rep(m, length(x))
}

# Note: total flux of newborns is calculated in the solver as
#       g(xb,E)u(xb,t) = integral(birthRate(x,t,E)*u(x,t))
birthRate = function(x,t,E){
  0.1/0.9*0.0838*x^0.7134*(1-E)
}

calcEnv <- function(x,U){
  X = x[-1]-diff(x)/2         
  h = diff(x)
  E=min(sum(0.396*X^0.749*U*h/10000),1)
  E = max(E,0)
  # print(E)
  E
}

calcEnv_CM <- function(x,u){
  
  w = 0.396*x^0.749/10000
  if (length(x) < 2){
    E = 1
  }
  else{
    E = min(integrate_trapezium(x, w*u),1)
  }
  E
}

calcEnv_EBT <- function(pi0,xint,N0,Nint){
  x0 = xb+pi0/(N0+1e-12)
  
  x = c(x0, xint)
  N = c(N0, Nint)
  
  w = 0.396*x^0.749/10000  
  E = sum(w*N)
  min(E,1)
}


# Analytical calculation of equilibrium distribution, if available
Ueq = function(x){
  a0<- 0.396
  phiA<- 0.749
  g0<-0.0838
  phiG<-0.7134
  m0<-1
  mort<-0.035
  mu0<-mort*m0/g0
  alpha<-0.10
  temp<- mu0/(1-phiG)
  coverage<-1-(1-alpha)/alpha*mu0/((mu0/(1-phiG))^(phiG/(phiG-1))*exp(mu0/(1-phiG))*as.numeric(gammainc(mu0/(1-phiG),phiG/(1-phiG)+1)[2]))
  Neq<-coverage/a0/(temp^(phiA/(phiG-1))*exp(temp)*as.numeric(gammainc(temp,phiA/(1-phiG)+1))[2])
  n0<-Neq*mort/g0
  return(n0*(x/m0)^(-phiG)*exp(mu0/(1-phiG)*(1-(x/m0)^(1-phiG)))*10000)
}


# calcEnv_int <- function(X,U,h){
#   uint = splinefun(x = X, y=U, method = "natural")
#   
#   w <- function(x1){
#     out = 0
#     if (x1<1/3) out = 1*uint(x1)
#     else if (x1 >=1/3 & x1 < 2/3) out = (2-3*x1)^3 * (54*x1^2 - 27*x1 +4) * uint(x1)
#     else out = 0 
#     out
#   }
#   
#   E = integrate(f = function(x){sapply(x,w)}, lower = 0, upper = 1)$value
#   E
# }



### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
