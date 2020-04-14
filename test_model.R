#### ~~~~~ Model description ~~~~~~~~~~~ ####

calcIC <- function(x){
  IC = (1-x)^2/(1+x)^4 + (1-x)/(1+x)^3
}

growthRate = function(x,t,E){
  a = 0.16+0.22*exp(-0.225*t*t)
  g = 0.225*(1-x*x)*(E/(1+E*E))*t*(1+a*a)/a
  g
}

mortalityRate = function(x,t,E){
  a = 0.16+0.22*exp(-0.225*t*t)
  m = 1.35*t*E/a
  rep(m, length(x))
}

# Note: total flux of newborns is calculated in the solver as
#       g(xb,E)u(xb,t) = integral(birthRate(x,t,E)*u(x,t))
birthRate = function(x,t,E){
  a = 0.16+0.22*exp(-0.225*t*t)
  oneplusa = 1.16+0.22*exp(-0.225*t*t)
  n1 = 0.225*t*x*x*(1-x)*(1-x)*E/(1+E)/(1+E)*oneplusa*oneplusa/a
  n2 = (1+exp(-0.225*t*t))/(61-88*log(2)+(38*log(2)-79.0/3)*exp(-0.225*t*t))
  n1*n2
}

calcEnv <- function(x,U){
  w = (2-3*x)^3 * (54*x^2 - 27*x +4)
  w[which(x <= 1/3)] = 1
  w[which(x >  2/3)] = 0
  
  h = diff(x)
  w = w[-1]-diff(w)/2
  
  E = sum(h*w*U)
  E
}

calcEnv_CM <- function(x,u){
  w = (2-3*x)^3 * (54*x^2 - 27*x +4)
  w[which(x <= 1/3)] = 1
  w[which(x >  2/3)] = 0

  if (length(x) < 2){
    E = 0
  }
  else{
    E = integrate_trapezium(x, w*u)
  }
  E
}

calcEnv_EBT <- function(pi0,xint,N0,Nint){
  x0 = xb+pi0/(N0+1e-12)
  
  x = c(x0, xint)
  N = c(N0, Nint)
  
  w = (2-3*x)^3 * (54*x^2 - 27*x +4)
  w[which(x <= 1/3)] = 1
  w[which(x >  2/3)] = 0
  
  E = sum(w*N)
  E
}


# Analytical calculation of equilibrium distribution, if available
Ueq = function(x){
  (1-x)^2/(1+x)^4
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
