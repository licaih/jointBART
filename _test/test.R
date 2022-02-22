#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'


f = function(x){
  10*sin(pi*x[,1]*x[,2]) + 20*(x[,3]-.5)^2+10*x[,4]+5*x[,5]
}
sigma = 1.0  #y = f(x) + sigma*z , z~N(0,1)
n = 100      #number of observations
set.seed(99)
x=matrix(runif(n*10),n,10) #10 variables, only first 5 matter
Ey = f(x)
y=Ey+sigma*rnorm(n)

temp = bartModelMatrix(x, 100, usequants=F,
                       cont=F, xinfo=matrix(0.0,0,0), rm.const=T)
temp


n = c(40, 60, 100)

set.seed(2022)
x  = list(matrix(runif(n[1]*10),n[1],10),
          matrix(runif(n[2]*10),n[2],10),
          matrix(runif(n[3]*10),n[3],10))
x
y  = list(rnorm(n[1]), rnorm(n[2]), rnorm(n[3]))
y
xp = list(matrix(runif(n[1]*10),n[1],10),
          matrix(runif(n[2]*10),n[2],10),
          matrix(runif(n[3]*10),n[3],10))
xp

K = length(y)
K
xinfo=vector("list", K)

for(k in 1:K){

  temp = bartModelMatrix(x[[k]], 10, usequants=F,
                         cont=F, xinfo=xinfo[[k]], rm.const=T)
  x[[k]] = t(temp$X)
  numcut = temp$numcut
  xinfo[[k]] = temp$xinfo
  if(length(xp)>0) {  # ignored
    xp[[k]] = bartModelMatrix(xp[[k]])
    xp[[k]] = t(xp[[k]][ , temp$rm.const])
  }
  rm.const <- temp$rm.const
  grp <- temp$grp
  rm(temp)
}

w = lapply(n, function(n1) rep(1,n1))
w
## --------------------------
## Test
## --------------------------
set.seed(2022)
JointBart(n = n,
          p = 10,
          np = n, #number of observations in test data
          x = x,   #pxn training data x
          y = y,   #pxn training data x
          xp = xp,   #p*np test data x
          m = 20,
          nc = numcut,
          nd = 50,
          burn = 50,
          mybeta = 2.0,
          alpha = 0.95,
          tau = rep(0.05, 3),
          nu =  rep(0.05, 3),
          lambda = rep(0.1, 3),
          sigma = c(0.1,0.15,0.2),
          w = w,
          dart = F,
          theta = 0,
          omega = 1.0,
          igrp = 1:10,
          a = 0.5,
          b = 1.0,
          rho = 3.0,
          aug = T,
          iXinfo = xinfo
          )


