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

temp = bartModelMatrix(x, 100, usequants=F,
                       cont=F, xinfo=matrix(0.0,0,0), rm.const=T)
temp


n = c(100, 100, 100)

set.seed(2022)
x  = list(matrix(runif(n[1]*10),n[1],10),
          matrix(runif(n[2]*10),n[2],10),
          matrix(runif(n[3]*10),n[3],10))
x
y  = list(f(x[[1]]) + rnorm(n[1]),
          f(x[[2]]) + rnorm(n[2]),
          f(x[[2]]) + rnorm(n[3]))
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
## Test C
## --------------------------
set.seed(2022)
res=JointBart(n = n,
          p = 10,
          np = n, #number of observations in test data
          x = x,   #pxn training data x
          y = y,   #pxn training data x
          xp = xp,   #p*np test data x
          m = 20,
          nc = numcut,
          nd = 15000,
          burn = 5000,
          mybeta = 2.0,
          alpha = 0.95,
          tau = rep(10, 3),
          nu =  rep(10, 3),
          lambda = rep(10, 3),
          sigma = c(1,1,1),
          w = w,
          dart = F,
          theta = 0,
          omega = 1.0,
          igrp = 1:10,
          a = 0.5,
          b = 1.0,
          rho = 1.0,
          aug = T,
          iXinfo = xinfo
          )

plot(res$sigma[,1])
plot(res$sigma[,2])
plot(res$sigma[,3])
apply(res$varcount[,,3],2, mean)
apply(res$varcount[,,2],2, mean)
apply(res$varcount[,,1],2, mean)

## -------------------------------
## Test R wrapper
## -------------------------------
# x.train = x; y.train = y; x.test=xp;
# sparse=FALSE;theta=0;omega=1; a=0.5;b=1;augment=FALSE; rho=NULL;
# xinfo=vector("list", K); usequants=FALSE; cont=FALSE;
# rm.const=TRUE;sigest=NA;sigdf=3; sigquant=.90;bk=2.0;
# power=2.0;base=.95;sigmaf=NA;lambda=NA;fmean=lapply(y.train, mean);
# w=lapply(n, function(n1) rep(1,n1));ntree=20L;numcut=30;ndpost=5000L;nskip=5000L
# transposed=FALSE

f = function(x){
  10*sin(pi*x[,1]*x[,2]) - 10*(x[,3]-.5)^2+10*x[,4]-5*x[,5]
}
sigma = 1.0  #y = f(x) + sigma*z , z~N(0,1)
n = c(100, 100, 100)
p = 50
# set.seed(2022)
x  = list(matrix(runif(n[1]*p),n[1],p),
          matrix(runif(n[2]*p),n[2],p),
          matrix(runif(n[3]*p),n[3],p))
x
y  = list(f(x[[1]]) + rnorm(n[1]),
          f(x[[2]]) + rnorm(n[2]),
          f(x[[3]]) + rnorm(n[3]))
y

K = length(y)
K

x.train = x; y.train = y;
Theta = matrix(0, 3, 3)
Theta[2,3] = Theta[3,2] = 1
Theta[1,3] = Theta[3,1] = 1
Theta[1,2] = Theta[2,1] = 1
adj   = matrix(rbinom(p*K, 1, 0), p, 3)
graph_nu = rgamma(p, 4, 2)
res = JointWBart(x.train, y.train,
                 Theta, adj, graph_nu, graph_alpha = 0.5,
                 graph_beta=0.5, my_w = 0.1, graph_a = 4, graph_b = 2,
                 x.test=vector("list", length(y.train)),
  sparse=FALSE, theta=0, omega=1, a=0.5, b=1, augment=FALSE, rho=NULL,
  xinfo=vector("list", length(y.train)), usequants=FALSE, cont=FALSE,
  rm.const=TRUE, sigest=NA, sigdf=3, sigquant=.90, bk=2.0, power=2.0, base=.95,
  sigmaf=NA, lambda=NA, fmean=lapply(y.train, mean),
  w=lapply(n, function(n1) rep(1,n1)),  ntree=10L, numcut=100, ndpost=1000L,
  nskip=30000L, transposed=FALSE
)

pip = apply(res$varcount>0, c(3,2), mean)
col = c(rep(2, 5), rep(1, p-5))
plot(pip[1,],  ylim =c(0,1), col = col, pch = 21)
points(pip[2,], col = col, pch = 22)
points(pip[3,], col = col, pch = 23)
abline(h = 0.5, col = 4, lty = 3)
legend("topright", legend = 1:3, pch = 21:23)

pip
apply(res$varcount, c(3,2), mean)

plot(res$sigma[,1], type = "l")
plot(res$sigma[,2], type = "l")
plot(res$sigma[,3], type = "l")

plot(y[[1]], apply(res$yhat.train[,1:n[1],1],2,mean)+res$mu[[1]])
abline(0,1,col = 2)
plot(y[[2]], apply(res$yhat.train[,1:n[2],2],2,mean)+res$mu[[2]])
abline(0,1,col = 2)
plot(y[[3]], apply(res$yhat.train[,1:n[3],3],2,mean)+res$mu[[3]])
abline(0,1,col = 2)


