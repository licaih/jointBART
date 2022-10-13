
## BART: Bayesian Additive Regression Trees
## Copyright (C) 2017 Robert McCulloch and Rodney Sparapani

## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program; if not, a copy is available at
## https://www.R-project.org/Licenses/GPL-2

JointPBart=function(
    x.train, # list of xk matrix dim: n1 X p, n2 X p, ..., nK X p
    y.train, # list of yk vector dim: n1, n2, n3
    Theta,
    adj,
    graph_nu,
    graph_alpha,
    graph_beta,
    my_w,
    graph_a,
    graph_b,
    x.test=vector("list", length(y.train)), # ignored
    sparse=FALSE,
    theta=0,
    omega=1,
    a=0.5,
    b=1,
    augment=FALSE,
    rho=NULL,
    xinfo=vector("list", length(y.train)),
    usequants=FALSE,
    cont=FALSE,
    rm.const=TRUE,
    bk=2.0, # k
    power=2.0,
    base=.95,
    binaryOffset=NULL,
    ntree=200L,
    numcut=100L,
    ndpost=1000L,
    nskip=100L,
    transposed=FALSE,
    adj_alpha0=0.05,
    adj_alpha1=1.,
    Joint = T
)
{
#--------------------------------------------------
#data
    n = sapply(y.train, length)
    K = length(y.train)

    B = data.matrix(expand.grid(replicate(K, 0:1, simplify = FALSE)))


    if(length(binaryOffset)==0){
        binaryOffset = rep(NA,K)
        for(k in 1:K)
            binaryOffset[[k]]=qnorm(mean(y.train[[k]]))
    }

    if(!transposed) {
        for(k in 1:K){

            temp = bartModelMatrix(x.train[[k]], numcut, usequants=usequants,
                                   cont=cont, xinfo=xinfo[[k]], rm.const=rm.const)
            x.train[[k]] = t(temp$X)
            xinfo[[k]] = temp$xinfo
            if(length(x.test[[k]])>0) {  # ignored
                x.test[[k]] = bartModelMatrix(x.test[[k]])
                x.test[[k]] = t(x.test[[k]][ , temp$rm.const])
            }
        }
        numcut = temp$numcut
        rm.const <- temp$rm.const
        grp <- temp$grp
        rm(temp)

    }
    else {
        rm.const <- NULL
        grp <- NULL
    }

    if(!all(n==sapply(x.train, ncol)))
        stop('The length of y.train and the number of rows in x.train must be identical')

    p  = sapply(x.train, nrow)

    if(!all(p==p[1]))
        stop('The number of variables in x.train are not identical')

    p = p[1]

    np = try(sapply(x.test, ncol))
    if(class(np) != "integer") np = rep(0, K)
    if(length(rho)==0) rho <- p
    if(length(rm.const)==0) rm.const <- 1:p
    if(length(grp)==0) grp <- 1:p

#prior

#--------------------------------------------------
    ptm <- proc.time()

    #call
    res = JointBartB(
        n,  #number of observations in training data
        p,  #dimension of x
        np, #number of observations in test data
        x.train,   #pxn training data x
        y.train,   #pxn training data x
        x.test,   #p*np test data x
        ntree,
        numcut,
        ndpost,
        nskip,
        power,
        base,
        binaryOffset,
        rep(3/(bk*sqrt(ntree)),K), #tau
        sparse,
        theta,
        omega,
        grp,
        a,
        b,
        rho,
        augment,
        xinfo,
        Theta,
        adj,
        graph_nu,
        B,
        graph_alpha,
        graph_beta,
        my_w,
        graph_a,
        graph_b,
        adj_alpha0,
        adj_alpha1,
        Joint
    )
    res$proc.time <- proc.time()-ptm
    res$prob.train.mean = res$prob.train = res$prob.test = res$prob.test.mean= list()
    for(k in 1:K){
        res$yhat.train[,1:n[k],k] = res$yhat.train[,1:n[k],k]+binaryOffset[k]
        res$prob.train[[k]] = pnorm(res$yhat.train[,1:n[k],k])
        res$prob.train.mean[[k]] =  apply(res$prob.train[[k]], 2, mean)
        if(np[k]!=0){
            res$yhat.test[,1:np[k],k] = res$yhat.test[,1:np[k],k]+binaryOffset[k]
            res$prob.test[[k]]= pnorm(res$yhat.test[,1:np[k],k])
            res$prob.test.mean[[k]] = apply(res$prob.test[[k]], 2, mean)
        }
    }

    res$varcount.mean <- apply(res$varcount, c(3,2), mean)
    res$varprob.mean <- apply(res$varprob, c(3,2), mean)
    res$rm.const <- rm.const
    res$binaryOffset=binaryOffset
    attr(res, 'class') <- 'pbart'
    return(res)
}
