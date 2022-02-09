
JointWBart=function(
  x.train, # list of xk matrix dim: n1 X p, n2 X p, ..., nK X p
  y.train, # list of yk vector dim: n1, n2, n3
  x.test=matrix(0.0,0,0), # ignored
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
  sigest=NA,
  sigdf=3,
  sigquant=.90,
  bk=2.0, # k
  power=2.0,
  base=.95,
  sigmaf=NA,
  lambda=NA,
  fmean=lapply(y.train, mean),  # vector of yk mean dim: n1, n2, n3
  w=lapply(y.train, function(y) rep(1,length(y))),
  ntree=200L,
  numcut=100L,
  ndpost=1000L,
  nskip=100L,
  keepevery=1L,
  nkeeptrain=ndpost,
  nkeeptest=ndpost,
  nkeeptestmean=ndpost,
  nkeeptreedraws=ndpost,
  printevery=100L,
  transposed=FALSE
)
{
  #--------------------------------------------------
  #data
  n = sapply(y.train, length)
  K = length(y.train)
  if(!transposed) {
    for(k in 1:length){

    temp = bartModelMatrix(x.train[[k]], numcut, usequants=usequants,
                           cont=cont, xinfo=xinfo[[k]], rm.const=rm.const)
    x.train[[k]] = t(temp$X)
    numcut = temp$numcut
    xinfo[[k]] = temp$xinfo
    if(length(x.test)>0) {  # ignored
      x.test[[k]] = bartModelMatrix(x.test[[k]])
      x.test[[k]] = t(x.test[[k]][ , temp$rm.const])
    }
    rm.const <- temp$rm.const
    grp <- temp$grp
    rm(temp)
    }
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

  np = sapply(x.text, ncol)
  if(length(rho)==0) rho=p
  if(length(rm.const)==0) rm.const <- 1:p
  if(length(grp)==0) grp <- 1:p

  ##if(p>1 & length(numcut)==1) numcut=rep(numcut, p)

  for(k in 1:K){
    y.train[[1]] = y.train[[1]]-fmean[[1]]
  }
  #--------------------------------------------------
  #set nkeeps for thinning
  if((nkeeptrain!=0) & ((ndpost %% nkeeptrain) != 0)) {
    nkeeptrain=ndpost
    cat('*****nkeeptrain set to ndpost\n')
  }
  if((nkeeptest!=0) & ((ndpost %% nkeeptest) != 0)) {
    nkeeptest=ndpost
    cat('*****nkeeptest set to ndpost\n')
  }
  if((nkeeptestmean!=0) & ((ndpost %% nkeeptestmean) != 0)) {
    nkeeptestmean=ndpost
    cat('*****nkeeptestmean set to ndpost\n')
  }
  if((nkeeptreedraws!=0) & ((ndpost %% nkeeptreedraws) != 0)) {
    nkeeptreedraws=ndpost
    cat('*****nkeeptreedraws set to ndpost\n')
  }
  #--------------------------------------------------
  #prior: change the dimension
  #--------------------------------------------------

  nu=sigdf
  if(is.na(lambda)) {
    lambda = rep(NA, K)
    if(is.na(sigest)) {
      sigest = rep(NA, K)
      for(k in 1:K){
        if(p < n[k]) {
          df = data.frame(t(x.train[[k]]),y.train[[k]])
          lmf = lm(y.train~.,df)
          sigest[k] = summary(lmf)$sigma
        } else {
          sigest[k] = sd(y.train[[k]])
        }
      }
    }
    qchi = qchisq(1.0-sigquant,nu)
    lambda = (sigest*sigest*qchi)/nu #lambda parameter for sigma prior
  } else {
    sigest=sqrt(lambda)
  }

  if(is.na(sigmaf)) {
    sigmaf = tau = rep(NA, K)
    for(k in 1:K){
      tau[k]=(max(y.train[[k]])-min(y.train[[k]]))/(2*bk*sqrt(ntree))
    }
  } else {
    tau = sigmaf/sqrt(ntree)
  }
  #--------------------------------------------------
  ptm <- proc.time()
  #call
  res = .Call("JointBart",
              n,  #number of observations in training data
              p,  #dimension of x
              np, #number of observations in test data
              x.train,   #pxn training data x
              y.train,   #pxn training data x
              x.test,   #p*np test data x
              ntree,
              numcut,
              ndpost*keepevery,
              nskip,
              power,
              base,
              tau,
              nu,
              lambda,
              sigest,
              w,
              sparse,
              theta,
              omega,
              grp,
              a,
              b,
              rho,
              augment,
              nkeeptrain,
              nkeeptest,
              nkeeptestmean,
              nkeeptreedraws,
              printevery,
              xinfo
  )

  res$proc.time <- proc.time()-ptm

  res$mu = fmean
  res$yhat.train.mean = res$yhat.train.mean+fmean
  res$yhat.train = res$yhat.train+fmean
  res$yhat.test.mean = res$yhat.test.mean+fmean
  res$yhat.test = res$yhat.test+fmean
  if(nkeeptreedraws>0)
    names(res$treedraws$cutpoints) = dimnames(x.train)[[1]]
  dimnames(res$varcount)[[2]] = as.list(dimnames(x.train)[[1]])
  dimnames(res$varprob)[[2]] = as.list(dimnames(x.train)[[1]])
  ##res$nkeeptreedraws=nkeeptreedraws
  res$varcount.mean <- apply(res$varcount, 2, mean)
  res$varprob.mean <- apply(res$varprob, 2, mean)
  res$rm.const <- rm.const
  attr(res, 'class') <- 'wbart'
  return(res)
}
