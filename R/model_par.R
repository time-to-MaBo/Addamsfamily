model.par <- function(theta,num.a,num.g,J,alpha.idx,gamma.idx,a.idx,PVF.idx,IG.idx,
                      aNBp.idx,aNB.idx,cP.idx,H.idx,aB.idx=FALSE,N=NA){

  if(num.a > 0){
    alpha <- theta[1:num.a]
    gamma <- theta[-(1:num.a)]
  }else{
    alpha <- NA
    gamma <- theta
  }

  model_par <- matrix(NA, nrow = J, ncol = 2)
  colnames(model_par) <- c('alpha', 'gamma')
  model_par[alpha.idx, 'alpha'] <- alpha
  model_par[gamma.idx, 'gamma'] <- gamma

  model_par[,'gamma'] <- exp(model_par[,'gamma'])
  model_par[a.idx,'alpha'] <- exp(model_par[a.idx, 'alpha'])
  model_par[PVF.idx,'alpha'] <- exp(model_par[PVF.idx, 'alpha']) - 1
  model_par[IG.idx, 'alpha'] <- -0.5

  model_par[aNB.idx, 'alpha'] <- exp(model_par[aNB.idx, 'alpha']) # neu
  model_par[aNB.idx, 'gamma'] <- model_par[aNB.idx, 'gamma'] + model_par[aNB.idx, 'alpha'] # neu
  model_par[aNBp.idx, 'alpha'] <- -exp(model_par[aNBp.idx, 'alpha'])
  #model_par[aNB.idx, 'alpha'] <- model_par[aNB.idx, 'gamma'] * exp(-exp(model_par[aNB.idx, 'alpha']))
  #model_par[aNB.idx, 'alpha'] <- model_par[aNB.idx, 'gamma'] * exp(-exp(model_par[aNB.idx, 'alpha']))
  model_par[aB.idx, 'alpha'] <- model_par[aB.idx, 'gamma'] +1/N ##add discrete
  model_par[cP.idx, 'alpha'] <- model_par[cP.idx, 'alpha']+1
  model_par[H.idx, 'alpha'] <- -exp(-(model_par[H.idx, 'alpha'] + 1))

  model_par[model_par == 0] <- .Machine$double.xmin

  return(model_par)

}
