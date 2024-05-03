piecewise_haz <- function(dt,tr.dt,m,d,psi_i,idx.lam,numh,num.lam,breaks,strata,Qstar){

  par <- rep(NA, num.lam)
  lambda <- list()

  nominator <- t(t(m)*d)
  denominator <- t(psi_i*t(dt-tr.dt))

  for(q in 1:Qstar){
    nominator.temp <- rowSums(nominator[(idx.lam[q]+1):idx.lam[q+1],,drop=FALSE])
    denominator.temp <- rowSums(denominator[(idx.lam[q]+1):idx.lam[q+1],,drop=FALSE])

    par[(idx.lam[q]+1):idx.lam[q+1]] <- nominator.temp/denominator.temp
  }

  par[par == 0] <- .Machine$double.xmin

  for(q in 1:Qstar){
    lambda[[q]] <- par[(idx.lam[q]+1):(idx.lam[q+1])]

    lam_name <- rep(NA, numh[q])
    for(i in 1:numh[q]){
      lam_name[i] <- paste('(', breaks[[q]][i], ',', breaks[[q]][i+1] ,']', sep = '')
    }
    names(lambda[[q]]) <- lam_name
  }

  names(lambda) <- strata

  return(lambda)

}
