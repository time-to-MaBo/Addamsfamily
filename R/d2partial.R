d2partial.db2 <- function(z,X,K,beta,d,y,tr,dt,tr.dt,m,n,scale.cond){
  # hier d?rfen keine GG-Elemente rein!

  R_i <- risk.set_i(z = z, X = X, beta = beta, dt = dt, tr.dt = tr.dt,
                    m = m, d = d, n = n)

  psi_i <- as.vector(z * exp(X%*%beta))
  R <- t(psi_i*t(dt-tr.dt))

  dR_i1 <- t(apply(X = m[,d==1,drop=FALSE], MARGIN = 2, FUN = function(x) t(x)%*%R))
  dR_i <- t(apply(X = dR_i1, MARGIN = 1, FUN = function(r){
    colSums(apply(X = X, MARGIN = 2,
                  FUN = function(x) x*r ))
  }))
  if(K == 1) dR_i <- t(dR_i)
  rownames(dR_i) <- paste('dR(', y[d==1],')/',sep='')
  colnames(dR_i) <- paste('db',1:K,sep='')
  R1 <- dR_i/R_i
  R1 <- apply(X = R1, MARGIN = 1, FUN = list)
  R1 <- lapply(X = R1, FUN = unlist)
  R1 <- lapply(X = R1, FUN = function(x) x%*%t(x))

  dR <- replicate(list(), n=K*(1+K)/2)
  idx <- 1
  for(k in 1:K){
    tmp <-  t(t(dR_i1)*X[,k])
    for(k2 in k:K){
      tmp2 <- t(t(tmp)*X[,k2])
      colnames(tmp2) <- paste('j=',1:n, sep = '')
      rownames(tmp2) <- paste('R_{i,j}(t_i=',y[d==1],')', sep = '')
      dR[[idx + (k2-k)]] <- tmp2
    }
    names(dR)[idx:(idx+K-(k))] <- paste('dR(t)/db_', k, 'db_', k:K, sep = '')
    idx <- idx  + K - (k-1)
  }
  # columns are all individuals, rows are uncensored individuals
  # list elemens are derivatives od individualcontributions to risk set (see names).

  d2R_i <- sapply(X = dR, FUN = rowSums, simplify = TRUE)

  nstar <- sum(d)
  R2 <- replicate(list(matrix(NA, nrow = K, ncol = K)), n = nstar)
  idx <- lower.tri(R2[[1]])
  diag(idx) <- TRUE
  for(i in 1:nstar){
    R2[[i]][idx] <- d2R_i[i,] / R_i[i]
    R2[[i]][upper.tri(R2[[i]])] <- R2[[i]][lower.tri(R2[[i]])]
    rownames(R2[[i]]) <- paste('db_', 1:K , sep = '')
    colnames(R2[[i]]) <- paste('db_', 1:K , sep = '')
  }
  names(R2) <- paste('d2R(',y[d==1],')/',sep='')

  d2l.db2 <- 0
  for(i in 1:nstar){
    d2l.db2 <- d2l.db2 + R1[[i]] - R2[[i]]
  }

  return(d2l.db2*scale.cond)
}
