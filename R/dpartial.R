dpartial.db <- function(beta,z,y,d,X,K,dt,tr.dt,m,n,scale.cond = -1){
  # hier d?rfen keine GG-Elemente rein!

  psi_i <- as.vector(z * exp(X%*%beta))
  R <- t(psi_i*t(dt-tr.dt))

  dR_i <- t(apply(X = m[,d==1,drop=FALSE], MARGIN=2, FUN = function(x) t(x)%*%R))
  dR_i <- t(apply(X = dR_i, MARGIN = 1, FUN = function(r){
    colSums(apply(X = X, MARGIN = 2,
                  FUN = function(x) x*r ))
  }))
  if(K == 1) dR_i <- t(dR_i)
  rownames(dR_i) <- paste('dR(', y[d==1],')',sep='')
  colnames(dR_i) <- paste('/db',1:K,sep='')
  # in den Spalten sind die verschiedenen x,
  ## die Zeilen beziehen sich auf die versch. unzensierten Personen
  R_i <- risk.set_i(z = z, X = X, beta = beta, dt = dt, tr.dt = tr.dt,
                    m = m, d = d, n = n)

  dl.db <-  colSums( X[d==1,,drop=FALSE] - dR_i/R_i )

  return(dl.db*scale.cond)
}
