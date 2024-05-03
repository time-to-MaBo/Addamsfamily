partial.loglihood <- function(beta,z,y,d,X,K,dt,tr.dt,m,n,scale.cond = -1){
  #hier d?rfen keine GG-Elemente rein!!!

  R_i <- risk.set_i(X = X, beta = beta, z = z, dt = dt, m = m,
                    d = d, tr.dt = tr.dt, n = n)
  loglihood <- sum(  (X[d==1,,drop=FALSE]%*%beta - log(R_i)) )
  # truncation is considered via the risk set R_i

  return(loglihood*scale.cond)
}
