dcs.dblam <- function(par,X,d,K,z=NA,z2=NA,num.lam.p,dt_1,scale.cond,n_xx){
  # for univariate data only !!!!
  beta <- par[1:K]
  psi_i <- as.vector(exp(X%*%beta))
  lambda <- exp(par[(K+1):(K+num.lam.p)])
  lambda[lambda == 0] <- .Machine$double.xmin

  zL_i <- psi_i * colSums(lambda*dt_1)
  lambdadt <- t(psi_i*t(lambda*dt_1))
  db <- dlam <- 0
  for(i in 1:length(d)){
    db <-  -zL_i[i] * X[i,] * (1 - d[i]/(1-exp(-zL_i[i])))*n_xx[i] + db
    dlam <- -lambdadt[,i] * (1 - d[i]/(1-exp(-zL_i[i])))*n_xx[i] + dlam
  }

  diff <- c(db,dlam)
  return(scale.cond*diff)
}
