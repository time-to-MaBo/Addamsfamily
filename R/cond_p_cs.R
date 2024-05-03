cond.p.cs <- function(par,X,d,z,z2,K,num.lam.p,dt_1,scale.cond,n_xx){

  beta <- par[1:K]
  psi_i <- as.vector(z*exp(X%*%beta))
  lambda <- exp(par[(K+1):(K+num.lam.p)])
  lambda[lambda == 0] <- .Machine$double.xmin
  Lambda <- colSums(lambda*dt_1)

  zL_i <- psi_i * Lambda
  approx_term <- -0.5*(z2-z^2)*(exp(X%*%beta)*Lambda)^2 * exp(-zL_i)/(1-exp(-zL_i))^2

  ll <- sum( -(1-d)*zL_i*n_xx + d*( log(1-exp(-zL_i))*n_xx + approx_term*n_xx )  )
  if(is.infinite(ll) | is.nan(ll)){
    ll <- -sqrt(.Machine$double.xmax)
    warning('conditional loglihood approached infinity. sqrt(.Machine$double.xmax) imputed.
            Consider setting -1 < scale.cond < 0 or to tackle this via control of optimax.')
  }

  return(scale.cond*ll)

}
