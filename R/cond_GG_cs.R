cond.GG.cs <- function(par,X,d,z,z2,K,num.lam,num.lam.p,y,Q,n,stratum,lam.idx,n.GG,partial.idx,
                       partial_i,GG_i,dt_1,GG.op,GG.sc,GG.special,scale.cond,n_xx){

  beta <- par[1:K]
  psi_i <- as.vector(z*exp(X%*%beta))
  lnlambda.p <- par[(K+1):(K+num.lam+num.lam.p)]
  lnlambda <- lnlambda.p[lam.idx]
  lnlambda.p[lam.idx] <- -Inf

  if(Q>0){
    tmp <- seq(from=1,to=3*Q,by=3)
  }else{
    tmp <- NULL
  }
  if(GG.special){
    lambda <- lnlambda[GG.sc]
    lambda[GG.op=='one'] <- 1
    lambda[GG.op=='negone'] <- -1
    lambda[GG.op=='zero'] <- 0
    lambda[GG.op=='exp'] <- exp(lambda[GG.op=='exp'])
    lambda[GG.op=='negexp'] <- -exp(lambda[GG.op=='negexp'])
    lambda[GG.op=='invexp'] <- exp(-lambda[GG.op=='invexp'])
    lambda[GG.op=='neginvexp'] <- -exp(-lambda[GG.op=='neginvexp'])
    lambda[GG.op=='sqrt2'] <- sqrt(2)
    lambda[GG.op=='sqrt2^-1'] <- 1/sqrt(2)
  }else{
    lambda <- lnlambda
    lambda[tmp] <- exp(lnlambda[tmp])
  }

  lambda[tmp][lambda[tmp]==0] <- .Machine$double.xmin
  sigma <- lambda[tmp]
  eta <- lambda[tmp+1]
  nu <- lambda[tmp+2]

  Lambda <- rep(NA, n)
  Lambda[GG_i] <- GG.haz(s = sigma, e = eta, nu = nu, y = y, str = stratum, n = n.GG)[,1]

  if(any(partial.idx)){
    lambda.p <- exp(lnlambda.p)
    lambda.p[lambda.p==0 & !lam.idx] <- .Machine$double.xmin
    Lambda[partial_i] <- colSums(dt_1*(lambda.p[!lam.idx]))
  }

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
