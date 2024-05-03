GG.loglihood <- function(beta,lnlambda,X,z,y,tr,d,GG.idx,Qstar,n,stratum,GG.op,GG.sc,GG.special){

  tmp <- seq(from=1,to=3*Qstar,by=3)
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

  Lambda <- GG.haz(s = sigma, e = eta, nu = nu, n = n, y = y,
                   str = stratum)
  tr.Lambda <- GG.haz(s = sigma, e = eta, nu = nu, n = n, y = tr,
                      str = stratum)[,1]

  xb <- X%*%beta
  psi_i <- as.vector(z*exp(xb))
  #lambda_i <- psi_i*Lambda[,2] # formerly took the role of xb in loglihood
  s <- psi_i*Lambda[,1]
  tr.s <- psi_i*tr.Lambda
  #lambda_i[lambda_i==0 & d == 0] <- 1 # z does not need to be included! (constant)

  loglihood <- sum(d*(xb+log(Lambda[,2]))-s+tr.s)
  if(is.infinite(loglihood)){
    loglihood <- -sqrt(.Machine$double.xmax)
    warning('conditional loglihood approached infinity. -sqrt(.Machine$double.xmax) imputed.
            Consider setting -1 < scale.cond < 0 or to tackle this via control of optimax.')
  }

  return(loglihood)
}
