dcs.multi <- function(par,K,num.lam,num.lam.p,X,y,n,n_C,n.GG,
                      partial.idx,GG.idx,n_xx,sign.info,overdisp,n.date,
                      partial_i,ag_i,a_i,g_i,PVF_i,ag.idx,a.idx,g.idx,PVF.idx,Q,
                      IG.idx,tr.idx,num.a,num.g,J,alpha.idx,gamma.idx, dt_1, stratum,
                      idx,sum.d_i,max.n_i,GG_i,lam.idx,
                      scale.marginal, ID, d.marg, frail.thresh,
                      GG.op,GG.sc,GG.special,sc.inv,aNBp.idx,aNB.idx,aB.idx=FALSE,N=NA,cP.idx,H.idx,
                      num.p,p.idx## add discrete ----
){
  # caution Q = Q.GG, stratum = GG.stratum, y = y.GG!
  theta <- par[1:(num.a+num.g)]
  p <- rep(0,J) ## add discrete ----
  if(num.p>0) p[p.idx] <- exp(par[(num.a+num.g+1):(num.a+num.g+num.p)]) ## add discrete ----

  beta <- par[(num.a+num.g+num.p+1):(num.a+num.g+num.p+K)]## add discrete ----
  exb <- as.vector(exp(X%*%beta))
  exb[exb==0] <- .Machine$double.xmin # it sometime occured that in early iterations hazard parameter approached infinity: 0*Inf=NaN

  lnlambda.p <- par[(num.a+num.g+num.p+K+1):(num.a+num.g+num.p+K+num.lam+num.lam.p)]
  lnlambda <- lnlambda.p[lam.idx]
  lnlambda.p[lam.idx] <- -Inf
  Lambda <- rep(NA, n)

  if(any(GG.idx)){
    tmp <- seq(from=1,to=3*Q,by=3) # caution Q is Q.GG
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

    Lambda[GG_i] <- GG.haz(s = sigma, e = eta, nu = nu, y = y, str = stratum, n = n.GG)[,1]
  }

  lambda.p <- exp(lnlambda.p)
  if(any(partial.idx)){
    lambda.p[lambda.p==0 & !lam.idx] <- .Machine$double.xmin
    tmp <- dt_1*(lambda.p[!lam.idx])
    tmp[is.nan(tmp)] <- 0 # it sometimes occured that in early iterations hazard parameter approached infinity: 0*Inf=NaN
    Lambda[partial_i] <- colSums(tmp)
  }

  Lambda_i <- Lambda*exb
  Lambda_i <- tapply(X = Lambda_i, INDEX = ID, FUN = function(x) c(x,rep(0,max.n_i-length(x))),
                     simplify = FALSE)
  Lambda_i <- t(sapply(X = Lambda_i, FUN = function(x) unlist(x)))
  Lambda_i[is.infinite(Lambda_i)] <- .Machine$double.xmax

  s <- matrix( 0, nrow = n_C, ncol = 1+sum( choose(max(sum.d_i), 1:max(sum.d_i)) ) )
  s[,1] <- rowSums(Lambda_i*(1-d.marg))
  sapply(X = 1:n_C, FUN = function(i){
    d_i <- sum.d_i[i]
    if(d_i>0){
      tmp.Lambda <- Lambda_i[i,d.marg[i,] == 1]
      sapply(X = 1:d_i, FUN = function(j){
        tmp <- as.matrix(combn(x = d_i, m = j))
        tmp2 <- apply(X = tmp, MARGIN = 2, FUN = function(x) sum(tmp.Lambda[x]))
        s[i,(2+sum(choose(d_i, 1:(j-1))*(j>1)) ) : (1+sum(choose(d_i, 1:(j))))] <<- tmp2
      })
    }
  })
  # for(i in 1:n_C){
  #   d_i <- sum.d_i[i]
  #   if(d_i>0){
  #     tmp.Lambda <- Lambda_i[i,d.marg[i,] == 1]
  #     for(j in 1:d_i){
  #       tmp <- as.matrix(utils::combn(x = d_i, m = j))
  #       tmp2 <- apply(X = tmp, MARGIN = 2, FUN = function(x) sum(tmp.Lambda[x]))
  #       s[i,(2+sum(choose(d_i, 1:(j-1))*(j>1)) ) : (1+sum(choose(d_i, 1:(j))))] <- tmp2
  #     }
  #   }
  # }
  if(any(is.infinite(s))){
    warning('cumulative hazard rate approached infinity.')
    s[is.infinite(s)] <- sqrt(.Machine$double.xmax)
  }

  if(overdisp){
    tau <- exp(par[(num.a+num.g+num.p+K+num.lam+num.lam.p+1):length(par)])
  }else{
    tau <- NULL
  }

  ll <- marginal.cs(par = theta, tau=tau, s = s, num.a = num.a, num.g = num.g,
                    n = n_C,aNBp.idx=aNBp.idx,aNB.idx=aNB.idx,aB.idx=aB.idx,N=N,
                    cP.idx=cP.idx,H.idx=H.idx,n_xx=n_xx,overdisp=overdisp,
                    ag.idx = ag.idx, a.idx = a.idx, g.idx = g.idx,
                    IG.idx = IG.idx, PVF.idx = PVF.idx, n.date=n.date,
                    ag_i = ag_i, a_i = a_i, g_i = g_i, frail.thresh=frail.thresh,
                    PVF_i = PVF_i, J = J, sign.info = sign.info,
                    idx = idx, alpha.idx = alpha.idx, gamma.idx = gamma.idx,
                    scale.marginal = scale.marginal,
                    p=p## add discrete ----
  )

  return(ll)
}
