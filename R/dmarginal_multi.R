dmarginal.multi <- function(par,X, partial_i, GG_i,GG.idx,fac,frail.thresh,
                            K, num.lam, cum.numh.p,n,numh.p, num.lam.p, breaks, p.strata,
                            tr,Q.p,Q.GG,y,d,n.GG,n.p,stratum.GG,dt_1,tr.dt_1,m_1,X.tilde,X.tilde_C,
                            ID, scale.marginal,slope,slope.model,sum.d_i, d.marg, GG.strata, max.n_i,v,N,
                            partial.idx, num.a, num.g, J, n_C, tr.idx,gamma.rank, alpha.rank,
                            ag.idx, a.idx, g.idx, IG.idx, PVF.idx,aNBp.idx,aNB.idx,cP.idx,H.idx,
                            ag_idx, a_idx, g_idx, IG_idx, PVF_idx,z,idx, frailty, alpha.idx, gamma.idx,
                            ag.type, a.type, g.type,IG.type, PVF.type, GG.op, GG.sc, GG.special){

  z <- as.vector(z*exp(v*X.tilde))
  z.p <- z[partial_i]
  z.GG <- z[GG_i]
  if(K>0){
    beta <- par[(num.a+num.g+1):(num.a+num.g+K)]
  }else{
    beta <- 0
  }
  exb <- as.vector(exp(X%*%beta))
  psi_i <- as.vector(z*exb)

  Lambda <- matrix(NA, nrow = n, ncol = 2)
  tr.Lambda <- Lambda
  colnames(Lambda) <- c('Lambda', 'lambda')
  partial.lambda <- NULL # Pre-implementation for case of all(!partial.idx) = TRUE

  if(any(GG.idx)){
    lnlambda <- par[(num.a+num.g+K+1):(num.a+num.g+K+num.lam)]
    tmp <- seq(from=1,to=3*Q.GG,by=3)
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
    names(sigma) <- GG.strata
    names(eta) <- GG.strata
    names(nu) <- GG.strata

    Lambda[GG_i,] <- GG.haz(s = sigma, e = eta, nu = nu, n = n.GG,
                            str = stratum.GG, y = y[GG_i])
    tr.Lambda[GG_i,1] <- GG.haz(s = sigma, e = eta, nu = nu, n = n.GG,
                                str = stratum.GG, y = tr[GG_i])[,1]
  }

  if(any(partial.idx)){
    partial.lambda <- piecewise_haz(dt=dt_1, tr.dt=tr.dt_1, m=m_1,
                                    d=d[partial_i],psi_i=psi_i[partial_i],
                                    idx.lam=cum.numh.p,numh=numh.p,
                                    num.lam=num.lam.p,breaks=breaks,strata=p.strata,
                                    Qstar=Q.p)

    Lambda[partial_i,1] <- colSums(dt_1*unlist(partial.lambda))
    Lambda[partial_i,2] <- (m_1*unlist(partial.lambda))[m_1==1]
    tr.Lambda[partial_i,1] <- colSums(tr.dt_1*unlist(partial.lambda))
  }
  Lambda_i <- Lambda*exb
  tr.Lambda_i <- tr.Lambda*exb

  # marginal model
  s <- as.vector(tapply(X = Lambda_i[,1], INDEX = ID, FUN = sum, simplify = TRUE))
  tr.s <- as.vector(tapply(X = tr.Lambda_i[,1], INDEX = ID, FUN = sum, simplify = TRUE))
  lambda_i <- tapply(X = Lambda_i[,2], INDEX = ID, FUN = function(x) c(x, rep(0, max.n_i-length(x))), simplify = TRUE)
  lambda_i <- matrix(unlist(lambda_i), ncol = max.n_i, byrow = TRUE)
  lambda_i <- lambda_i*d.marg
  lambda_i[lambda_i == 0] <- 1
  lnlambda_i <- rowSums(log(lambda_i))

  theta <- par[1:(num.a+num.g)]

  if(!slope){
    ll <- marginal.multi(par = theta,
                         s = s, tr.s = tr.s, num.a = num.a, num.g = num.g,
                         sum.d_i = sum.d_i, tr.idx = tr.idx, n_C = n_C,
                         ag.idx = ag.idx, a.idx = a.idx, g.idx = g.idx,
                         IG.idx = IG.idx, PVF.idx = PVF.idx,
                         ag_i = ag_idx, a_i = a_idx, g_i = g_idx,
                         IG_i = IG_idx, PVF_i = PVF_idx, J = J,fac = fac,
                         aNBp.idx=aNBp.idx,aNB.idx=aNB.idx,cP.idx=cP.idx,H.idx=H.idx,
                         idx = as.numeric(frailty), alpha.idx = alpha.idx, gamma.idx = gamma.idx,
                         rank_a = alpha.rank, rank_g = gamma.rank,
                         type_ag = ag.type, type_a = a.type, type_g = g.type,frail.thresh=frail.thresh,
                         type_IG = IG.type, type_PVF = PVF.type, scale.marginal = scale.marginal)
  }else{
    Lambda_i <- tapply(X = Lambda_i[,1], INDEX = ID, FUN = function(x) c(x, rep(0, max.n_i-length(x))), simplify = TRUE)
    Lambda_i <- matrix(unlist(Lambda_i), ncol = max.n_i, byrow = TRUE)
    tr.Lambda_i <- tapply(X = tr.Lambda_i[,1], INDEX = ID, FUN = function(x) c(x, rep(0, max.n_i-length(x))), simplify = TRUE)
    tr.Lambda_i <- matrix(unlist(tr.Lambda_i), ncol = max.n_i, byrow = TRUE)
    Lambda_i2 <- tr.Lambda_i2 <- c()
    for(i in 1:(N+1)){
      Lambda_i2 <- rbind(Lambda_i2,Lambda_i)
      tr.Lambda_i2 <- rbind(tr.Lambda_i2,tr.Lambda_i)
    }
    Lambda_i <- Lambda_i2
    tr.Lambda_i <- tr.Lambda_i2
    rm(list = c('Lambda_i2','tr.Lambda_i2'))
    rownames(Lambda_i) <- rep(unique(ID),N+1)
    rownames(tr.Lambda_i) <- rep(unique(ID),N+1)
    lnlambda_i <- matrix(rep(lnlambda_i,N+1),ncol=N+1)
    rownames(lnlambda_i) <- unique(ID)

    lnpi <- par[num.a+num.g+K+num.lam+1]
    if(slope.model == 'free-Willy'){
      lnphipsi <- par[(num.a+num.g+K+num.lam+2):(num.a+num.g+K+num.lam+N+1)]
    }else{
      lnphipsi <- par[(num.a+num.g+K+num.lam+2):(num.a+num.g+K+num.lam+3)]
    }

    ll <- slope.opt(par = c(theta,lnpi,lnphipsi),Lambda_i=Lambda_i,
                    tr.Lambda_i=tr.Lambda_i, n_C=n_C,X.tilde_C=X.tilde_C,d.marg=d.marg, num.a=num.a, num.g=num.g, N=N,
                    J=J, alpha.idx = alpha.idx,gamma.idx = gamma.idx, a.idx = a.idx,tr.idx=tr.idx,
                    PVF.idx = PVF.idx, IG.idx = IG.idx, aNBp.idx=aNBp.idx,aNB.idx=aNB.idx,
                    cP.idx=cP.idx,H.idx=H.idx,sum.d_i=sum.d_i,idx=idx,fac = fac,
                    ag_i=ag_idx,g_i=g_idx,a_i=a_idx,PVF_i=PVF_idx,IG_i=IG_idx,ag.idx=ag.idx,g.idx=g.idx,slope.model=slope.model,
                    scale.marginal=scale.marginal,frail.thresh=frail.thresh)
  }

  return(ll+sum(lnlambda_i))
}
