GG.multi <- function(theta, lnpi, lnphipsi, n, z, v, X.tilde,X.tilde_C, X, N,
                     partial_i, GG_i, beta, lnlambda,GG.op,GG.sc,GG.special,
                     K, num.lam, cum.numh.p,numh.p, num.lam.p,breaks,p.strata, d.marg, n_i, max.n_i,
                     scale.cond, scale.marginal,ID,tr.idx,Q.p,Q.GG,y,tr,d,n.GG,n.p,stratum.GG,dt_1,tr.dt_1,m_1,
                     max.cond = TRUE, iter = 0, marginal = list(value = Inf), GG.strata, partial.idx,
                     num.a, num.g, J, cond = marginal, gamma.rank, alpha.rank, n_C, sum.d_i,Addams.idx,
                     ag.idx, a.idx, g.idx, IG.idx, PVF.idx, ag_idx, a_idx, g_idx, IG_idx, PVF_idx,
                     aNBp.idx,aNB.idx,cP.idx,H.idx, slope,slope.model, fac = fac,frail.thresh,
                     idx, frailty, alpha.idx, gamma.idx, ag.type, a.type, g.type, IG.type, PVF.type,
                     c.options,m.options,R.options, converge, prnt, iterlim,initial.z,start_dist){

  c.method <- c.options$method
  m.method <- m.options$method
  Richardson <- m.options$Richardson
  if(Richardson){
    c.gr <- function(par,...){
      numDeriv::grad(func = cond.loglihood, x = par, method = 'Richardson',
           side = NULL, method.args = R.options, ...)
    }
    c.Hessian <- function(par,...){
      numDeriv::hessian(func = cond.loglihood, x = par,
              method = 'Richardson', method.args = R.options, ...)
    }
    func <- ifelse(slope,slope.opt,marginal.multi)
    m.gr <- function(par,...){
      numDeriv::grad(func = func, x = par, method = 'Richardson',
           side = NULL, method.args = R.options, ...)
    }
    m.Hessian <- function(par,...){
      numDeriv::hessian(func = func, x = par,
              method = 'Richardson', method.args = R.options, ...)
    }
  }else{
    m.gr <- c.gr <- m.Hessian <- c.Hessian <- NULL
  }

  pi <- NULL
  v.Domain <- NULL

  while(max.cond){

    iter <- iter + 1
    old.marginal <- marginal
    old.cond <- cond

    z <- z*exp(v*X.tilde)
    z.p <- z[partial_i]
    z.GG <- z[GG_i]

    old.marginal <- marginal
    old.cond <- cond

    cond <- optimax(par = c(beta,lnlambda), fn = cond.loglihood, gr = c.gr,
                    Hessian = c.Hessian, options = c.options, Richardson = Richardson,
                    R.options = R.options, method = c.method,
                    X.GG = X[GG_i,,drop=FALSE], y.GG = y[GG_i], tr.GG = tr[GG_i], z.GG = z.GG,
                    d.GG = d[GG_i],n.GG = n.GG, X.p=X[partial_i,,drop=FALSE], z.p=z.p, d.p=d[partial_i],
                    n.p=n.p, Q.GG = Q.GG, stratum.GG = stratum.GG, dt = dt_1, m = m_1, tr.dt = tr.dt_1,
                    num.lam=num.lam,K=K,GG.op=GG.op,
                    GG.sc=GG.sc,GG.special=GG.special,scale.cond = scale.cond)

    beta <- cond$par[1:K]
    exb <- as.vector(exp(X%*%beta))
    psi_i <- as.vector(z*exb)
    Lambda <- matrix(NA, nrow = n, ncol = 2)
    tr.Lambda <- Lambda
    colnames(Lambda) <- c('Lambda', 'lambda')
    partial.lambda <- NULL # Pre-implementation for case of all(!partial.idx) = TRUE

    lnlambda <- cond$par[(K+1):(K+num.lam)]
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

    if(!slope){
      marginal <- optimax(par = theta, fn = marginal.multi, gr = m.gr, Hessian = m.Hessian,
                          options = m.options, Richardson = Richardson, R.options = R.options,
                          method = m.method,
                          s = s, tr.s = tr.s, num.a = num.a, num.g = num.g,
                          sum.d_i = sum.d_i, tr.idx = tr.idx, n_C = n_C,
                          ag.idx = ag.idx, a.idx = a.idx, g.idx = g.idx,
                          IG.idx = IG.idx, PVF.idx = PVF.idx,
                          ag_i = ag_idx, a_i = a_idx, g_i = g_idx,
                          IG_i = IG_idx, PVF_i = PVF_idx, J = J,
                          idx = as.numeric(frailty), alpha.idx = alpha.idx, gamma.idx = gamma.idx,
                          aNBp.idx=aNBp.idx,aNB.idx=aNB.idx,cP.idx=cP.idx,H.idx=H.idx,
                          rank_a = alpha.rank, rank_g = gamma.rank, fac = fac,
                          type_ag = ag.type, type_a = a.type, type_g = g.type,
                          type_IG = IG.type, type_PVF = PVF.type, scale.marginal = scale.marginal,
                          frail.thresh=frail.thresh)
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

      marginal <- optimax(par = c(theta,lnpi,lnphipsi), fn = slope.opt, gr = m.gr, Hessian = m.Hessian,
                          options = m.options, Richardson = Richardson, R.options = R.options, method = m.method,
                          Lambda_i=Lambda_i, tr.Lambda_i=tr.Lambda_i, n_C=n_C,X.tilde_C=X.tilde_C,d.marg=d.marg, num.a=num.a, num.g=num.g, N=N,
                          J=J, alpha.idx = alpha.idx,gamma.idx = gamma.idx, a.idx = a.idx,tr.idx=tr.idx,
                          PVF.idx = PVF.idx, IG.idx = IG.idx, aNBp.idx=aNBp.idx,aNB.idx=aNB.idx,
                          cP.idx=cP.idx,H.idx=H.idx,sum.d_i=sum.d_i,idx=idx, fac = fac,
                          ag_i=ag_idx,g_i=g_idx,a_i=a_idx,PVF_i=PVF_idx,IG_i=IG_idx,ag.idx=ag.idx,g.idx=g.idx,slope.model=slope.model,
                          frail.thresh=frail.thresh,scale.marginal=scale.marginal)

      lnpi <- marginal$par[num.a+num.g+1]
      pi <- exp(-exp(lnpi))
      if(slope.model == 'free-Willy'){
        lnphipsi <- marginal$par[(num.a+num.g+2):(num.a+num.g+N+1)]
        phi <- cumsum(-exp(lnphipsi[1:(N/2)]))[(N/2):1]
        psi <- cumsum(exp(lnphipsi[(N/2+1):N]))
      }else{
        lnphipsi <- marginal$par[(num.a+num.g+2):(num.a+num.g+3)]
        phi <- -exp(lnphipsi[1])*((N/2):1)
        psi <- exp(lnphipsi[2])*(1:(N/2))
      }
      v.Domain <- c(phi,0,psi)
    }

    marginal$value.nolam <- marginal$value
    marginal$value <- marginal$value + scale.marginal*sum(lnlambda_i)

    theta <- marginal$par[1:(num.a+num.g)]
    model_par <- model.par(theta, num.a, num.g, J, alpha.idx, gamma.idx, a.idx, PVF.idx, IG.idx,
                           aNBp.idx,aNB.idx,cP.idx,H.idx)
    if(iter == 1){
      uni.logli <- cond$value[1]/scale.cond
      if(initial.z){
        # This condition implements z with strating values of parameters if
        ## estimated frailty parameters are close to a no-frailty model after
        ## first iteration optimization
        model_par1 <- model.par(start_dist, num.a, num.g, J, alpha.idx, gamma.idx, a.idx, PVF.idx, IG.idx,
                                aNBp.idx,aNB.idx,cP.idx,H.idx)

        if(any(model_par[,'alpha'] > model_par[,'gamma'] & ag.idx, na.rm = TRUE)){
          tmp <- (model_par[,'alpha'] > model_par[,'gamma']) & ag.idx
          t.idxa <- alpha.rank[tmp]
          t.idxg <- gamma.rank[tmp]
          model_par[tmp,] <- model_par1[tmp,]
          theta[1:num.a][t.idxa] <- start_dist[1:num.a][t.idxa]
          theta[(num.a+1):(num.a+num.g)][t.idxg] <- start_dist[(num.a+1):(num.a+num.g)][t.idxg]
          warning('In 1st global iteration frailty stratum ', which(tmp), ' has invalid parameter constellation and
                  was set to starting values again (in order to compute E[Z|data]).')
        }else if(any(model_par[,'gamma'] < frail.thresh*1000 & (ag.idx | g.idx), na.rm = TRUE)){
          tmp <- model_par[,'gamma'] < frail.thresh*1000 & (ag.idx | g.idx)
          t.idxg <- gamma.rank[tmp]
          model_par[tmp,] <- model_par1[tmp,]
          theta[(num.a+1):(num.a+num.g)][t.idxg] <- start_dist[(num.a+1):(num.a+num.g)][t.idxg]
          if(num.a>0){
            theta[1:num.a][ag.idx & tmp] <- start_dist[1:num.a][ag.idx & tmp]
          }
          warning('In 1st global iteration gamma of stratum ', which(tmp), ' is close to zero and
                  was set to starting values again (in order to compute E[Z|data]).')
        }else if(any(model_par[,'alpha'] < frail.thresh*1000 & a.idx, na.rm = TRUE)){
          tmp <- model_par[,'alpha'] < frail.thresh*1000 & a.idx
          t.idxa <- alpha.rank[tmp]
          model_par[tmp,] <- model_par1[tmp,]
          theta[1:num.a][t.idxa] <- start_dist[1:num.a][t.idxa]
          warning('In 1st global iteration alpha of stratum ', which(tmp), ' is close to zero and
                  was set to starting values again (in order to compute E[Z|data]).')
        }else if(any(model_par[,'gamma']^-1 < frail.thresh*1000 & (PVF.idx | IG.idx), na.rm = TRUE)){
          tmp <- model_par[,'gamma']^-1 < frail.thresh*1000 & (PVF.idx | IG.idx)
          t.idxa <- alpha.rank[tmp]
          t.idxg <- gamma.rank[tmp]
          model_par[tmp,] <- model_par1[tmp,]
          theta[1:num.a][t.idxa] <- start_dist[1:num.a][t.idxa]
          theta[(num.a+1):(num.a+num.g)][t.idxg] <- start_dist[(num.a+1):(num.a+num.g)][t.idxg]
          warning('In 1st global iteration gamma^-1 of stratum ', which(tmp), ' is close to zero and
                  was set to starting values again (in order to compute E[Z|data]).')
        }else if(any(model_par[,'alpha']+1 < frail.thresh*1000 & PVF.idx, na.rm = TRUE)){
          tmp <- model_par[,'alpha']+1 < frail.thresh*1000 & PVF.idx
          t.idxa <- alpha.rank[tmp]
          t.idxg <- gamma.rank[tmp]
          model_par[tmp,] <- model_par1[tmp,]
          theta[1:num.a][t.idxa] <- tstart_dist[1:num.a][t.idxa]
          theta[(num.a+1):(num.a+num.g)][t.idxg] <- start_dist[(num.a+1):(num.a+num.g)][t.idxg]
          warning('In 1st global iteration alpha of stratum ', which(tmp), ' is close to zero and
                  was set to starting values again (in order to compute E[Z|data]).')
        }
      }
    }
    alpha <- model_par[,'alpha']
    gamma <- model_par[,'gamma']

    # Expectation step: Z
    if(any(Addams.idx)){
      g.idx[Addams.idx] <- abs(model_par[Addams.idx,'alpha'])< 1e-4
      a.idx[Addams.idx] <- abs(model_par[Addams.idx,'alpha']-model_par[Addams.idx,'gamma'])<1e-4 & !g.idx[Addams.idx]
      ag.idx[Addams.idx] <- !a.idx & !g.idx & Addams.idx
      ag_idx <- ag.idx[frailty]
      a_idx <- a.idx[frailty]
      g_idx <- g.idx[frailty]
    }

    if(!slope){
      z <- z.multi(alpha = alpha, gamma = gamma, ag_i = ag_idx, g_i = g_idx, a_i = a_idx, PVF_i = PVF_idx,
                   s = s, frailty = frailty, n = n_C, ID = ID, sum.d_i = sum.d_i, fac = fac, iter = iter)$z
    }else{
      z <- vz.frail(alpha=alpha, gamma=gamma, ag_i=ag_idx, g_i=g_idx, a_i=a_idx, PVF_i=PVF_idx,
                    d.marg=d.marg, frailty = frailty, n_C = n_C, ID = ID, sum.d_i = sum.d_i, iter = iter,
                    pi=pi, v.Domain=v.Domain, N=N, X.tilde_C=X.tilde_C, fac = fac, Lambda_i = Lambda_i)
      v <- z$v
      z <- z$z
    }

    # convergence condition met?
    diff.marginal <- (marginal$value - old.marginal$value)/scale.marginal
    marginal.worsened <- diff.marginal < 0
    diff.cond <- (cond$value - old.cond$value)/scale.cond

    max.cond <- diff.marginal >= converge | marginal.worsened
    if(prnt) {
      print('model par = ')
      print(model_par)
      print('GG par = ')
      print(lambda)
      print(paste('loglihood = ', marginal$value))
    }
    if(iter == iterlim & max.cond){iter <- iter+1; break}
  }


  return(list(marginal = marginal, cond = cond, beta = beta, iter = iter, uni.logli=uni.logli,
              diff.cond = diff.cond, diff.marginal = diff.marginal, z = z, v=v,v.Domain=v.Domain,
              lambda_i = lambda_i, theta = theta, model_par = model_par, lnpi=lnpi, pi=pi,lnphipsi=lnphipsi,
              alpha = alpha, gamma = gamma,  lambda = lambda, lnlambda = lnlambda,
              partial.lambda = partial.lambda, s = s, tr.s = tr.s,Lambda_i=Lambda_i,
              sigma = sigma, eta = eta, nu = nu))
}
