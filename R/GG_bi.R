GG.bi <- function(theta, n, z, X, partial_i, GG_i, beta, lnlambda,
                  K, num.lam, cum.numh.p,numh.p, num.lam.p,breaks,p.strata, d.marg, d_1, d_2, y, tr,
                  scale.cond, scale.marginal, ID, tr.idx, n.p, n.GG, d, Q.p, Q.GG, stratum.GG, dt_1, tr.dt_1, m_1,
                  max.cond = TRUE, iter = 0, marginal = list(value = Inf), GG.strata, partial.idx,
                  num.a, num.g, J, cond = marginal, gamma.rank, alpha.rank, n_C,Addams.idx,
                  aNBp.idx,aNB.idx,cP.idx,H.idx,frail.thresh,initial.z,start_dist,
                  ag.idx, a.idx, g.idx, IG.idx, PVF.idx, ag_idx, a_idx, g_idx, IG_idx, PVF_idx,
                  idx, frailty, alpha.idx, gamma.idx, ag.type, a.type, g.type, IG.type, PVF.type,
                  GG.op,GG.sc,GG.special,c.options,m.options,R.options, converge, prnt, iterlim){

  Richardson <- m.options$Richardson
  m.method <- m.options$method
  c.method <- c.options$method
  if(Richardson){
    c.gr <- function(par,...){
      numDeriv::grad(func = cond.loglihood, x = par, side = NULL,
           method = 'Richardson', method.args = R.options, ...)
    }
    c.Hessian <- function(par,...){
      numDeriv::hessian(func = cond.loglihood, x = par,
              method = 'Richardson', method.args = R.options, ...)
    }
    m.Hessian <- function(par,...){
      numDeriv::jacobian(func = dmarginal.dtheta, x = par, method = 'Richardson',
               side = NULL, method.args = R.options, ...)
    }
  }else{
    c.gr <- c.Hessian <- m.Hessian <- NULL
  }


  while(max.cond){

    iter <- iter + 1
    old.marginal <- marginal
    old.cond <- cond

    z.p <- z[partial_i]
    z.GG <- z[GG_i]

    old.marginal <- marginal
    old.cond <- cond

    cond <- optimax(par = c(beta,lnlambda), fn = cond.loglihood, gr = c.gr, Hessian = c.Hessian,
                    options = c.options, method = c.method, Richardson = Richardson, R.options = R.options,
                    num.lam=num.lam, K=K,
                    X.GG = X[GG_i,,drop=FALSE], y.GG=y[GG_i], tr.GG = tr[GG_i], z.GG=z.GG,d.GG=d[GG_i],
                    n.GG=n.GG,X.p=X[partial_i,drop=FALSE], z.p=z.p,d.p=d[partial_i],n.p=n.p,Q.GG=Q.GG,
                    stratum.GG=stratum.GG,dt=dt_1,m=m_1,tr.dt=tr.dt_1,
                    scale.cond = scale.cond,GG.op=GG.op,GG.sc=GG.sc,GG.special=GG.special)

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
      Lambda[partial_i,2] <- m_1*unlist(partial.lambda)
      tr.Lambda[partial_i,1] <- colSums(tr.dt_1*unlist(partial.lambda))
    }
    Lambda_i <- Lambda*exb
    tr.Lambda_i <- tr.Lambda*exb

    # marginal model
    s <- as.vector(tapply(X = Lambda_i[,1], INDEX = ID, FUN = sum, simplify = TRUE))
    tr.s <- as.vector(tapply(X = tr.Lambda_i[,1], INDEX = ID, FUN = sum, simplify = TRUE))
    lambda_i <- tapply(X = Lambda_i[,2], INDEX = ID, FUN = function(x) x, simplify = FALSE)
    lambda_i <- lapply(X = lambda_i, FUN = function(x) if(length(x)==2){x}else{c(x,0)})
    lambda_i <- matrix(unlist(lambda_i), ncol = 2, byrow = TRUE)

    marginal <- optimax(par = theta, fn = marginal.bi, gr = dmarginal.dtheta, Hessian = m.Hessian,
                        options = m.options, method = m.method, Richardson = Richardson, R.options = R.options,
                        s = s, tr.s = tr.s,  num.a = num.a, num.g = num.g,
                        d_1 = d_1, d_2 = d_2, tr.idx = tr.idx, n = n_C,
                        ag.idx = ag.idx, a.idx = a.idx, g.idx = g.idx,
                        IG.idx = IG.idx, PVF.idx = PVF.idx,
                        ag_i = ag_idx, a_i = a_idx, g_i = g_idx,
                        IG_i = IG_idx, PVF_i = PVF_idx, J = J,
                        aNBp.idx=aNBp.idx,aNB.idx=aNB.idx,cP.idx=cP.idx,H.idx=H.idx,
                        idx = as.numeric(frailty), alpha.idx = alpha.idx, gamma.idx = gamma.idx,
                        rank_a = alpha.rank, rank_g = gamma.rank,
                        type_ag = ag.type, type_a = a.type, type_g = g.type,
                        type_IG = IG.type, type_PVF = PVF.type, scale.marginal = scale.marginal,
                        frail.thresh=frail.thresh)
    marginal$value.nolam <- marginal$value
    marginal$value <- marginal$value + scale.marginal*sum(log(lambda_i[,1])*d_1 + log(lambda_i[,2])*d_2)

    theta <-  marginal$par[1:(num.a+num.g)]
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

    z <- lnz_EM(alpha = alpha, gamma = gamma, ag_i = ag_idx, g_i = g_idx,
                a_i = a_idx, PVF_i = PVF_idx, s = s, d = d.marg, frailty = frailty, n = n_C,
                ID = ID, iter = iter)$z

    # convergence condition met?
    diff.marginal <- (marginal$value - old.marginal$value)/scale.marginal
    marginal.worsened <- diff.marginal < 0
    diff.cond <- (cond$value - old.cond$value)/scale.cond

    max.cond <- diff.marginal >= converge | marginal.worsened
    if(prnt) {print(model_par); print(paste('loglihood = ', marginal$value))}
    if(iter == iterlim & max.cond){iter <- iter+1; break}
  }

  return(list(marginal = marginal, cond = cond, beta = beta, iter = iter,
              diff.cond = diff.cond, diff.marginal = diff.marginal, z = z,
              lambda_i = lambda_i, theta = theta, model_par = model_par, uni.logli=uni.logli,
              alpha = alpha, gamma = gamma,  lambda = lambda, lnlambda = lnlambda,
              partial.lambda = partial.lambda, s = s, tr.s = tr.s,
              sigma = sigma, eta = eta, nu = nu))
}
