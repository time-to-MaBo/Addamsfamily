profile.bi <- function(iter = 0, max.cond = TRUE, marginal = list(value = Inf), cond = list(value = Inf),
                       X, cond.scale, K,
                       tr, y, d, stratum, z, beta, theta,
                       scale.cond,dt_1, tr.dt_1, m_1,cum.numh.p, numh.p, num.lam.p, breaks,p.strata,
                       Q.partial, n, n.p, partial.stratum,
                       num.a, num.g,d_1, d_2, d.marg, tr.idx,n_C,
                       ag.idx, a.idx, g.idx, IG.idx, PVF.idx, Addams.idx,aNBp.idx,aNB.idx,cP.idx,H.idx,
                       ag_idx, a_idx, g_idx, IG_idx, PVF_idx, J,
                       idx, frailty, alpha.idx, gamma.idx,initial.z,start_dist,
                       alpha.rank, gamma.rank, ag.type, a.type, g.type, ID,iterlim,frail.thresh,
                       IG.type, PVF.type, scale.marginal,c.options,m.options,R.options,converge,prnt){

  c.method <- c.options$method
  m.method <- m.options$method
  Richardson <- m.options$Richardson
  if(Richardson){
    c.Hessian <- function(par,...){
      numDeriv::jacobian(func = dpartial.db, x = par, method = 'Richardson',
               side = NULL, method.args = R.options, ...)
    }
    m.Hessian <- function(par,...){
      numDeriv::jacobian(func = dmarginal.dtheta, x = par, method = 'Richardson',
               side = NULL, method.args = R.options, ...)
    }
  }else{
    m.Hessian <- c.Hessian <- NULL
  }

  while(max.cond){

    iter <- iter + 1
    old.marginal <- marginal
    old.cond <- cond

    cond <- optimax(par = beta, fn = partial.loglihood, gr = dpartial.db,
                    Hessian = c.Hessian, options = c.options, Richardson = FALSE,
                    R.options = NULL, method = c.method,
                    z=z,y=y,d=d,X=X,dt=dt_1,tr.dt=tr.dt_1,m=m_1,n=n,K=K,
                    scale.cond = scale.cond)

    beta <- cond$par[1:K]
    exb <- as.vector(exp(X%*%beta))

    psi_i <- as.vector(z*exb)

    partial.lambda <- piecewise_haz(dt=dt_1,tr.dt=tr.dt_1,m=m_1,
                                    d=d,psi_i=psi_i,
                                    idx.lam=cum.numh.p,numh=numh.p,
                                    num.lam=num.lam.p,breaks=breaks,strata=p.strata,
                                    Qstar=Q.partial)

    lambda <- unlist(partial.lambda)
    Lambda <- colSums(dt_1*lambda)
    Lambda <- cbind(Lambda,(m_1*lambda)[m_1==1])
    tr.Lambda <- colSums(tr.dt_1*lambda)

    Lambda_i <- Lambda*exb
    tr.Lambda_i <- tr.Lambda*exb

    # marginal model
    s <- as.vector(tapply(X = Lambda_i[,1], INDEX = ID, FUN = sum, simplify = TRUE))
    tr.s <- as.vector(tapply(X = tr.Lambda_i, INDEX = ID, FUN = sum, simplify = TRUE))
    lambda_i <- tapply(X = Lambda_i[,2], INDEX = ID, FUN = function(x) x, simplify = FALSE)
    lambda_i <- lapply(X = lambda_i, FUN = function(x) if(length(x)==2){x}else{c(x,0)})
    lambda_i <- matrix(unlist(lambda_i), ncol = 2, byrow = TRUE)

    marginal <- optimax(par = theta, fn = marginal.bi, gr = dmarginal.dtheta,
                        Hessian = m.Hessian, options = m.options, Richardson = FALSE,
                        R.options = NULL, method = m.method,
                        s = s, tr.s = tr.s,
                        num.a = num.a, num.g = num.g,
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
      uni.logli <- sum(d_1*log(lambda_i[,1]) + d_2*log(lambda_i[,2]) - s + tr.s)
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
    marginal.worsened = diff.marginal < 0
    diff.cond <- (cond$value - old.cond$value)/scale.cond

    max.cond <- diff.marginal >= converge | marginal.worsened
    if(prnt) {print(model_par); print(paste('loglihood = ', marginal$value))}
    if(iter == iterlim & max.cond){iter <- iter+1; break}
  }

  return(list(marginal = marginal, cond = cond, iter = iter, beta = beta,
              diff.cond = diff.cond, diff.marginal = diff.marginal, z = z,
              lambda_i = lambda_i, theta = theta, model_par = model_par, uni.logli = uni.logli,
              alpha = alpha, gamma = gamma, s = s, tr.s = tr.s, sigma = NULL, eta = NULL, nu = NULL,
              partial.lambda = partial.lambda, lambda = NULL, lnlambda = NULL))
}
