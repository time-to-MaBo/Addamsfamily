partial.multi <- function(iter = 0, max.cond = TRUE, marginal = list(value = Inf), cond = list(loglik = Inf),
                          X, iterlim, cox_control = Surv::coxph.control(iter.max = iterlim),
                          tr, y, d, cox.surv = Surv::Surv(tr,y,d), stratum, lnz = log(z), z, v, beta, theta,
                          dt_1, tr.dt_1, m_1,cum.numh.p, numh.p, num.lam.p, breaks,p.strata,
                          Q.partial, n.p, partial.stratum,X.tilde,X.tilde_C,initial.z,start_dist,
                          num.a, num.g, d.marg, tr.idx,n_C, K, n_i, max.n_i, sum.d_i,
                          ag.idx, a.idx, g.idx, IG.idx, PVF.idx, Addams.idx,aNBp.idx,aNB.idx,cP.idx,H.idx,
                          ag_idx, a_idx, g_idx, IG_idx, PVF_idx, J, slope.model,
                          idx, frailty, alpha.idx, gamma.idx, slope,lnphipsi,N,lnpi, fac,
                          alpha.rank, gamma.rank, ag.type, a.type, g.type, ID, prnt,frail.thresh,
                          IG.type, PVF.type, scale.marginal, m.options, R.options, converge){

  m.method <- m.options$method
  Richardson <- m.options$Richardson
  if(Richardson & !(m.method %in% c('nlm','optimize'))){
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
    m.gr <- m.Hessian <- NULL
  }

  v.Domain <- NULL
  pi <- NULL
  while(max.cond){

    iter <- iter + 1
    old.marginal <- marginal
    old.cond <- cond

    offset <- lnz + v*X.tilde
    z <- exp(lnz + v*X.tilde)

    if(K>0){
      cond <- survival::agreg.fit(x = X, y = cox.surv,  strata = stratum,
                        offset = lnz, init = beta, control = cox_control,
                        method = 'breslow', rownames = NULL)

      beta <- cond$coefficients[1:K]
      exb <- as.vector(exp(X%*%beta))
    }else{
      exb <- as.vector(X+1)
    }
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

    s <- as.vector(tapply(X = Lambda_i[,1], INDEX = ID, FUN = sum, simplify = TRUE))
    tr.s <- as.vector(tapply(X = tr.Lambda_i, INDEX = ID, FUN = sum, simplify = TRUE))
    lambda_i <- tapply(X = Lambda_i[,2], INDEX = ID, FUN = function(x) c(x, rep(0, max.n_i-length(x))), simplify = TRUE)
    lambda_i <- matrix(unlist(lambda_i), ncol = max.n_i, byrow = TRUE)
    lambda_i <- lambda_i*d.marg
    lambda_i[lambda_i == 0] <- 1
    lnlambda_i <- rowSums(log(lambda_i))

    if(!slope){
      marginal <- optimax(par = theta, fn = marginal.multi, gr = m.gr, Hessian = m.Hessian,
                          options = m.options, Richardson = Richardson, R.options = R.options,
                          method = m.method,s = s, tr.s = tr.s,
                          num.a = num.a, num.g = num.g,
                          sum.d_i = sum.d_i, tr.idx = tr.idx, n_C = n_C,
                          ag.idx = ag.idx, a.idx = a.idx, g.idx = g.idx,
                          IG.idx = IG.idx, PVF.idx = PVF.idx,
                          ag_i = ag_idx, a_i = a_idx, g_i = g_idx,
                          IG_i = IG_idx, PVF_i = PVF_idx, J = J,
                          idx = idx, alpha.idx = alpha.idx, gamma.idx = gamma.idx,
                          rank_a = alpha.rank, rank_g = gamma.rank,
                          aNBp.idx=aNBp.idx,aNB.idx=aNB.idx,cP.idx=cP.idx,H.idx=H.idx,
                          type_ag = ag.type, type_a = a.type, type_g = g.type, fac = fac,
                          type_IG = IG.type, type_PVF = PVF.type, scale.marginal = scale.marginal,
                          frail.thresh=frail.thresh)
    }else{
      Lambda_i <- tapply(X = Lambda_i[,1], INDEX = ID, FUN = function(x) c(x, rep(0, max.n_i-length(x))), simplify = TRUE)
      Lambda_i <- matrix(unlist(Lambda_i), ncol = max.n_i, byrow = TRUE)
      tr.Lambda_i <- tapply(X = tr.Lambda_i, INDEX = ID, FUN = function(x) c(x, rep(0, max.n_i-length(x))), simplify = TRUE)
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

      marginal <- optimax(par = c(theta,lnpi,lnphipsi), fn = slope.opt,
                          gr = m.gr, Hessian = m.Hessian, method = m.method, options = m.options,
                          Richardson = Richardson, R.options = R.options, Lambda_i=Lambda_i,
                          tr.Lambda_i=tr.Lambda_i, n_C=n_C,X.tilde_C=X.tilde_C,d.marg=d.marg, num.a=num.a, num.g=num.g, N=N,
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
      uni.logli <- sum(lnlambda_i - s + tr.s)
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
                   s = s, frailty = frailty, n = n_C, ID = ID, sum.d_i = sum.d_i, iter = iter, fac = fac)
    }else{
      z <- vz.frail(alpha=alpha, gamma=gamma, ag_i=ag_idx, g_i=g_idx, a_i=a_idx, PVF_i=PVF_idx,
                    d.marg=d.marg, frailty = frailty, n_C = n_C, ID = ID, sum.d_i = sum.d_i, iter = iter,
                    pi=pi, v.Domain=v.Domain, N=N, X.tilde_C=X.tilde_C, Lambda_i = Lambda_i, fac = fac)
      v <- z$v
    }
    # ln{E[Z|data]} instead of E[ln{z}|data] because of link to full conditional Likelihood
    lnz <- z$lnz
    z <- z$z

    # convergence condition met?
    diff.marginal <- (marginal$value - old.marginal$value)/scale.marginal
    marginal.worsened <- diff.marginal < 0
    diff.cond <- cond$loglik[1] - old.cond$loglik[1]

    max.cond <- diff.marginal >= converge | marginal.worsened
    if(prnt) {print(model_par); print(paste('marginal loglihood = ', marginal$value))}
    if(iter == iterlim & max.cond){iter <- iter+1; break}
  }

  cond$value <- cond$loglik[1]

  return(list(marginal = marginal, cond = cond, iter = iter, beta = beta,
              diff.cond = diff.cond, diff.marginal = diff.marginal, z = z, uni.logli=uni.logli,
              v=v, v.Domain=v.Domain,lnphipsi=lnphipsi,
              lambda_i = lambda_i, theta = theta, model_par = model_par, Lambda_i=Lambda_i, lnpi=lnpi, pi=pi,
              alpha = alpha, gamma = gamma, s = s, tr.s = tr.s, sigma = NULL, eta = NULL, nu = NULL,
              partial.lambda = partial.lambda, lambda = NULL, lnlambda = NULL))
}
