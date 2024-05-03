GG.cs <- function(beta,lnlambda,theta,tau,K,num.lam,num.lam.p,X,d,y,tr,z,z2=z,n,n_C,n.GG,breaks,strata,cum.numh,partial.idx,
                  partial.rank,partial_i,ag_i,a_i,g_i,PVF_i,ag.idx,a.idx,g.idx,PVF.idx,Q,numh.p,overdisp,n.date,
                  IG.idx,num.a,num.g,J,alpha.idx,gamma.idx, dt_1, stratum,alpha.rank,gamma.rank,
                  idx,sum.d_i,max.n_i, cond = list(value = Inf),Addams.idx,GG_i,lam.idx,GG.idx,
                  scale.cond, scale.marginal, ID, d.marg, converge,initial.z,start_dist,
                  GG.op,GG.sc,GG.special,sc.inv,aNBp.idx,aNB.idx,aB.idx,N,cP.idx,H.idx,frail.thresh,
                  n_xx,c.n_xx,sign.info,marginal = cond,c.options,m.options,R.options,prnt,
                  iterlim,cs.options,iter=0,
                  num.p,p.idx## add discrete ----
){
  # caution Q = Q.GG, stratum = GG.stratum, y = y.GG!

  p <- rep(0,J)
  lnp <- rep(0,num.p)

  EM.cs <- max.cond <- cs.options$EM.cs
  Richardson <- m.options$Richardson
  m.method <- m.options$method
  if(EM.cs){
    if(is.null(attr(z,'z2'))){z2 <- z}else{z2 <- attr(z,'z2')}
    c.method <- c.options$method
    if(Richardson & !(c.method %in% c('nlm','optimize'))){
      c.gr <- function(par,...){
        numDeriv::grad(func = cond.GG.cs, x = par, method = 'Richardson',
             side = NULL, method.args = R.options, ...)
      }
      c.Hessian <- function(par,...){
        numDeriv::hessian(func = cond.GG.cs, x = par, method = 'Richardson',
                method.args = R.options, ...)
      }
    }else{
      c.gr <- NULL
      c.Hessian <- NULL
    }
  }
  if(Richardson & !(m.method %in% c('nlm','optimize'))){
    func <- ifelse(EM.cs,marginal.cs,dcs.multi)
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

  if(!EM.cs){
    if(cs.options$univariate){
      if(all(is.na(cs.options$uni.model))){
        cond <- loglihood.uni(X = X, iterlim = iterlim, tr = y*0, y = y, d = d, stratum = NULL, stratum.GG = stratum,
                              beta = beta, dt_1 = dt_1, tr.dt_1 = dt_1*0, m_1 = NULL, GG_i = GG_i,
                              cum.numh.p = NULL, numh.p = numh.p, num.lam.p = num.lam.p, num.lam = num.lam,
                              breaks = breaks, p.strata = NULL, Q.partial = NULL, Q.GG = Q,
                              K = K, no.piece = NULL, partial.cond = any(partial_i), n.p = n-n.GG, n.GG = n.GG,
                              lnlambda = lnlambda, GG.op = GG.op, n_xx = c.n_xx,
                              partial.rank=partial.rank, cum.numh=cum.numh,
                              sc.inv=sc.inv, current.status=TRUE,lam.idx=lam.idx,
                              GG.sc = GG.sc, GG.special = GG.special, c.options = c.options, R.options = R.options,
                              partial.idx = partial.idx, partial_i = partial_i,
                              scale.cond = scale.cond, GG.strata = NULL) # caution Q,stratum is x.GG! NULLs not needed in cs case
        cond$univariate <- 'Caution!!! As you have chosen to optimize the marginal loglihood directly and not via the EM algorithm,
                            this is the univariate model and NOT a conditional model, i.e. z=1!!!'
      }else{
        cond <- cs.options$uni.model
      }
    }else{
      cond <- list(univariate = 'Not calculated as you have chosen univariate=FALSE in cs.options list.',
                   cond = list(value=NA))
      if(cs.options$uni.initial){
        warning('You have chosen not to calculate a univariate model (univariate=FALSE) and also not to
                supply one(uni.model=NA). Nevertheless you want the parameters to be
                initialised with the univariate estimates (uni.initial=TRUE). This makes
                obviously no sense. Hence, uni.initial was set to FALSE! Maybe you wanted
                the marginal loglihood to be optimized by the EM Algorithm?
                If so, set EM.cs = TRUE and re-estimate the model.
                (This refers to input list cs.options.)')
        cs.options$uni.initial <- FALSE
      }
    }

    uni.logli <- cond$cond$value[1]/scale.cond
    if(cs.options$uni.initial){
      beta <- cond$beta
      lnlambda <- cond$lnlambda
    }

    marginal <- optimax(par = c(theta,beta,lnlambda,tau['lntau',]), fn = dcs.multi, method = m.method,
                        options = m.options, Richardson = Richardson, R.options = R.options,
                        gr = m.gr, Hessian = m.Hessian,n_xx=n_xx,overdisp=overdisp,n.date=n.date,
                        K=K,num.lam=num.lam,num.lam.p=num.lam.p,X=X,y=y,n=n,n_C=n_C,n.GG=n.GG,
                        partial.idx=partial.idx,GG.idx=GG.idx, sign.info = sign.info,
                        partial_i=partial_i,ag_i=ag_i,a_i=a_i,g_i=g_i,PVF_i=PVF_i,
                        ag.idx=ag.idx,a.idx=a.idx,g.idx=g.idx,PVF.idx=PVF.idx,Q=Q,
                        IG.idx=IG.idx,tr.idx=tr.idx,num.a=num.a,num.g=num.g,J=J,alpha.idx=alpha.idx,gamma.idx=gamma.idx, dt_1=dt_1, stratum=stratum,
                        idx=idx,sum.d_i=sum.d_i,max.n_i=max.n_i,GG_i=GG_i,lam.idx=lam.idx,
                        frail.thresh=frail.thresh,scale.marginal=scale.marginal, ID=ID, d.marg=d.marg,
                        GG.op=GG.op,GG.sc=GG.sc,GG.special=GG.special,sc.inv=sc.inv,aNBp.idx=aNBp.idx,aNB.idx=aNB.idx,aB.idx=aB.idx,N=N,
                        cP.idx=cP.idx,H.idx=H.idx,
                        num.p=num.p,p.idx=p.idx## add discrete ----
    )

    theta <-  marginal$par[1:(num.a+num.g)]
    model_par <- model.par(theta = theta, num.a = num.a, num.g = num.g, J = J, alpha.idx = alpha.idx,
                           gamma.idx = gamma.idx, a.idx = a.idx, PVF.idx = PVF.idx, IG.idx = IG.idx,
                           aNBp.idx=aNBp.idx,aNB.idx=aNB.idx,aB.idx=aB.idx,N=N,cP.idx=cP.idx,H.idx=H.idx)
    alpha <- model_par[,'alpha']
    gamma <- model_par[,'gamma']
    if(num.p>0){
      lnp <- marginal$par[(num.a+num.g+1):(num.a+num.g+num.p)]
      p[p.idx] <- exp(lnp)
    } ## add discrete ----
    beta <-  marginal$par[(num.a+num.g+num.p+1):(num.a+num.g+num.p+K)]## add discrete ----
    exb <- as.vector(exp(X%*%beta))
    lnlambda.p <-   marginal$par[(num.a+num.g+num.p+K+1):(num.a+num.g+num.p+K+num.lam+num.lam.p)]## add discrete ----
    lnlambda <- lnlambda.p[lam.idx]
    lnlambda.p[lam.idx] <- -Inf

    if(Q>0){
      tmp <- seq(from=1,to=3*Q,by=3) # caution Q is Q.GG
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

    Lambda <- tr.Lambda <- rep(NA, n)
    Lambda[GG_i] <- GG.haz(s = sigma, e = eta, nu = nu, y = y, str = stratum, n = n.GG)[,1]

    lambda.p <- exp(lnlambda.p)
    if(any(partial.idx)){
      lambda.p[lambda.p==0 & !lam.idx] <- .Machine$double.xmin
      Lambda[partial_i] <- colSums(dt_1*(lambda.p[!lam.idx]))
    }
    lnlambda.p[lam.idx] <- lnlambda
    lnlambda <- lnlambda.p

    Lambda_i <- Lambda*exb
    Lambda_i <- tapply(X = Lambda_i, INDEX = ID, FUN = function(x) c(x,rep(0,max.n_i-length(x))),
                       simplify = FALSE)
    Lambda_i <- t(sapply(X = Lambda_i, FUN = function(x) unlist(x)))

    # s <- matrix( 0, nrow = n_C, ncol = 1+sum( choose(max(sum.d_i), 1:max(sum.d_i)) ) )
    # sign.info <- s + 1
    # s[,1] <- rowSums(Lambda_i*(1-d.marg))
    # for(i in 1:n_C){
    #   d_i <- sum.d_i[i]
    #   if(d_i>0){
    #     tmp.Lambda <- Lambda_i[i,d.marg[i,] == 1]
    #     for(j in 1:d_i){
    #       tmp <- as.matrix(combn(x = d_i, m = j))
    #       tmp2 <- apply(X = tmp, MARGIN = 2, FUN = function(x) sum(tmp.Lambda[x]))
    #       s[i,(2+sum(choose(d_i, 1:(j-1))*(j>1)) ) : (1+sum(choose(d_i, 1:(j))))] <- tmp2
    #       if(j%%2!=0){
    #         sign.info[i,(2+sum(choose(d_i, 1:(j-1))*(j>1)) ) : (1+sum(choose(d_i, 1:j)))] <- -1
    #       }
    #     }
    #   }
    # }
    s <- matrix( 0, nrow = n_C, ncol = 1+sum( choose(max(sum.d_i), 1:max(sum.d_i)) ) )
    s[,1] <- rowSums(Lambda_i*(1-d.marg))
    for(i in 1:n_C){
      d_i <- sum.d_i[i]
      if(d_i>0){
        tmp.Lambda <- Lambda_i[i,d.marg[i,] == 1]
        for(j in 1:d_i){
          tmp <- as.matrix(utils::combn(x = d_i, m = j))
          tmp2 <- apply(X = tmp, MARGIN = 2, FUN = function(x) sum(tmp.Lambda[x]))
          s[i,(2+sum(choose(d_i, 1:(j-1))*(j>1)) ) : (1+sum(choose(d_i, 1:(j))))] <- tmp2
        }
      }
    }

    if(overdisp){
      tau['lntau',] <- marginal$par[(num.a+num.g+num.p+K+num.lam+num.lam.p+1):length(marginal$par)]## add discrete ----
      tau['tau',] <- exp(tau['lntau',])
      tau['omega',] <- 1/(1+tau['tau',])
    }

    # Expectation step: Z
    if(any(Addams.idx)){
      g.idx[Addams.idx] <- abs(model_par[Addams.idx,'alpha'])< 1e-4
      a.idx[Addams.idx] <- abs(model_par[Addams.idx,'alpha']-model_par[Addams.idx,'gamma'])<1e-4 & !g.idx[Addams.idx]
      ag.idx[Addams.idx] <- !a.idx & !g.idx & Addams.idx
      ag_i <- ag.idx[frailty]
      a_i <- a.idx[frailty]
      g_i <- g.idx[frailty]
    }

    # Addams model approaches gamma distribution (limiting case)
    tmp2 <- ag_i & abs(alpha[idx])<frail.thresh
    if(any(tmp2)){
      ag_i[tmp2] <- FALSE
      g_i[tmp2] <- TRUE
      warning('Addams model approached gamma distribution. Gamma formulas were used.
              Only implemented for current status data. Note that this is the final result!')
    }

    if(any(p.idx)){
      z <- NA # its not just z+p. Derivatives of Laplace with p>0 involve product rule. Hence, NA imputed.
    }else{
      z <- z.cs(alpha = alpha, gamma = gamma, ag_i = ag_i, g_i = g_i, a_i = a_i, PVF_i = PVF_i,
                s = s, idx = idx, n_C = n_C, ID = ID, sign.info = sign.info)$z
    }

    # convergence condition met?
    diff.marginal <- NULL
    marginal.worsened <- NULL
    diff.cond <- NULL

  }

  while(max.cond){

    iter <- iter + 1
    old.marginal <- list(marginal = marginal)
    old.cond <- list(cond = cond)

    if(iter>1){
      old.marginal$theta <- theta
      old.marginal$model_par <- model_par
      old.marginal$alpha <- alpha
      old.marginal$gamma <- gamma
      old.cond$beta <- beta
      old.cond$lnlambda <- lnlambda
      old.cond$lambda <- lambda
      old.cond$lambda.p <- lambda.p
      old.cond$sigma <- sigma
      old.cond$eta <- eta
      old.cond$nu <- nu
      old.cond$z <- z
      old.cond$s <- s
    }

    cond <- optimax(par = c(beta,lnlambda), fn = cond.GG.cs, gr = c.gr, Hessian = c.Hessian,
                    options = c.options, method = c.method, Richardson = Richardson, R.options = R.options,
                    d = d, y = y,  stratum = stratum, Q = Q, n = n, n_xx = c.n_xx,
                    X = X, z = z, z2 = z2, K = K, num.lam = num.lam, lam.idx = lam.idx, n.GG = n.GG, partial.idx = partial.idx,
                    partial_i = partial_i, GG_i = GG_i, dt_1 = dt_1, num.lam.p = num.lam.p, GG.op = GG.op,
                    GG.sc = GG.sc, GG.special = GG.special, scale.cond = scale.cond)

    beta <- cond$par[1:K]
    exb <- as.vector(exp(X%*%beta))
    lnlambda.p <-  cond$par[(K+1):(K+num.lam+num.lam.p)]
    lnlambda <- lnlambda.p[lam.idx]
    lnlambda.p[lam.idx] <- -Inf

    if(Q>0){
      tmp <- seq(from=1,to=3*Q,by=3) # caution Q is Q.GG
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

    Lambda <- tr.Lambda <- rep(NA, n)
    Lambda[GG_i] <- GG.haz(s = sigma, e = eta, nu = nu, y = y, str = stratum, n = n.GG)[,1]

    lambda.p <- exp(lnlambda.p)
    if(any(partial.idx)){
      lambda.p[lambda.p==0 & !lam.idx] <- .Machine$double.xmin
      Lambda[partial_i] <- colSums(dt_1*(lambda.p[!lam.idx]))
    }
    lnlambda.p[lam.idx] <- lnlambda
    lnlambda <- lnlambda.p

    Lambda_i <- Lambda*exb
    Lambda_i <- tapply(X = Lambda_i, INDEX = ID, FUN = function(x) c(x,rep(0,max.n_i-length(x))),
                       simplify = FALSE)
    Lambda_i <- t(sapply(X = Lambda_i, FUN = function(x) unlist(x)))

    s <- matrix( 0, nrow = n_C, ncol = 1+sum( choose(max(sum.d_i), 1:max(sum.d_i)) ) )
    s[,1] <- rowSums(Lambda_i*(1-d.marg))
    for(i in 1:n_C){
      d_i <- sum.d_i[i]
      if(d_i>0){
        tmp.Lambda <- Lambda_i[i,d.marg[i,] == 1]
        for(j in 1:d_i){
          tmp <- as.matrix(combn(x = d_i, m = j))
          tmp2 <- apply(X = tmp, MARGIN = 2, FUN = function(x) sum(tmp.Lambda[x]))
          s[i,(2+sum(choose(d_i, 1:(j-1))*(j>1)) ) : (1+sum(choose(d_i, 1:(j))))] <- tmp2
        }
      }
    }

    marginal <- optimax(par = theta, fn = marginal.cs, gr = m.gr, Hessian = m.Hessian,
                        options = m.options, method = m.method, Richardson = Richardson, R.options = R.options,
                        s = s,num.a = num.a, num.g = num.g, n_xx=n_xx,
                        n = n_C,frail.thresh=frail.thresh,tau=NULL,overdisp=FALSE,n.date=NULL,
                        ag.idx = ag.idx, a.idx = a.idx, g.idx = g.idx,
                        IG.idx = IG.idx, PVF.idx = PVF.idx,
                        aNBp.idx=aNBp.idx,aNB.idx=aNB.idx,cP.idx=cP.idx,H.idx=H.idx,
                        ag_i = ag_i, a_i = a_i, g_i = g_i,
                        PVF_i = PVF_i, J = J, sign.info = sign.info,
                        idx = idx, alpha.idx = alpha.idx, gamma.idx = gamma.idx,
                        scale.marginal = scale.marginal, lnp = rep(-Inf,J))

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
      ag_i <- ag.idx[frailty]
      a_i <- a.idx[frailty]
      g_i <- g.idx[frailty]
    }

    z <- z.cs(alpha = alpha, gamma = gamma, ag_i = ag_i, g_i = g_i, a_i = a_i, PVF_i = PVF_i,
              s = s, idx = idx, n_C = n_C, ID = ID, sign.info = sign.info)
    z2 <- z$z2
    z <- z$z

    # convergence condition met?
    diff.marginal <- (marginal$value - old.marginal$marginal$value)/scale.marginal
    marginal.worsened <- diff.marginal < 0
    diff.cond <- (cond$value - old.cond$cond$value)/scale.cond

    max.cond <- diff.marginal >= converge | marginal.worsened
    if(prnt) {print(model_par); print(paste('loglihood = ', marginal$value))}
    if(iter == iterlim & max.cond){iter <- iter+1; break}

    if(iter>1 & marginal.worsened){
      warning(paste('Marginal loglihood decreases! This happened in global iteration ',
                    iter, '. Results of former iteration are chosen as optimaxum.',
                    'Maybe try to reestimate model with current estimates and
                     altered step length or something like that.', sep=''))
      iter <- iter-1
      marginal <- old.marginal$marginal
      cond <- old.cond$cond
      theta <- old.marginal$theta
      model_par <- old.marginal$model_par
      alpha <- old.marginal$alpha
      gamma <- old.marginal$gamma
      beta <- old.cond$beta
      sigma <- old.cond$sigma
      eta <- old.cond$eta
      nu <- old.cond$nu
      lnlambda <- old.cond$lnlambda
      lambda <- old.cond$lambda
      lambda.p <- old.cond$lambda.p
      z <- old.cond$z
      s <- old.cond$s
      break
    }

  }

  partial.lambda <- list()
  if(any(partial.idx)){
    for(q in which(partial.idx)){
      partial.lambda[[partial.rank[q]]] <- lambda.p[(cum.numh[q]+1):cum.numh[q+1]]
      lam_name <- rep(NA, numh.p[partial.rank[q]])
      for(i in 1:numh.p[partial.rank[q]]){
        lam_name[i] <- paste('(', breaks[[partial.rank[q]]][i], ',', breaks[[partial.rank[q]]][i+1] ,']', sep = '')
      }
      names(partial.lambda[[partial.rank[q]]]) <- lam_name
    }
    names(partial.lambda) <- strata[partial.idx]
  }

  if(any(GG.idx)) lambda.p[lam.idx] <- lambda[sc.inv]
  lambda <- lambda.p

  return(list(marginal = marginal, cond = cond, iter = iter, beta = beta,
              diff.cond = diff.cond, diff.marginal = diff.marginal, z = z, uni.logli=uni.logli,
              lambda_i = NULL, theta = theta, model_par = model_par, p = p, lnp = lnp,
              alpha = alpha, gamma = gamma, s = s, tau = tau,
              partial.lambda = partial.lambda, lambda = lambda, lnlambda = lnlambda,
              sigma = sigma, eta = eta, nu = nu, sign.info = sign.info))

}
