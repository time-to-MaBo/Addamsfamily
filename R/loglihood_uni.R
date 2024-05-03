loglihood.uni <- function(X, iterlim, cox_control = survival::coxph.control(iter.max = iterlim), current.status, lam.idx,
                          tr, y, d, cox.surv = survival::Surv(tr,y,d), stratum, stratum.GG, beta, num.lam,
                          dt_1, tr.dt_1, m_1,cum.numh, cum.numh.p, numh.p, num.lam.p, breaks,p.strata,
                          Q.partial, Q.GG, K, no.piece, partial.cond, n.p, n.GG, GG_i, n_xx,
                          GG.op, GG.sc, sc.inv, GG.special, lnlambda, c.options, R.options,
                          partial.idx,partial.rank, scale.cond, GG.strata, partial_i){

  Richardson <- c.options$Richardson
  method <- c.options$method

  z <- rep(1,n.p+n.GG)
  z.GG <- z[GG_i]
  z.p <- z[partial_i]

  if(current.status){
    if(partial.cond){ # no GG + cs

      if(Richardson & method %in% c('nlm','nlminb')){
        # as gradient is available analytically, Hessian is handled not within optimax for nlm!
        Hessian <- function(par,...){
          numDeriv::jacobian(func = dcs.dblam, x = par, method = 'Richardson',
                   side = NULL, method.args = R.options, ...)
        }
      }else{
        Hessian <- NULL
      }

      sigma <- eta <- nu <- NULL
      cond <- optimax(par = c(beta,lnlambda), fn = cond.p.cs, gr = dcs.dblam, Hessian = Hessian,
                      method = method, options = c.options, Richardson = FALSE, R.options = NULL,
                      z=z, z2=z, X=X,d=d,K=K,num.lam.p=num.lam.p,dt_1=dt_1,scale.cond=scale.cond, n_xx = n_xx)

      lnlambda <- cond$par[(K+1):(K+num.lam+num.lam.p)]
      lambda <- lambda.p <- exp(lnlambda)
    }else{ # GG hazard + cs

      if(Richardson & !(method %in% c('nlm','optimize'))){
        gr <- function(par,...){
          grad(func = cond.GG.cs, x = par, method = 'Richardson',
               side = NULL, method.args = R.options, ...)
        }
        Hessian <- function(par,...){
          hessian(func = cond.GG.cs, x = par, method = 'Richardson',
                  method.args = R.options, ...)
        }
      }else{
        gr <- NULL
        Hessian <- NULL
      }
      cond <- optimax(par = c(beta,lnlambda), fn = cond.GG.cs, gr = gr, Hessian = Hessian,
                      method = method, options = c.options, R.options = R.options, Richardson = Richardson,
                      d = d, y = y[GG_i], stratum = stratum.GG, Q = Q.GG, n = n.p+n.GG,
                      X = X, z = z, z2 = z, K = K, num.lam = num.lam, lam.idx = lam.idx, n.GG = n.GG, partial.idx = partial.idx,
                      partial_i = partial_i, GG_i = GG_i, dt_1 = dt_1, num.lam.p = num.lam.p, GG.op = GG.op,
                      GG.sc = GG.sc, GG.special = GG.special, scale.cond = scale.cond, n_xx = n_xx)

      lnlambda.p <-  cond$par[(K+1):(K+num.lam+num.lam.p)]
      lnlambda <- lnlambda.p[lam.idx]
      lnlambda.p[lam.idx] <- -Inf

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

      lambda.p <- exp(lnlambda.p)

      lnlambda.p[lam.idx] <- lnlambda
      lnlambda <- lnlambda.p

    }
  }else{
    if(partial.cond){ # no GG
      lambda <- lnlambda <- sigma <- eta <- nu <- NULL
      if(no.piece){ # coxph
        scale.cond <- 1
        cond <- survival::agreg.fit(x = X, y = cox.surv,  strata = stratum,
                          init = beta, control = cox_control,
                          method = 'breslow', rownames = NULL)
        cond$par <- cond$coefficients
        cond$value <- cond$loglik[1]
      }else{ # profile loglihood
        if(Richardson & !(method %in% c('optimize'))){
          Hessian <- function(par,...){
            numDeriv::jacobian(func = dpartial.db, x = par, side = NULL,
                     method = 'Richardson', method.args = R.options, ...)
          }
        }else{
          Hessian <- NULL
        }
        cond <- optimax(par = beta, fn = partial.loglihood, gr = dpartial.db, Hessian = Hessian,
                        method=method, options = c.options, R.options = R.options, Richardson = Richardson,
                        z=rep(1, n.p),y=y,d=d,X=X,dt=dt_1,tr.dt=tr.dt_1,m=m_1,n=n.p,K=K,
                        scale.cond = scale.cond)
      }
    }else{ # fully parametric loglihood

      if(Richardson & !(method %in% c('nlm','optimize'))){
        gr <- function(par,...){
          numDeriv::grad(func = cond.loglihood, x = par, method = 'Richardson',
               side = NULL, method.args = R.options, ...)
        }
        Hessian <- function(par,...){
          numDeriv::hessian(func = cond.loglihood, x = par,
                  method = 'Richardson', method.args = R.options, ...)
        }
      }else{
        gr <- NULL
        Hessian <- NULL
      }
      cond <- optimax(par = c(beta,lnlambda), fn = cond.loglihood, gr = gr, Hessian = Hessian,
                      method = method, options = c.options, R.options = R.options, Richardson = Richardson,
                      num.lam=num.lam, K=K,
                      X.GG = X[GG_i,,drop=FALSE], y.GG=y[GG_i], tr.GG = tr[GG_i], z.GG=z.GG,d.GG=d[GG_i],
                      n.GG=n.GG,X.p=X[partial_i,,drop=FALSE], z.p=z.p,d.p=d[partial_i],n.p=n.p,Q.GG=Q.GG,
                      stratum.GG=stratum.GG,dt=dt_1,m=m_1,tr.dt=tr.dt_1,
                      scale.cond = scale.cond,GG.op=GG.op,GG.sc=GG.sc,GG.special=GG.special)

      partial.lambda <- NULL

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
    }
  }

  beta <- cond$par[1:K]
  exb <- as.vector(exp(X%*%beta))[partial_i]

  if(any(partial.idx) & !current.status){
    partial.lambda <- piecewise_haz(dt=dt_1,tr.dt=tr.dt_1,m=m_1,
                                    d=d[partial_i],psi_i=exb,
                                    idx.lam=cum.numh.p,numh=numh.p,
                                    num.lam=num.lam.p,breaks=breaks,strata=p.strata,
                                    Qstar=Q.partial)
  }
  if(current.status){
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
      names(partial.lambda) <- p.strata
    }
    lambda.p[lam.idx] <- lambda[sc.inv]
    lambda <- lambda.p
  }

  return(list(marginal = NULL, cond = cond, iter = 1, beta = beta,
              theta = NULL, model_par = NULL, alpha = NULL, gamma = NULL,
              partial.lambda = partial.lambda, lambda = lambda, lnlambda = lnlambda, history = NULL,
              sigma = sigma, eta = eta, nu = nu, scale.cond = scale.cond, z = z))
}
