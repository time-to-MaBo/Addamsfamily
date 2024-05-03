marginal.cs <- function(par,tau=NULL,s,ag_i,a_i,g_i,PVF_i,ag.idx,a.idx,g.idx,PVF.idx,
                        IG.idx,n_C,num.a,num.g,J,alpha.idx,gamma.idx,
                        aNBp.idx,aNB.idx,aB.idx=FALSE,N=NA,cP.idx,H.idx,n_xx,overdisp,n.date=NULL,
                        scale.marginal,sign.info,idx,frail.thresh,
                        p ## add discrete
){

  model_par <- model.par(theta = par, num.a = num.a, num.g = num.g, J = J,
                         alpha.idx = alpha.idx, gamma.idx = gamma.idx, a.idx = a.idx,
                         PVF.idx = PVF.idx, IG.idx = IG.idx,
                         aNBp.idx=aNBp.idx,aNB.idx=aNB.idx,aB.idx=aB.idx,N=N,cP.idx=cP.idx,H.idx=H.idx)
  alpha <- model_par[idx,'alpha']
  gamma <- model_par[idx,'gamma']

  p <- p[idx]

  pi0 <- matrix(NA, nrow = n_C, ncol = dim(s)[2])

  # Estimates tend towards no frailty model if
  ## PVF: alpha -> -1, gamma -> infinity
  ## Addams: alpha -> -infty, gamma <- 0
  ## aP: alpha -> 0
  ## tmp handles no frailty model
  tmp <- (gamma < frail.thresh & (ag_i | g_i)) |
    (abs(alpha^-1) < frail.thresh & alpha < 0 &  ag_i) |
    (alpha < frail.thresh & a_i) |
    (alpha + 1 < frail.thresh &  PVF_i) |
    (gamma^-1 < frail.thresh &  (PVF_i|ag_i)) # neu: ag_i, erst ab diskrete Modelle

  if(any(tmp)){
    warning('Estimates approach no-frailty model.
            For numerical stability full conditional loglihood was calculated for marginal loglihood.')
    PVF_i[tmp] <- ag_i[tmp] <- g_i[tmp] <- a_i[tmp] <- FALSE
    uni <- unique(idx[tmp])
    ag.idx[uni] <- a.idx[uni] <- g.idx[uni] <- PVF.idx[uni] <- FALSE
  }

  # Addams model approaches gamma distribution (limiting case)
  tmp2 <- ag_i & abs(alpha)<frail.thresh
  if(any(tmp2)){
    ag_i[tmp2] <- FALSE
    g_i[tmp2] <- TRUE
    warning('Addams model approached gamma distribution. Gamma formulas were used. Only implemented for current status data.')
  }
  # Addams model approaches poisson distribution (limiting case)
  tmp2 <- ag_i & abs(alpha-gamma)<frail.thresh
  if(any(tmp2)){
    ag_i[tmp2] <- FALSE
    a_i[tmp2] <- TRUE
    warning('Addams model approached poisson distribution. Poisson formulas were used. Only implemented for current status data.')
  }

  pi0[,1] <- exp( (lnS(alpha = alpha, gamma = gamma, s = s[,1], n = n_C,
                       idx_ag = ag_i, idx_a = a_i, idx_g = g_i, idx_PVF = PVF_i) - p*s[,1] )*(s[,1]>0) ) ## add discrete ----
  pi0[tmp,1] <- exp(-s[tmp,1]*(1+p[tmp])*(s[tmp,1]>0)) ## add discrete ----

  pi0[,2:dim(s)[2]]  <- apply(X = s[,2:dim(s)[2],drop=FALSE], MARGIN = 2, function(H0){
    exp(lnS(alpha = alpha, gamma = gamma, s = H0+s[,1], n = n_C,
            idx_ag = ag_i, idx_a = a_i, idx_g = g_i, idx_PVF = PVF_i)-p*(H0+s[,1]))*(H0>0)
  })
  if(any(tmp)){
    pi0[tmp,2:dim(s)[2]]  <- apply(X = s[tmp,2:dim(s)[2],drop=FALSE], MARGIN = 2, function(H0){
      exp(-(s[tmp,1]+H0[tmp])*(1+p[tmp]))*(H0[tmp]>0)
    })
  }
  # for(k in 2:dim(s)[2]){
  #   pi0[,k] <- exp(lnS(alpha = alpha, gamma = gamma, s = s[,k]+s[,1], n = n_C,
  #                      idx_ag = ag_i, idx_a = a_i, idx_g = g_i, idx_PVF = PVF_i)-p*(s[,k]+s[,1]))*(s[,k]>0)
  #   pi0[tmp,k] <- exp(-(s[tmp,1]+s[tmp,k])*(1+p[tmp]))*(s[tmp,k]>0)
  # }

  # special handling of infinite variation model (tmp with new meaning)
  ## ag: alpha -> Inf, gamma -> Inf and alpha !-> -Inf
  ## a: alpha -> Inf
  ## g: gamma -> Inf
  ## PVF: alpha -> Inf and gamma !-> Inf, gamma -> 0
  ## only implemented for current staus data yet!!!
  tmp <- cbind(
    (PVF_i & abs(alpha^-1) < frail.thresh)*1,
    (PVF_i & gamma < frail.thresh & abs(alpha^-1) >= frail.thresh)*2,
    (g_i & gamma^-1 < frail.thresh)*2,
    ((ag_i|a_i) & alpha>0 & abs(alpha)^-1 < frail.thresh)*1,
    (ag_i & abs(gamma)^-1 < frail.thresh & abs(alpha^-1) >= frail.thresh)*2)
  tmp <- rowSums(tmp)
  if(any(as.logical(tmp))){
    warning(paste('Estimates approach infinite variation frailty model.',
                  'Loglihood approaches negative infinity (-sqrt(.Machine$double.xmax) imputed*val)*(.97,1).'),
            'val=max(frail.thresh,1e-5) for PVF & gamma -> 0, 1 else.(Otherwise difference value is non-finite).',
            'Only considered in current status data so far!')
    idx <- idx[as.logical(tmp)]
    tmp <- tmp[as.logical(tmp)]
    alpha <- model_par[idx,'alpha']
    gamma <- model_par[idx,'gamma']
    alpha[tmp==2] <- gamma[tmp==2]
    alpha <- unique(alpha)
    alpha[is.infinite(alpha)] <- .Machine$double.xmax
    marg <- (-sqrt(.Machine$double.xmax)) *
      ( sum(abs(log(alpha)))/(sum(abs(log(alpha)))+length(alpha)) )

    return(scale.marginal*marg)
  }

  # back to regular case
  pi0 <- rowSums(pi0*sign.info)

  if(any(pi0<=0)){
    clusters <- paste0(which(pi0<0),sep='\',',collapse = '')
    warning(paste('Likelihood had negative values. This is theoretically impossible.
            This was the case for contributions of the ',
                  clusters, 'th ordered cluster (caution: only = cluster ID if ID = 1:n_C).',
                  '
            -.Machine$double.xmin was imputed. Likely reason:',
                  'This usually(!) happens if the hazards are extremely small (initialised)
            and if there are more than two dependent variates. In such a case it is virtually
            impossible to have an event, let alone more than one or two,
            and for those clusters the true Likelihood is
            basically zero. It occurs that the computed value is slightly below zero.
            I advise you to countercheck your results with initialized ln{hazard}-parameters}
            that are closer to zero than the standard initialization.',sep=''))
    pi0[pi0<=0] <- .Machine$double.xmin
  }

  if(overdisp){
    tau <- tau[idx]
    tau2 <- tau[n.date[,1]]
    tmp <- tau^-1<frail.thresh
    tau2[tau2^-1<frail.thresh] <- 1
    n.date[tau2^-1<frail.thresh,2] <- 0
    marg <- lgamma(n_xx + pi0*tau) - lgamma(pi0*tau)
    if(any(tmp)){
      warning('Overdispersion model approached multinomial model!')
      marg[tmp] <- sum(log(pi0[tmp])*n_xx[tmp])
    }
    marg <- sum(marg) + sum(lgamma(tau2) - lgamma(n.date[,2] + tau2))
  }else{
    pi0 <- log(pi0)
    marg <- sum(pi0*n_xx)
  }

  return(scale.marginal*marg)

}
