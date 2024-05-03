marginal.multi <- function(par,num.a, num.g, tr.s = 0, s,
                           sum.d_i, tr.idx, n_C, J, alpha.idx, gamma.idx,
                           aNBp.idx,aNB.idx,cP.idx,H.idx, ag.idx, a.idx, g.idx, fac,
                           IG.idx, PVF.idx,ag_i, a_i, g_i, IG_i, PVF_i, idx, rank_a, rank_g,
                           type_ag, type_a, type_g, type_PVF, type_IG, scale.marginal = -1,frail.thresh){

  model_par <- model.par(theta = par, num.a = num.a, num.g = num.g, J = J, alpha.idx = alpha.idx,
                         gamma.idx = gamma.idx, a.idx = a.idx, PVF.idx = PVF.idx, IG.idx = IG.idx,
                         aNBp.idx=aNBp.idx,aNB.idx=aNB.idx,cP.idx=cP.idx,H.idx=H.idx)

  alpha <- model_par[idx,'alpha']
  gamma <- model_par[idx,'gamma']

  l_i <- tr.lnS_i <- rep(NA, n_C)

  # Estimates tend towards no frailty model if
  ## PVF: alpha -> -1, gamma -> infinity
  ## Addams: alpha -> -infty, gamma <- 0
  ## aP: alpha -> 0
  tmp <- (gamma < frail.thresh & (ag_i | g_i)) |
    (abs(alpha^-1) < frail.thresh & alpha < 0 &  ag_i) |
    (alpha < frail.thresh & a_i) |
    (alpha + 1 < frail.thresh &  PVF_i) |
    (gamma^-1 < frail.thresh &  PVF_i)
  if(any(tmp)){
    warning('Estimates approach no-frailty model.
            For numerical stability full conditional loglihood was calculated for marginal loglihood.')
    PVF_i[tmp] <- ag_i[tmp] <- g_i[tmp] <- a_i[tmp] <- FALSE
    uni <- unique(idx[tmp])
    ag.idx[uni] <- a.idx[uni] <- g.idx[uni] <- PVF.idx[uni] <- FALSE
    l_i[tmp] <- exp(-s[tmp])*(-1)^sum.d_i[tmp]
  }

  if(any(ag.idx|g.idx|a.idx|PVF.idx)){
    tr.lnS_i <- lnS( alpha = alpha, gamma = gamma, s = tr.s, n = n_C,
                     idx_ag = ag_i,idx_a = a_i, idx_g = g_i, idx_PVF = PVF_i)*tr.idx
  }
  tr.lnS_i[tmp] <- -tr.s[tmp]

  if(any(ag.idx)){
    l_i[ag_i] <- ag.Laplace(a = alpha[ag_i], g = gamma[ag_i], s = s[ag_i], order = sum.d_i[ag_i], fac = fac)
  }
  if(any(a.idx)){
    l_i[a_i] <- a.Laplace(a = alpha[a_i], order = sum.d_i[a_i], s = s[a_i], idx = a_i[a_i])
  }
  if(any(g.idx)){
    l_i[g_i] <- g.Laplace(g = gamma[g_i], order = sum.d_i[g_i], s = s[g_i] )
  }
  if(any(PVF.idx)){
    l_i[PVF_i] <- PVF.Laplace(a = alpha[PVF_i], g = gamma[PVF_i], s = s[PVF_i],
                              order = sum.d_i[PVF_i], idx = PVF_i[PVF_i]) # habe hier und bei alpha bei idx noch [PVF_i] hinzugef?gt. Bei z checken! -----
  }

  l_i <- log(l_i*(-1)^sum.d_i)

  ll <- sum(l_i - tr.lnS_i)

  # special handling of infinite variation model (tmp with new meaning)
  ## ag: alpha -> Inf, gamma -> Inf and alpha !-> -Inf
  ## a: alpha -> Inf
  ## g: gamma -> Inf
  ## PVF: alpha -> Inf and gamma !-> Inf, gamma -> 0
  tmp <- cbind(
    (PVF_i & abs(alpha^-1) < frail.thresh)*1,
    (PVF_i & gamma < frail.thresh & abs(alpha^-1) >= frail.thresh)*2,
    (g_i & gamma^-1 < frail.thresh)*2,
    ((ag_i|a_i) & alpha>0 & abs(alpha)^-1 < frail.thresh)*1,
    (ag_i & abs(gamma)^-1 < frail.thresh & abs(alpha^-1) >= frail.thresh)*2)
  tmp <- rowSums(tmp)
  if(any(as.logical(tmp))){
    warning('Estimates approach infinite variation frailty model.
            Loglihood approaches negative infinity (-sqrt(.Machine$double.xmax) imputed)*(.97,1).')

    idx <- idx[as.logical(tmp)]
    tmp <- tmp[as.logical(tmp)]
    alpha <- model_par[idx,'alpha']
    gamma <- model_par[idx,'gamma']
    alpha[tmp==2] <- gamma[tmp==2]
    alpha <- unique(alpha)
    alpha[is.infinite(alpha)] <- .Machine$double.xmax
    ll <- -sqrt(.Machine$double.xmax) * ( sum(abs(log(alpha)))/(sum(abs(log(alpha)))+1) )

    return(scale.marginal*ll)
  }

  if(is.nan(ll)){
    if( any(is.infinite(l_i[sum.d_i == 0]) & is.infinite(tr.lnS_i[sum.d_i == 0]) ) ){
      ll <- sum(l_i[sum.d_i>0]) + sum(l_i[sum.d_i == 0][is.finite(l_i[sum.d_i == 0])]) -
        sum(tr.lnS_i[is.finite(tr.lnS_i)])
      warning('ln(S) and ln(tr.S) are infinite. Terms omitted.')
    }
    if(any( is.infinite(l_i[sum.d_i>0]) |
            any(is.infinite(l_i[sum.d_i==0]) & is.finite(tr.lnS_i[sum.d_i==0]) ) )){
      ll <- sqrt(.Machine$double.xmax)*-1
      warning('Loglihood approached infinity')
    }
  }

  return(scale.marginal*ll)
}
