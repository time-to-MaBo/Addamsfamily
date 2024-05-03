marginal.bi <- function(par,num.a, num.g,
                        tr.s = 0, s,
                        d_1, d_2, tr.idx, n,
                        J, alpha.idx, gamma.idx,
                        ag.idx, a.idx, g.idx,
                        aNBp.idx,aNB.idx,cP.idx,H.idx,
                        IG.idx, PVF.idx,ag_i, a_i, g_i,
                        IG_i, PVF_i, idx, rank_a, rank_g,
                        type_ag, type_a, type_g, type_PVF,
                        type_IG, scale.marginal = -1,frail.thresh){

  model_par <- model.par(theta = par, num.a = num.a, num.g = num.g, J = J, alpha.idx = alpha.idx,
                         gamma.idx = gamma.idx, a.idx = a.idx, PVF.idx = PVF.idx, IG.idx = IG.idx,
                         aNBp.idx=aNBp.idx,aNB.idx=aNB.idx,cP.idx=cP.idx,H.idx=H.idx)

  alpha <- model_par[idx,'alpha']
  gamma <- model_par[idx,'gamma']

  pi_00tr <- lnS(alpha = alpha, gamma = gamma, s = tr.s, n = n,
                 idx_ag = ag_i, idx_a = a_i, idx_g = g_i,
                 idx_PVF = PVF_i)*(tr.idx)
  pi_00 <- lnS(alpha = alpha, gamma = gamma, s = s, n = n,
               idx_ag = ag_i, idx_a = a_i, idx_g = g_i,
               idx_PVF = PVF_i)
  pop_haz1 <- lnbipop_haz(alpha = alpha, gamma = gamma,
                          s = s, lambda_i = rep(1,n), n = n,
                          idx_ag = ag_i, idx_a = a_i, idx_g = g_i,
                          idx_PVF = PVF_i)
  pop_haz2 <- lnbipop_haz(alpha = alpha, gamma = gamma,
                          s = s, lambda_i = rep(1,n), n = n,
                          idx_ag = ag_i, idx_a = a_i, idx_g = g_i,
                          idx_PVF = PVF_i)
  crf <- lnCRF(alpha = alpha, gamma = gamma, s = s, n = n,
               idx_ag = ag_i, idx_a = a_i, idx_g = g_i,
               idx_PVF = PVF_i)

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
    crf[tmp] <- 0
    pi_00tr[tmp] <- -tr.s[tmp]
    pi_00[tmp] <- -s[tmp]
    pop_haz1[tmp] <- 0
    pop_haz2[tmp] <- 0
  }

  ll <- sum(pi_00 + d_1 * pop_haz1 + d_2 * pop_haz2 +
              d_1*d_2 * crf - pi_00tr)

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
    if( any(is.infinite(pi_00) & is.infinite(pi_00tr) ) ){
      ll <- sum(d_1 * pop_haz1 + d_2 * pop_haz2 +
                  d_1*d_2 * crf)
      warning('ln(S) and ln(tr.S) are infinite. Terms omitted.')
    }
    if(any( is.infinite(c(pop_haz1, pop_haz2, crf)) |
            any(is.infinite(pi_00) & is.finite(pi_00tr) ) )){
      ll <- sqrt(.Machine$double.xmax)*-1
      warning('Loglihood approached infinity')
    }
  }

  return(scale.marginal*ll)
}
