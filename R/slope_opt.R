slope.opt <- function(par,Lambda_i,tr.Lambda_i,tr.idx, X.tilde_C,d.marg,num.a,num.g,N,fac,
                      J,alpha.idx,gamma.idx,a.idx,PVF.idx,IG.idx,ag.idx,g.idx,aNBp.idx,aNB.idx,
                      ag_i,g_i,a_i,PVF_i,IG_i,cP.idx,H.idx,idx,sum.d_i,n_C,slope.model,scale.marginal,frail.thresh){

  pi <- exp(-exp(par[(num.a+num.g+1)]))
  if(slope.model == 'free-Willy'){
    phi <- par[(num.a+num.g+2):(num.a+num.g+N+1)]
    psi <- phi[(N/2+1):N]
    phi <- phi[1:(N/2)]
    phi <- cumsum(-exp(phi))[(N/2):1]
    psi <- cumsum(exp(psi))
  }else{
    phi <- -exp(par[num.a+num.g+2])
    phi <- ((N/2):1)*phi
    psi <- exp(par[num.a+num.g+3])
    psi <- (1:(N/2))*psi
  }

  v <- c(phi,0,psi)
  theta <- par[1:(num.a+num.g)]
  model_par <- model.par(theta = theta, num.a = num.a, num.g = num.g, J = J, alpha.idx = alpha.idx,
                         gamma.idx = gamma.idx, a.idx = a.idx, PVF.idx = PVF.idx, IG.idx = IG.idx,
                         aNBp.idx=aNBp.idx,aNB.idx=aNB.idx,cP.idx=cP.idx,H.idx=H.idx)


  alpha <- model_par[idx,'alpha']
  gamma <- model_par[idx,'gamma']

  evx <- exp(v%x%X.tilde_C)
  s.v <- evx * Lambda_i
  s.v <- rowSums(s.v)
  tr.s.v <- evx * tr.Lambda_i
  tr.s.v <- rowSums(tr.s.v)
  vxD <- matrix(v%x%rowSums(X.tilde_C*d.marg),ncol=N+1)

  l_i <- tr.lnS_i <- matrix(NA, nrow = n_C, ncol = N+1)

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
  }

  for(i in 1:(N+1)){
    prob <- stats::dbinom(x = i-1, size = N, prob = pi, log = TRUE)
    l_i[tmp,i] <- exp(-(s.v[((i-1)*n_C+1):(n_C*i)])[tmp])*(-1)^sum.d_i[tmp]
    if(any(ag.idx)){
      l_i[ag_i,i] <- ag.Laplace(a = alpha[ag_i], g = gamma[ag_i],
                                s = (s.v[((i-1)*n_C+1):(n_C*i)])[ag_i], order = sum.d_i[ag_i], fac = fac)
    }
    if(any(a.idx)){
      l_i[a_i,i] <- a.Laplace(a = alpha[a_i], order = sum.d_i[a_i],
                              s = (s.v[((i-1)*n_C+1):(n_C*i)])[a_i], idx = a_i)
    }
    if(any(g.idx)){
      l_i[g_i,i] <- g.Laplace(g = gamma[g_i], order = sum.d_i[g_i],
                              s = (s.v[((i-1)*n_C+1):(n_C*i)])[g_i] )
    }
    if(any(PVF.idx)){
      l_i[PVF_i,i] <- PVF.Laplace(a = alpha[PVF_i], g = gamma[PVF_i], s = (s.v[((i-1)*n_C+1):(n_C*i)])[PVF_i],
                                  order = sum.d_i[PVF_i], idx = PVF_i)
    }

    l_i[,i] <- log(l_i[,i]*(-1)^sum.d_i) + prob

    tr.lnS_i[,i] <- lnS(alpha = alpha, gamma = gamma, s = tr.s.v[((i-1)*n_C+1):(n_C*i)], n = n_C,
                        idx_ag = ag_i,idx_a = a_i, idx_g = g_i, idx_PVF = PVF_i)*tr.idx
    tr.lnS_i[tmp,i] <- -tr.s.v[((i-1)*n_C+1):(n_C*i)][tmp]
  }

  l_i <- l_i + vxD

  ll <- sum(l_i - tr.lnS_i)

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
    if( any(is.infinite(l_i[rep(sum.d_i == 0,N+1)]) & is.infinite(tr.lnS_i[rep(sum.d_i == 0,N+1)]) ) ){
      ll <- sum(l_i[rep(sum.d_i>0,N+1)]) + sum(l_i[rep(sum.d_i == 0,N+1)][is.finite(l_i[rep(sum.d_i == 0,N+1)])]) -
        sum(tr.lnS_i[is.finite(tr.lnS_i)])
      warning('ln(S) and ln(tr.S) are infinite. Terms omitted.')
    }
    if(any( is.infinite(l_i[rep(sum.d_i>0,N+1)]) |
            any(is.infinite(l_i[rep(sum.d_i==0,N+1)]) & is.finite(tr.lnS_i[rep(sum.d_i==0,N+1)] ) ) )){
      ll <- sqrt(.Machine$double.xmax)*-1
      warning('Loglihood approached infinity')
    }
  }

  return(scale.marginal*ll)

}
