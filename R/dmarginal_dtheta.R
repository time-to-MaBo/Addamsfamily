dmarginal.dtheta <- function(par, d_1, d_2, n, tr.idx,
                             s, tr.s,
                             idx, ag.idx, a.idx, g.idx, IG.idx, PVF.idx,
                             num.a, num.g, J, alpha.idx, gamma.idx,
                             aNBp.idx,aNB.idx,cP.idx,H.idx,frail.thresh,
                             ag_i, a_i, g_i, PVF_i, IG_i, rank_a, rank_g,
                             type_ag, type_a, type_g, type_PVF, type_IG,
                             scale.marginal = -1
){

  # Retrieving/ Transforming Model Parameters
  model_par <- model.par(theta = par, num.a = num.a, num.g = num.g, J = J, alpha.idx = alpha.idx,
                         gamma.idx = gamma.idx, a.idx = a.idx, PVF.idx = PVF.idx, IG.idx = IG.idx,
                         aNBp.idx=aNBp.idx,aNB.idx=aNB.idx,cP.idx=cP.idx,H.idx=H.idx)

  alpha <- model_par[,'alpha']
  gamma <- model_par[,'gamma']

  alpha_i <- alpha[idx]
  gamma_i <- gamma[idx]

  s_i <- s
  tr.s_i <- tr.s
  CRF_i <- exp( lnCRF(alpha = alpha_i, gamma = gamma_i, s = s, n = n,
                      idx_ag = ag_i, idx_a = a_i, idx_g = g_i,
                      idx_PVF = PVF_i ) )
  tr.CRF_i <- exp( lnCRF(alpha = alpha_i, gamma = gamma_i, s = tr.s,  n = n,
                         idx_ag = ag_i, idx_a = a_i, idx_g = g_i,
                         idx_PVF = PVF_i) )
  RFV_i <- CRF_i - 1
  tr.RFV_i <- tr.CRF_i - 1

  ### del ln{S_i} / del alpha or alpha* for poisson model, &
  ### del ln{S_i} / del gamma*
  dlS.da <- tr.dlS.da <- matrix(0, nrow = num.a, ncol = n)
  dlS.dg <- tr.dlS.dg <- matrix(0, nrow = num.g, ncol = n)
  ### del ln{h_i1} / del alpha or alpha* for poisson model, &
  ### del ln{h_i1} / del gamma*
  dlh1.da <- matrix(0, nrow = num.a, ncol = n)
  dlh1.dg <- matrix(0, nrow = num.g, ncol = n)
  ### del ln{h_i2} / del alpha or alpha* for poisson model, &
  ### del ln{h_i2} / del gamma*
  dlh2.da <- matrix(0, nrow = num.a, ncol = n)
  dlh2.dg <- matrix(0, nrow = num.g, ncol = n)
  # Aufpassen bei Typunterscheidung. Nicht alle alpha sind aus alpha gamma
  ### del ln{CRF} / del alpha or alpha* for poisson model, &
  ### del ln{CRF} / del gamma*
  dlcrf.da <- matrix(0, nrow = num.a, ncol = n)
  dlcrf.dg <- matrix(0, nrow = num.g, ncol = n)
  # Aufpassen bei Typunterscheidung. Nicht alle alpha sind aus alpha gamma

  if(any(ag.idx)){

    a <- alpha_i[ag_i]
    g <- gamma_i[ag_i]
    rfv <- RFV_i[ag_i]
    s <- s_i[ag_i]
    #### Survival Function
    S.a <- 1/(a-g)^2 * ( (a-g)/(a-g+rfv) * (s*(g-a) + g/a - rfv/a) -
                           log(g/rfv * (a-g+rfv)/a)  )
    S.g <-  1/(a-g)^2 * (
      (a-g)*(rfv-g)/(a-g+rfv) + g*log(g/rfv * (a-g+rfv)/a)
    )
    #### pop hazard function
    d1 <- d_1[ag_i]
    h1.a <- d1*(1/a*(rfv-g)-s*rfv)/(a-g+rfv)
    h1.g <- d1*(g-rfv)/(a-g+rfv)
    d2 <- d_2[ag_i]
    h2.a <- d2*(1/a*(rfv-g)-s*rfv)/(a-g+rfv)
    h2.g <- d2*(g-rfv)/(a-g+rfv)
    #### Cross Ratio Function
    crf.a <- rfv*s/(1+rfv)*d1*d2
    crf.g <- rfv/(1+rfv)*d1*d2
    #### truncation Survival Function
    rfv <- tr.RFV_i[ag_i]
    tidx <- tr.idx[ag_i]
    s <- tr.s_i[ag_i]

    tr.S.a <- (1/(a-g)^2 * ( (a-g)/(a-g+rfv) * (s*(g-a) + g/a - rfv/a) -
                               log(g/rfv * (a-g+rfv)/a)  ))*tidx

    tr.S.g <-  1/(a-g)^2 * (
      (a-g)*(rfv-g)/(a-g+rfv) + g*log(g/rfv * (a-g+rfv)/a)
    )*tidx

    for(j in type_ag){
      col_idx <- idx == j
      vec_idx <- col_idx[ag_i]
      dlS.da[rank_a[j], col_idx] <- S.a[vec_idx]
      dlS.dg[rank_g[j], col_idx] <- S.g[vec_idx]
      tr.dlS.da[rank_a[j], col_idx] <- tr.S.a[vec_idx]
      tr.dlS.dg[rank_g[j], col_idx] <- tr.S.g[vec_idx]
      dlh1.da[rank_a[j], col_idx] <- h1.a[vec_idx]
      dlh2.da[rank_a[j], col_idx] <- h2.a[vec_idx]
      dlh1.dg[rank_g[j], col_idx] <- h1.g[vec_idx]
      dlh2.dg[rank_g[j], col_idx] <- h2.g[vec_idx]
      dlcrf.da[rank_a[j], col_idx] <- crf.a[vec_idx]
      dlcrf.dg[rank_g[j], col_idx] <- crf.g[vec_idx]
    }

  }

  #### del alpha*
  if(any(a.idx)){

    a <- alpha_i[a_i]
    rfv <- RFV_i[a_i]
    s <- s_i[a_i]
    #### Survival Function
    S <- -(rfv^(-1)*(1+a*s) - 1/a)
    #### population hazard function
    d1 <- d_1[a_i]
    h1 <- -s*a*d1
    d2 <- d_2[a_i]
    h2 <- -s*a*d2
    #### cross ratio function
    crf <- rfv*(1+a*s)/(1+rfv)*d1*d2
    #### truncation survival function
    rfv <- tr.RFV_i[a_i]
    tidx <- tr.idx[a_i]
    s <- tr.s_i[a_i]

    tr.S <- -tidx*(rfv^(-1)*(1+a*s) - 1/a)

    for(j in type_a){
      col_idx <- idx == j
      vec_idx <- col_idx[a_i]
      dlS.da[rank_a[j], col_idx] <- S[vec_idx]
      tr.dlS.da[rank_a[j], col_idx] <- tr.S[vec_idx]
      dlh1.da[rank_a[j], col_idx] <- h1[vec_idx]
      dlh2.da[rank_a[j], col_idx] <- h2[vec_idx]
      dlcrf.da[rank_a[j], col_idx] <- crf[vec_idx]
    }

  }

  #### del gamma*
  if(any(g.idx)){

    g <- gamma_i[g_i]
    rfv <- RFV_i[g_i]
    s <- s_i[g_i]
    #### Survival Function
    S <- 1/g * log(1+g*s) - s/(1+g*s)
    #### population hazard function
    d1 <- d_1[g_i]
    h1 <- -g*s/(1+g*s)*d1
    d2 <- d_2[g_i]
    h2 <- -g*s/(1+g*s)*d2
    #### cross ratio function
    crf <- d1*d2*g/(1+g)
    #### truncation survival function
    tidx <- tr.idx[g_i]
    rfv <- tr.RFV_i[g_i]
    s <- tr.s_i[g_i]
    tr.S<- (1/g * log(1+g*s) - s/(1+g*s))*tidx

    for(j in type_g){
      col_idx <- idx == j
      vec_idx <- col_idx[g_i]
      dlS.dg[rank_g[j], col_idx] <- S[vec_idx]
      tr.dlS.dg[rank_g[j], col_idx] <- tr.S[vec_idx]
      dlh1.dg[rank_g[j], col_idx] <- h1[vec_idx]
      dlh2.dg[rank_g[j], col_idx] <- h2[vec_idx]
      dlcrf.dg[rank_g[j], col_idx] <- crf[vec_idx]
    }

  }

  if(any(PVF.idx)){

    a <- alpha_i[PVF_i]
    g <- gamma_i[PVF_i]
    rfv <- RFV_i[PVF_i]
    s <- s_i[PVF_i]
    #### del alpha*
    S.a <-(g*(1 - (g/(g + s))^a)/a^2 +
             g*(g/(g + s))^a*log(g/(g + s))/a)*(a + 1)
    #### del gamma*
    S.g <- (-(1 - (g/(g + s))^a)/a +
              (g/(g + s))^a*(1/(g + s) - g/(g + s)^2)*(g + s))*g
    #### population hazard function
    d1 <- d_1[PVF_i]
    h1.a <- d1*(a + 1)*log(g/(g + s))
    h1.g <- (a + 1)*(1/(g + s) - g/(g + s)^2)*(g + s)*d1
    d2 <- d_2[PVF_i]
    h2.a <- d2*(a + 1)*log(g/(g + s))
    h2.g <- (a + 1)*(1/(g + s) - g/(g + s)^2)*(g + s)*d2
    #### cross ratio function
    crf.a <-d1*d2*(((1 + s/g)^a/g +
                      (a + 1)*(1 + s/g)^a*log(1 + s/g)/g)*
                     (a + 1)/(1 + (a + 1)*(1 + s/g)^a/g))
    crf.g <- d1*d2*((-(a + 1)*(1 + s/g)^a/g^2 -
                       (a + 1)*(1 + s/g)^a*a*s/(g^3*(1 + s/g)))*
                      g/(1 + (a + 1)*(1 + s/g)^a/g))
    #### truncation survival function
    tidx <- tr.idx[PVF_i]
    rfv <- tr.RFV_i[PVF_i]
    s <- tr.s_i[PVF_i]

    tr.S.a <- tidx*(g*(1 - (g/(g + s))^a)/a^2 +
                      g*(g/(g + s))^a*log(g/(g + s))/a)*(a + 1)
    tr.S.g <- tidx*(-(1 - (g/(g + s))^a)/a +
                      (g/(g + s))^a*(1/(g + s) - g/(g + s)^2)*(g + s))*g

    for(j in c(type_IG,type_PVF)){
      col_idx1 <- idx == j
      col_idx <- col_idx1 & !IG_i
      vec_idx <- col_idx[PVF_i]
      vec_idx1 <- col_idx1[PVF_i]
      dlS.da[rank_a[j], col_idx] <- S.a[vec_idx]
      tr.dlS.da[rank_a[j], col_idx] <- tr.S.a[vec_idx]
      dlS.dg[rank_g[j], col_idx1] <- S.g[vec_idx1]
      tr.dlS.dg[rank_g[j], col_idx1] <- tr.S.g[vec_idx1]
      dlh1.da[rank_a[j], col_idx] <- h1.a[vec_idx]
      dlh2.da[rank_a[j], col_idx] <- h2.a[vec_idx]
      dlh1.dg[rank_g[j], col_idx1] <- h1.g[vec_idx1]
      dlh2.dg[rank_g[j], col_idx1] <- h2.g[vec_idx1]
      dlcrf.da[rank_a[j], col_idx] <- crf.a[vec_idx]
      dlcrf.dg[rank_g[j], col_idx1] <- crf.g[vec_idx1]
    }

  }

  dlS.da <- rowSums(dlS.da)
  dlS.dg <- rowSums(dlS.dg)
  tr.dlS.da <- rowSums(tr.dlS.da)
  tr.dlS.dg <- rowSums(tr.dlS.dg)

  dlh1.da <- rowSums(dlh1.da)
  dlh1.dg <- rowSums(dlh1.dg)

  dlh2.da <- rowSums(dlh2.da)
  dlh2.dg <- rowSums(dlh2.dg)

  dlcrf.da <- rowSums(dlcrf.da)
  dlcrf.dg <- rowSums(dlcrf.dg)

  Score <- c( -tr.dlS.da + dlS.da + dlh1.da + dlh2.da + dlcrf.da,
              -tr.dlS.dg + dlS.dg + dlh1.dg + dlh2.dg + dlcrf.dg )

  if(any(is.nan(Score))){
    warning(paste('Element', which(is.nan(Score)), ' of Score is NaN!' ))
    cond1 <- is.infinite( c(dlS.da,dlS.dg) ) & is.infinite( c(tr.dlS.da,tr.dlS.dg) )
    if( any(cond1) ){
      Score[which(cond1)] <- c(dlh1.da + dlh2.da + dlcrf.da,
                               dlh1.dg + dlh2.dg + dlcrf.dg)[cond1]
      warning('dln{S}/d alpha or gamma and dln{tr.S}/d alpha or gamma are infinite.
              Terms omitted.')
    }
  }
  cond2 <- is.infinite(Score)
  if(any(cond2)){
    warning('Loglihood approached infinity. sqrt(.Machine$double.xmax) inserted.')
    Score[cond2] <- sqrt(.Machine$double.xmax) * sign(Score[cond2])
  }

  return(scale.marginal*Score)
}
