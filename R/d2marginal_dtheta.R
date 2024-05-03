d2marginal.dtheta2 <- function(alpha, gamma, d_1, d_2, n, tr.idx,
                               s, tr.s,
                               idx, ag.idx, a.idx, g.idx, IG.idx, PVF.idx,
                               aNBp.idx,aNB.idx,cP.idx,H.idx,
                               num.a, num.g, J, alpha.idx, gamma.idx,
                               ag_i, a_i, g_i, PVF_i, IG_i,
                               rank_a, rank_g, rank_ag,
                               type_ag, type_a, type_g, type_PVF, type_IG){

  # Retrieving/ Transforming Model Parameters

  s_i <- s
  tr.s_i <- tr.s
  alpha_i <- alpha[idx]
  gamma_i <- gamma[idx]
  CRF_i <- exp( lnCRF(alpha = alpha_i, gamma = gamma_i, s = s, n = n, idx_ag = ag_i,
                      idx_a = a_i, idx_g = g_i, idx_PVF = PVF_i) )
  tr.CRF_i <- exp( lnCRF(alpha = alpha_i, gamma = gamma_i, s = tr.s,  n = n, idx_ag = ag_i,
                         idx_a = a_i, idx_g = g_i, idx_PVF = PVF_i) )
  RFV_i <- CRF_i - 1
  tr.RFV_i <- tr.CRF_i - 1
  if(any(g.idx)){
    S <- exp(lnS(alpha = NULL, gamma = gamma_i[g_i], s = s[g_i],  n = sum(g_i), idx_ag = FALSE,
                 idx_a = FALSE, idx_g = TRUE, idx_PVF = FALSE)) # nur f?r gamma!!!!
    tr.S <- exp(lnS(alpha = NULL, gamma = gamma_i[g_i], s = tr.s[g_i],  n = sum(g_i), idx_ag = FALSE,
                    idx_a = FALSE, idx_g = TRUE, idx_PVF = FALSE)) # nur f?r gamma!!!!
  }

  # 2nd Order Derivatives ----

  ## d^2l/dtheta^2 ----
  num.ag <- J - sum(is.na(rank_ag))
  d2lS.da2 <- tr.d2lS.da2 <- matrix(0, nrow = num.a, ncol = n)
  d2lS.dg2 <- tr.d2lS.dg2 <- matrix(0, nrow = num.g, ncol = n)
  d2lS.dgda <- tr.d2lS.dgda <- matrix(0, nrow = num.ag, ncol = n)
  d2lh1.da2 <- matrix(0, nrow = num.a, ncol = n)
  d2lh1.dg2 <- matrix(0, nrow = num.g, ncol = n)
  d2lh1.dgda <- matrix(0, nrow = num.ag,ncol = n)
  d2lh2.da2 <- matrix(0, nrow = num.a, ncol = n)
  d2lh2.dg2 <- matrix(0, nrow = num.g, ncol = n)
  d2lh2.dgda <- matrix(0, nrow = num.ag,ncol = n)
  d2lcrf.da2 <- matrix(0, nrow = num.a, ncol = n)
  d2lcrf.dg2 <- matrix(0, nrow = num.g, ncol = n)
  d2lcrf.dgda <- matrix(0, nrow = num.ag,ncol = n)
  # Aufpassen bei Typunterscheidung. Nicht alle alpha sind aus alpha gamma

  ### d^2 ln{S_i} / d alpha^2 od d alpha*^2
  if(any(ag.idx)){

    a <- alpha_i[ag_i]
    g <- gamma_i[ag_i]
    s <- s_i[ag_i]
    rfv <- RFV_i[ag_i]
    # Survival
    S.a <- 2*log((1 - g/a)*exp(-a*s) + g/a)/(a - g)^3 -
      2*(g*exp(-a*s)/a^2 -
           (1 - g/a)*s*exp(-a*s) - g/a^2)/((a - g)^2*((1 - g/a)*exp(-a*s) + g/a)) +
      (-2*g*exp(-a*s)/a^3 - 2*g*s*exp(-a*s)/a^2 + (1 - g/a)*s^2*exp(-a*s) +
         2*g/a^3)/((a - g)*((1 - g/a)*exp(-a*s) + g/a)) +
      - (g*exp(-a*s)/a^2 - (1 - g/a)*s*exp(-a*s) -
           g/a^2)^2/((a - g)*((1 - g/a)*exp(-a*s) + g/a)^2)
    S.g <-  2*g^2*log((1 - g/a)*exp(-a*s) + g/a)/(a - g)^3 +
      g*log((1 - g/a)*exp(-a*s) + g/a)/(a - g)^2 +
      2*g*(-exp(-a*s)*g/a + g/a)/((a - g)^2*((1 - g/a)*exp(-a*s) + g/a)) +
      (-exp(-a*s)*g/a + g/a)/((a - g)*((1 - g/a)*exp(-a*s) + g/a)) -
      (-exp(-a*s)*g/a + g/a)^2/((a - g)*((1 - g/a)*exp(-a*s) + g/a)^2)
    S.ga <- -2*g*log((1 - g/a)*exp(-a*s) + g/a)/(a - g)^3 + g*(g*exp(-a*s)/a^2 - (1 - g/a)*s*exp(-a*s) - g/a^2)/((a - g)^2*((1 - g/a)*exp(-a*s) + g/a)) - (-exp(-a*s)*g/a + g/a)/((a - g)^2*((1 - g/a)*exp(-a*s) + g/a)) + (s*exp(-a*s)*g/a + g*exp(-a*s)/a^2 - g/a^2)/((a - g)*((1 - g/a)*exp(-a*s) + g/a)) - (-exp(-a*s)*g/a + g/a)*(g*exp(-a*s)/a^2 - (1 - g/a)*s*exp(-a*s) - g/a^2)/((a - g)*((1 - g/a)*exp(-a*s) + g/a)^2)
    # pop hazard
    d1 <- d_1[ag_i]
    h1.a <- d1/(a-g+rfv) * (
      (1/a-s)*s*rfv - 1/a^2*(rfv-g) -
        (1/a*(rfv-g)-s*rfv)*(1+s*rfv)/(a-g+rfv)
    )
    h1.g <- d1*(g-rfv)/(a-g+rfv)*(1+(g-rfv)/(a-g+rfv))
    h1.ga <- d1*(1/a*(rfv-g) - s*rfv)/(a-g+rfv) *
      (1-(rfv-g)/(a-g+rfv))
    d2 <- d_2[ag_i]
    h2.a <- d2/(a-g+rfv) * (
      (1/a-s)*s*rfv - 1/a^2*(rfv-g) -
        (1/a*(rfv-g)-s*rfv)*(1+s*rfv)/(a-g+rfv)
    )
    h2.g <- d2*(g-rfv)/(a-g+rfv)*(1+(g-rfv)/(a-g+rfv))
    h2.ga <- d2*(1/a*(rfv-g) - s*rfv)/(a-g+rfv) *
      (1-(rfv-g)/(a-g+rfv))
    # CRF
    crf.a <- d1*d2*s^2*rfv/(1+rfv) * (1 - rfv/(1+rfv))
    crf.g <- d1*d2*rfv/(1+rfv) * (1 - rfv/(1+rfv))
    crf.ga <- d1*d2*rfv*s/(1+rfv) * (1- rfv/(1+rfv))
    # Survival truncation
    s <- tr.s_i[ag_i]
    rfv <- tr.RFV_i[ag_i]
    tidx <- tr.idx[ag_i]
    tr.S.a <- tidx*(
      2*log((1 - g/a)*exp(-a*s) + g/a)/(a - g)^3 -
        2*(g*exp(-a*s)/a^2 -
             (1 - g/a)*s*exp(-a*s) - g/a^2)/((a - g)^2*((1 - g/a)*exp(-a*s) + g/a)) +
        (-2*g*exp(-a*s)/a^3 - 2*g*s*exp(-a*s)/a^2 + (1 - g/a)*s^2*exp(-a*s) +
           2*g/a^3)/((a - g)*((1 - g/a)*exp(-a*s) + g/a)) +
        - (g*exp(-a*s)/a^2 - (1 - g/a)*s*exp(-a*s) -
             g/a^2)^2/((a - g)*((1 - g/a)*exp(-a*s) + g/a)^2)
    )
    tr.S.g <- tidx*(
      2*g^2*log((1 - g/a)*exp(-a*s) + g/a)/(a - g)^3 +
        g*log((1 - g/a)*exp(-a*s) + g/a)/(a - g)^2 +
        2*g*(-exp(-a*s)*g/a + g/a)/((a - g)^2*((1 - g/a)*exp(-a*s) + g/a)) +
        (-exp(-a*s)*g/a + g/a)/((a - g)*((1 - g/a)*exp(-a*s) + g/a)) -
        (-exp(-a*s)*g/a + g/a)^2/((a - g)*((1 - g/a)*exp(-a*s) + g/a)^2)
    )
    tr.S.ga <- tidx*(
      -2*g*log((1 - g/a)*exp(-a*s) + g/a)/(a - g)^3 + g*(g*exp(-a*s)/a^2 - (1 - g/a)*s*exp(-a*s) - g/a^2)/((a - g)^2*((1 - g/a)*exp(-a*s) + g/a)) - (-exp(-a*s)*g/a + g/a)/((a - g)^2*((1 - g/a)*exp(-a*s) + g/a)) + (s*exp(-a*s)*g/a + g*exp(-a*s)/a^2 - g/a^2)/((a - g)*((1 - g/a)*exp(-a*s) + g/a)) - (-exp(-a*s)*g/a + g/a)*(g*exp(-a*s)/a^2 - (1 - g/a)*s*exp(-a*s) - g/a^2)/((a - g)*((1 - g/a)*exp(-a*s) + g/a)^2)
    )

    for(j in type_ag){
      col_idx <- idx == j
      vec_idx <- col_idx[ag_i]
      d2lS.da2[rank_a[j], col_idx] <- S.a[vec_idx]
      tr.d2lS.da2[rank_a[j], col_idx] <- tr.S.a[vec_idx]
      d2lS.dg2[rank_g[j], col_idx] <- S.g[vec_idx]
      tr.d2lS.dg2[rank_g[j], col_idx] <- tr.S.g[vec_idx]
      d2lS.dgda[rank_ag[j], col_idx] <- S.ga[vec_idx]
      tr.d2lS.dgda[rank_ag[j], col_idx] <- tr.S.ga[vec_idx]
      d2lh1.da2[rank_a[j], col_idx] <- h1.a[vec_idx]
      d2lh2.da2[rank_a[j], col_idx] <- h2.a[vec_idx]
      d2lh1.dg2[rank_g[j], col_idx] <- h1.g[vec_idx]
      d2lh2.dg2[rank_g[j], col_idx] <- h2.g[vec_idx]
      d2lh1.dgda[rank_ag[j], col_idx] <- h1.ga[vec_idx]
      d2lh2.dgda[rank_ag[j], col_idx] <- h2.ga[vec_idx]
      d2lcrf.da2[rank_a[j], col_idx] <- crf.a[vec_idx]
      d2lcrf.dg2[rank_g[j], col_idx] <- crf.g[vec_idx]
      d2lcrf.dgda[rank_ag[j], col_idx] <- crf.ga[vec_idx]
    }

  }

  if(any(a.idx)){

    a <- alpha_i[a_i]
    g <- gamma_i[a_i]
    s <- s_i[a_i]
    rfv <- RFV_i[a_i]
    # Survival
    S.a <- 1/rfv * ( (1+a*s)^2 - s*a - rfv/a )
    # pop hazard
    d1 <- d_1[a_i]
    h1.a <- -s*a*d1
    d2 <- d_2[a_i]
    h2.a <- -s*a*d2
    # CRF
    crf.a <-  d1*d2*a/(1+rfv) *
      ( s*rfv + (1+a*s)*(rfv/a + s*rfv)*(1 - rfv/(1+rfv)) )
    # Survival truncation
    s <- tr.s_i[a_i]
    rfv <- tr.RFV_i[a_i]
    tidx <- tr.idx[a_i]
    tr.S.a <- tidx*(1/rfv * ( (1+a*s)^2 - s*a - rfv/a ))

    for(j in type_a){
      col_idx <- idx == j
      vec_idx <- col_idx[a_i]
      d2lS.da2[rank_a[j], col_idx] <- S.a[vec_idx]
      tr.d2lS.da2[rank_a[j], col_idx] <- tr.S.a[vec_idx]
      d2lh1.da2[rank_a[j], col_idx] <- h1.a[vec_idx]
      d2lh2.da2[rank_a[j], col_idx] <- h2.a[vec_idx]
      d2lcrf.da2[rank_a[j], col_idx] <- crf.a[vec_idx]
    }

  }

  if(any(PVF.idx)){

    a <- alpha_i[PVF_i]
    g <- gamma_i[PVF_i]
    s <- s_i[PVF_i]
    rfv <- RFV_i[PVF_i]
    # Survival
    S.a <- ((-2*g*(1 - (g/(g + s))^a)/a^3 -
               2*g*(g/(g + s))^a*log(g/(g + s))/a^2 +
               g*(g/(g + s))^a*log(g/(g + s))^2/a)*(a + 1) +
              g*(1 - (g/(g + s))^a)/a^2 + g*(g/(g + s))^a*log(g/(g + s))/a)*(a + 1)
    S.g <- (((g/(g + s))^a*(1/(g + s) - g/(g + s)^2)*(g + s)/g + (g/(g + s))^a*a*(1/(g + s) - g/(g + s)^2)^2*(g + s)^2/g +
               (g/(g + s))^a*(-2/(g + s)^2 + 2*g/(g + s)^3)*(g + s) + (g/(g + s))^a*(1/(g + s) - g/(g + s)^2))*g -
              (1 - (g/(g + s))^a)/a + (g/(g + s))^a*(1/(g + s) - g/(g + s)^2)*(g + s))*g
    S.ga <- ((1 - (g/(g + s))^a)/a^2 +
               (g/(g + s))^a*log(g/(g + s))/a +
               (g/(g + s))^a*(1/(g + s) -
                                g/(g + s)^2)*(g + s)*log(g/(g + s)))*(a + 1)*g
    # pop hazard
    d1 <- d_1[PVF_i]
    h1.a <- d1*(a + 1)*log(g/(g + s))
    h1.g <- g*d1*((a + 1)*(-2/(g + s)^2 + 2*g/(g + s)^3)*(g + s) + (a + 1)*(1/(g + s) - g/(g + s)^2))
    h1.ga <- d1*((a + 1)*(1/(g + s) - g/(g + s)^2)*(g + s))
    d2 <- d_2[PVF_i]
    h2.a <- d2*(a + 1)*log(g/(g + s))
    h2.g <- g*d2*((a + 1)*(-2/(g + s)^2 + 2*g/(g + s)^3)*(g + s) + (a + 1)*(1/(g + s) - g/(g + s)^2))
    h2.ga <- d2*((a + 1)*(1/(g + s) - g/(g + s)^2)*(g + s))
    # CRF
    crf.a <- d1*d2*(
      ((2*(1 + s/g)^a*log(1 + s/g)/g +
          (a + 1)*(1 + s/g)^a*log(1 + s/g)^2/g)*(a + 1)/(1 + (a + 1)*(1 + s/g)^a/g) -
         ((1 + s/g)^a/g +
            (a + 1)*(1 + s/g)^a*log(1 + s/g)/g)^2*(a + 1)/(1 + (a + 1)*(1 + s/g)^a/g)^2 +
         ((1 + s/g)^a/g + (a + 1)*(1 + s/g)^a*log(1 + s/g)/g)/(1 + (a + 1)*(1 + s/g)^a/g))*(a + 1)
    )
    crf.g <- d1*d2*(
      ((2*(a + 1)*(1 + s/g)^a/g^3 +
          4*(a + 1)*(1 + s/g)^a*a*s/(g^4*(1 + s/g)) +
          (a + 1)*(1 + s/g)^a*a^2*s^2/(g^5*(1 + s/g)^2) -
          (a + 1)*(1 + s/g)^a*a*s^2/(g^5*(1 + s/g)^2))*g/(1 + (a + 1)*(1 + s/g)^a/g) - (-(a + 1)*(1 + s/g)^a/g^2 - (a + 1)*(1 + s/g)^a*a*s/(g^3*(1 + s/g)))^2*g/(1 + (a + 1)*(1 + s/g)^a/g)^2 + (-(a + 1)*(1 + s/g)^a/g^2 -
                                                                                                                                                                                                   (a + 1)*(1 + s/g)^a*a*s/(g^3*(1 + s/g)))/(1 + (a + 1)*(1 + s/g)^a/g))*g
    )
    crf.ga <- d1*d2*(
      ((-(1 + s/g)^a/g^2 - (a + 1)*(1 + s/g)^a*log(1 + s/g)/g^2 - (1 + s/g)^a*a*s/(g^3*(1 + s/g)) - (a + 1)*(1 + s/g)^a*log(1 + s/g)*a*s/(g^3*(1 + s/g)) - (a + 1)*(1 + s/g)^a*s/(g^3*(1 + s/g)))*g/(1 + (a + 1)*(1 + s/g)^a/g) - (-(a + 1)*(1 + s/g)^a/g^2 - (a + 1)*(1 + s/g)^a*a*s/(g^3*(1 + s/g)))*g*((1 + s/g)^a/g + (a + 1)*(1 + s/g)^a*log(1 + s/g)/g)/(1 + (a + 1)*(1 + s/g)^a/g)^2)*(a + 1)
    )
    # Survival truncation
    s <- tr.s_i[PVF_i]
    rfv <- tr.RFV_i[PVF_i]
    tidx <- tr.idx[PVF_i]
    tr.S.a <- tidx*(
      ((-2*g*(1 - (g/(g + s))^a)/a^3 -
          2*g*(g/(g + s))^a*log(g/(g + s))/a^2 +
          g*(g/(g + s))^a*log(g/(g + s))^2/a)*(a + 1) +
         g*(1 - (g/(g + s))^a)/a^2 + g*(g/(g + s))^a*log(g/(g + s))/a)*(a + 1)
    )
    tr.S.g <- tidx*(
      (((g/(g + s))^a*(1/(g + s) - g/(g + s)^2)*(g + s)/g + (g/(g + s))^a*a*(1/(g + s) - g/(g + s)^2)^2*(g + s)^2/g +
          (g/(g + s))^a*(-2/(g + s)^2 + 2*g/(g + s)^3)*(g + s) + (g/(g + s))^a*(1/(g + s) - g/(g + s)^2))*g -
         (1 - (g/(g + s))^a)/a + (g/(g + s))^a*(1/(g + s) - g/(g + s)^2)*(g + s))*g
    )
    tr.S.ga <- tidx*(
      ((1 - (g/(g + s))^a)/a^2 +
         (g/(g + s))^a*log(g/(g + s))/a +
         (g/(g + s))^a*(1/(g + s) -
                          g/(g + s)^2)*(g + s)*log(g/(g + s)))*(a + 1)*g
    )

    for(j in c(type_IG,type_PVF)){
      col_idxIG <- idx == j
      col_idx <- col_idxIG & !IG_i
      vec_idxIG <- col_idxIG[PVF_i]
      vec_idx <- col_idx[PVF_i]
      d2lS.da2[rank_a[j], col_idx] <- S.a[vec_idx]
      tr.d2lS.da2[rank_a[j], col_idx] <- tr.S.a[vec_idx]
      d2lS.dg2[rank_g[j], col_idxIG] <- S.g[vec_idxIG]
      tr.d2lS.dg2[rank_g[j], col_idxIG] <- tr.S.g[vec_idxIG]
      d2lS.dgda[rank_ag[j], col_idx] <- S.ga[vec_idx]
      tr.d2lS.dgda[rank_ag[j], col_idx] <- tr.S.ga[vec_idx]
      d2lh1.da2[rank_a[j], col_idx] <- h1.a[vec_idx]
      d2lh2.da2[rank_a[j], col_idx] <- h2.a[vec_idx]
      d2lh1.dg2[rank_g[j], col_idxIG] <- h1.g[vec_idxIG]
      d2lh2.dg2[rank_g[j], col_idxIG] <- h2.g[vec_idxIG]
      d2lh1.dgda[rank_ag[j], col_idx] <- h1.ga[vec_idx]
      d2lh2.dgda[rank_ag[j], col_idx] <- h2.ga[vec_idx]
      d2lcrf.da2[rank_a[j], col_idx] <- crf.a[vec_idx]
      d2lcrf.dg2[rank_g[j], col_idxIG] <- crf.g[vec_idxIG]
      d2lcrf.dgda[rank_ag[j], col_idx] <- crf.ga[vec_idx]
    }

  }

  if(any(g.idx)){

    a <- alpha_i[g_i]
    g <- gamma_i[g_i]
    s <- s_i[g_i]
    rfv <- RFV_i[g_i]
    # Survival
    S.g <- log(S) + s/(1+g*s) + s^2*g/S^(-2*g)
    # pop hazard
    d1 <- d_1[g_i]
    h1.g <- d1*g*s/(1+g*s)*(g*s/(1+g*s) - 1 )
    d2 <- d_2[g_i]
    h2.g <- d2*g*s/(1+g*s)*(g*s/(1+g*s) - 1 )
    # CRF
    crf.g <- d1*d2*g/(1+g)*(1 - g/(1+g))
    # Survival truncation
    s <- tr.s_i[g_i]
    rfv <- tr.RFV_i[g_i]
    tidx <- tr.idx[g_i]
    tr.S.g <- tidx*(log(tr.S) + s/(1+g*s) + s^2*g/tr.S^(-2*g))

    for(j in type_g){
      col_idx <- idx == j
      vec_idx <- col_idx[g_i]
      d2lS.dg2[rank_g[j], col_idx] <- S.g[vec_idx]
      tr.d2lS.dg2[rank_g[j], col_idx] <- tr.S.g[vec_idx]
      d2lh1.dg2[rank_g[j], col_idx] <- h1.g[vec_idx]
      d2lh2.dg2[rank_g[j], col_idx] <- h2.g[vec_idx]
      d2lcrf.dg2[rank_g[j], col_idx] <- crf.g[vec_idx]
    }

  }

  d2lS.da2 <- rowSums(d2lS.da2)
  tr.d2lS.da2 <- rowSums(tr.d2lS.da2)

  d2lS.dg2 <- rowSums(d2lS.dg2)
  tr.d2lS.dg2 <- rowSums(tr.d2lS.dg2)

  d2lS.dgda <- rowSums(d2lS.dgda)
  tr.d2lS.dgda <- rowSums(tr.d2lS.dgda)

  d2lh1.da2 <- rowSums(d2lh1.da2)

  d2lh1.dg2 <- rowSums(d2lh1.dg2)

  d2lh1.dgda <- rowSums(d2lh1.dgda)

  d2lh2.da2 <- rowSums(d2lh2.da2)

  d2lh2.dg2 <- rowSums(d2lh2.dg2)

  d2lh2.dgda <- rowSums(d2lh2.dgda)

  d2lcrf.da2 <- rowSums(d2lcrf.da2)

  d2lcrf.dg2 <- rowSums(d2lcrf.dg2)

  d2lcrf.dgda <- rowSums(d2lcrf.dgda)

  H_theta <- matrix(0, nrow = num.a + num.g, ncol = num.a + num.g)
  diag(H_theta) <- c(-tr.d2lS.da2 + d2lS.da2 + d2lh1.da2 + d2lh2.da2 + d2lcrf.da2,
                     -tr.d2lS.dg2 + d2lS.dg2 + d2lh1.dg2 + d2lh2.dg2 + d2lcrf.dg2)
  if(num.ag>0){
    rank_g2 <- rank_g[sort(c(type_ag, type_PVF))]
    rank_a2 <- rank_a[sort(c(type_ag, type_PVF))]
    i <- 1
    for( j in rank_ag[!is.na(rank_ag)] ){
      H_theta[j, j+num.a] <- (-tr.d2lS.dgda + d2lS.dgda + d2lh1.dgda + d2lh2.dgda + d2lcrf.dgda)[i]
      H_theta[j+num.a, j] <- H_theta[j, j+num.a]
      i <- i + 1
    }

  }

  return(H_theta)
}
