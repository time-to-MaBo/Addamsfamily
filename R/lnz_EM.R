lnz_EM <- function(alpha, gamma, ag_i, g_i , a_i, PVF_i, s, d, frailty, n, ID, ind = TRUE, iter = NULL){

  alpha_i <- alpha[frailty]
  gamma_i <- gamma[frailty]

  # Survivors
  z <- lnbipop_haz(alpha = alpha_i, gamma = gamma_i, s = s,
                   lambda_i = rep(1, n), n = n, idx_ag = ag_i,
                   idx_a = a_i, idx_g = g_i, idx_PVF = PVF_i)

  # single event
  idx <- rowSums(d) >= 1
  z[idx] <- z[idx] +
    lnCRF(alpha = alpha_i[idx], gamma = gamma_i[idx], s = s[idx],
          n = sum(idx), idx_ag = ag_i[idx],
          idx_a = a_i[idx], idx_g = g_i[idx], idx_PVF = PVF_i[idx])

  # both dead
  idx <- rowSums(d) == 2
  z[idx] <- ln.d3Laplace(alpha = alpha_i[idx], gamma = gamma_i[idx], s = s[idx],
                         n = sum(idx), idx_ag = ag_i[idx],
                         idx_a = a_i[idx], idx_g = g_i[idx], idx_PVF = PVF_i[idx]) -
    (
      z[idx] +

        lnbipop_haz(alpha = alpha_i[idx], gamma = gamma_i[idx], s = s[idx],
                    lambda_i = rep(1, sum(idx)), n = sum(idx), idx_ag = ag_i[idx],
                    idx_a = a_i[idx], idx_g = g_i[idx], idx_PVF = PVF_i[idx]) +

        lnS(alpha = alpha_i[idx], gamma = gamma_i[idx], s = s[idx],
            n = sum(idx), idx_ag = ag_i[idx],
            idx_a = a_i[idx], idx_g = g_i[idx], idx_PVF = PVF_i[idx])
    )

  tmp <- is.nan(z)
  if(any(tmp)){
    g.app <- abs(alpha) < abs(gamma - alpha)
    P.app <- !g.app
    g_i2 <- g.app[frailty[tmp]]
    a_i2 <- P.app[frailty[tmp]]
    n2 <- sum(tmp)
    warning('Individuals with frailty stratum ', paste(unique(frailty[tmp]),sep = ','), ' had invalid parameter constellation
          and hence, either a poisson model or gamma was chosen to calculate frailties in one iteration.
          If those individuals are not from the alpha_gamma model this approach is total non-sense.
          The should not happen however. This happened in global iteration ', iter)
    z[tmp] <- lnz_EM(alpha, gamma, ag_i = rep(FALSE, n2), g_i = g_i2, a_i = a_i2, PVF_i = rep(FALSE,n2),
                     s = s[tmp], d = d[tmp,,drop = FALSE], frailty = frailty[tmp], n = n2, ID = ID[tmp], ind = FALSE, iter = iter)
  }

  if(ind){
    lnz <- as.vector(z)
    names(lnz) <- levels(ID)
    lnz <- lnz[ID]
    z <- exp(lnz)
    z <- list(z = z, lnz = lnz)
  }

  return(z)
}
