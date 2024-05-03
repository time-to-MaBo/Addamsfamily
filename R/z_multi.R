z.multi <- function(alpha, gamma, ag_i, g_i , a_i, PVF_i, s, frailty, n, ID, sum.d_i, iter, fac, ind = TRUE){

  alpha_i <- alpha[frailty]
  gamma_i <- gamma[frailty]

  del.L <- matrix(NA, nrow = n, ncol = 2)
  colnames(del.L) <- c('d^(d_i.)L', 'd^(d_i.+1)L')

  if(any(ag_i)){
    del.L[ag_i,] <- ag.Laplace(s = rep(s[ag_i], 2), a = rep(alpha_i[ag_i],2), fac = fac,
                               g = rep(gamma_i[ag_i],2), order = c(sum.d_i[ag_i], sum.d_i[ag_i]+1))
  }
  if(any(a_i)){
    del.L[a_i,] <- a.Laplace(s = rep(s[a_i], 2), a = rep(alpha_i[a_i],2),
                             order = c(sum.d_i[a_i], sum.d_i[a_i]+1), idx = rep(a_i,2))
  }
  if(any(g_i)){
    del.L[g_i,] <- g.Laplace(s = rep(s[g_i], 2), g = rep(gamma_i[g_i],2),
                             order = c(sum.d_i[g_i], sum.d_i[g_i]+1))
  }
  if(any(PVF_i)){
    del.L[PVF_i,] <- PVF.Laplace(s = rep(s[PVF_i], 2), a = rep(alpha_i[PVF_i],2),
                                 g = rep(gamma_i[PVF_i],2), order = c(sum.d_i[PVF_i], sum.d_i[PVF_i]+1),
                                 idx = rep(PVF_i,2))
  }

  z <- -del.L[,2]/del.L[,1]

  tmp <- z < 0
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
    z[tmp] <- z.multi(alpha = alpha, gamma = gamma, ag_i = rep(FALSE, n2), g_i = g_i2, a_i = a_i2, PVF_i = rep(FALSE,n2), fac = fac,
                      s = s[tmp], sum.d_i = sum.d_i[tmp], frailty = frailty[tmp], n = n2, ID = ID[tmp], ind = FALSE, iter = iter)
  }

  if(ind){
    z <- as.vector(z)
    names(z) <- levels(ID)
    z <- z[ID]
    lnz <- log(z)
    z <- list(z = z, lnz = lnz)
  }

  return(z)
}
