vz.frail <- function(alpha, gamma, pi, v.Domain, N, ag_i, g_i , a_i, PVF_i, X.tilde_C,fac,
                     d.marg, Lambda_i, frailty, n_C, ID, sum.d_i, iter, v.cond = TRUE, ind = TRUE){

  alpha_i <- alpha[frailty]
  gamma_i <- gamma[frailty]

  evx <- exp(v.Domain%x%X.tilde_C)
  s.v <- evx * Lambda_i
  s.v <- rowSums(s.v)
  evxD <- exp(matrix(v.Domain%x%rowSums(X.tilde_C*d.marg),ncol=N+1))

  prob <- stats::dbinom(0:N, size = N, prob = pi)

  del.L <- matrix(NA, nrow = n_C, ncol = 2*(N+1))
  #colnames(del.L) <- c('d^(d_i.)L', 'd^(d_i.+1)L')

  for(i in 1:(N+1)){
    if(any(ag_i)){
      del.L[ag_i,((i-1)*2+1):(i*2)] <- ag.Laplace(s = rep((s.v[((i-1)*n_C+1):(n_C*i)])[ag_i], 2),
                                                  a = rep(alpha_i[ag_i],2),g = rep(gamma_i[ag_i],2),
                                                  order = c(sum.d_i[ag_i], sum.d_i[ag_i]+1), fac = fac) *
        rep((evxD[((i-1)*n_C+1):(n_C*i)])[ag_i], 2) * prob[i]
    }
    if(any(a_i)){
      del.L[a_i,((i-1)*2+1):(i*2)] <- a.Laplace(s = rep((s.v[((i-1)*n_C+1):(n_C*i)])[a_i], 2), a = rep(alpha_i[a_i],2),
                                                order = c(sum.d_i[a_i], sum.d_i[a_i]+1), idx = rep(a_i,2)) *
        rep((evxD[((i-1)*n_C+1):(n_C*i)])[a_i], 2) * prob[i]
    }
    if(any(g_i)){
      del.L[g_i,((i-1)*2+1):(i*2)] <- g.Laplace(s = rep((s.v[((i-1)*n_C+1):(n_C*i)])[g_i], 2), g = rep(gamma_i[g_i],2),
                                                order = c(sum.d_i[g_i], sum.d_i[g_i]+1)) *
        rep((evxD[((i-1)*n_C+1):(n_C*i)])[g_i], 2) * prob[i]
    }
    if(any(PVF_i)){
      del.L[PVF_i,((i-1)*2+1):(i*2)] <- PVF.Laplace(s = rep((s.v[((i-1)*n_C+1):(n_C*i)])[PVF_i], 2), a = rep(alpha_i[PVF_i],2),
                                                    g = rep(gamma_i[PVF_i],2), order = c(sum.d_i[PVF_i], sum.d_i[PVF_i]+1),
                                                    idx = rep(PVF_i,2)) *
        rep((evxD[((i-1)*n_C+1):(n_C*i)])[PVF_i], 2) * prob[i]
    }
  }

  tmp.idx <- 1:(N+1)*2
  tmp2.idx <- tmp.idx-1
  U_bar <- rowSums(del.L[,tmp2.idx,drop=FALSE])
  z <- -rowSums(del.L[,tmp.idx,drop=FALSE])/U_bar

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
    z[tmp] <- vz.frail(alpha = alpha, gamma = gamma, v.Domain = v.Domain, Lambda_i = Lambda_i[rep(tmp,N+1),], X.tilde_C=X.tilde_C[tmp,,drop=FALSE], pi=pi,N=N,
                       d.marg=d.marg[tmp,,drop=FALSE], ag_i = rep(FALSE, n2), g_i = g_i2, a_i = a_i2, PVF_i = rep(FALSE,n2), v.cond = FALSE,
                       sum.d_i = sum.d_i[tmp], frailty = frailty[tmp], n_C = n2, ID = ID[tmp], ind = FALSE, iter = iter, fac = fac)
  }

  if(v.cond){
    for(i in 1:(N+1)){
      del.L[,tmp.idx[i]] <- del.L[,tmp.idx[i],drop=FALSE]*v.Domain[i]
    }

    v <- -rowSums(del.L[,tmp.idx,drop=FALSE])/U_bar
  }

  if(ind){
    z <- as.vector(z)
    names(z) <- levels(ID)
    z <- z[ID]
    if(v.cond){
      v <- as.vector(v)
      lnz <- log(z)
      v <- v[ID]
      names(v) <- levels(ID)
      z <- list(z = z, lnz = lnz, v = v)
    }else{
      z <- list(z = z, lnz = lnz)
    }
  }

  return(z)
}
