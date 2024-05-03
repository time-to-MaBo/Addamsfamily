z.cs <- function(alpha,gamma,idx,s,ag_i,a_i,g_i,PVF_i,n_C,sign.info,ID,ind = TRUE, iter = NULL){

  alpha_i <- alpha[idx]
  gamma_i <- gamma[idx]

  del <- del2 <- condi <- s*NA
  condi[,1] <- exp(lnS(alpha = alpha_i, gamma = gamma_i, s = s[,1], n = n_C,idx_ag = ag_i,
                       idx_a = a_i, idx_g = g_i, idx_PVF = PVF_i)*(s[,1]>0))
  del[,1] <- -exp(lnbipop_haz(alpha = alpha_i, gamma = gamma_i, s = s[,1], lambda_i = rep(1,n_C),
                              n = n_C,idx_ag = ag_i, idx_a = a_i, idx_g = g_i, idx_PVF = PVF_i)) *
    condi[,1]
  del2[,1] <- exp(lnCRF(alpha = alpha_i, gamma = gamma_i, s = s[,1],
                        n = n_C,idx_ag = ag_i, idx_a = a_i, idx_g = g_i, idx_PVF = PVF_i))

  for(k in 2:(dim(s)[2])){
    condi[,k] <- exp(lnS(alpha = alpha_i, gamma = gamma_i, s = s[,k] + s[,1], n = n_C,idx_ag = ag_i,
                         idx_a = a_i, idx_g = g_i, idx_PVF = PVF_i))*(s[,k]>0)
    del[,k] <- -exp(lnbipop_haz(alpha = alpha_i, gamma = gamma_i, s = s[,k] + s[,1], n = n_C,
                                lambda_i = rep(1,n_C), idx_ag = ag_i, idx_a = a_i, idx_g = g_i, idx_PVF = PVF_i)) *
      condi[,k]
    del2[,k] <- exp(lnCRF(alpha = alpha_i, gamma = gamma_i, s = s[,k] + s[,1], n = n_C,idx_ag = ag_i,
                          idx_a = a_i, idx_g = g_i, idx_PVF = PVF_i))
  }
  del2 <- del2*del^2
  del2[condi>0] <- del2[condi>0]/condi[condi>0]
  del2[condi==0] <- 0
  del <- del*sign.info
  condi <- condi*sign.info
  del2 <- del2*sign.info

  z <- -rowSums(del)/rowSums(condi)
  z2 <- rowSums(del2)/rowSums(condi)
  tmp <- z < 0
  if(any(tmp)){
    g.app <- abs(alpha) < abs(gamma - alpha)
    P.app <- !g.app
    g_i2 <- g.app[idx[tmp]]
    a_i2 <- P.app[idx[tmp]]
    n2 <- sum(tmp)
    warning('Individuals with frailty stratum ', paste(unique(idx[tmp]),sep = ','), ' had invalid parameter constellation
          and hence, either a poisson model or gamma was chosen to calculate frailties in one iteration.
          If those individuals are not from the alpha_gamma model this approach is total non-sense.
          The should not happen however. This happened in global iteration ', iter, 'Not well tested for cs data yet anyway!!!')
    z.tmp <- z.cs(alpha = alpha[P.app], gamma = gamma[g.app], ag_i = rep(FALSE, n2), g_i = g_i2, a_i = a_i2, PVF_i = rep(FALSE,n2),
                  s = s[tmp,,drop=FALSE], idx = idx[tmp], n = n2, ID = ID[tmp], ind = FALSE, iter = iter,
                  sign.info[tmp,,drop=FALSE])
    z[tmp] <- z.tmp$z
    z2[tmp] <- z.tmp$z2
  }

  if(ind){
    names(z) <- names(z2)<-  levels(ID)
    z <- z[ID]
    z2 <- z2[ID]
    attr(z,'z2') <- z2
  }

  z <- list(z=z,z2=z2)
  return(z)
}
