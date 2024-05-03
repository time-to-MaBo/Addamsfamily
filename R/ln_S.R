lnS <- function(alpha = NULL, gamma = NULL, s, n,
                idx_ag, idx_a, idx_g, idx_PVF){

  pi_00 <- rep(NA, n)

  if(any(idx_ag)){

    pi_00[idx_ag] <- 1/(alpha[idx_ag]  - gamma[idx_ag] ) *
      log( (1-gamma[idx_ag]/alpha[idx_ag]) * exp(-alpha[idx_ag]*s[idx_ag]) +
             gamma[idx_ag]/alpha[idx_ag] )

  }

  if(any(idx_a)){

    pi_00[idx_a] <- 1/alpha[idx_a] * (exp(-alpha[idx_a]*s[idx_a]) - 1)

  }

  if(any(idx_g)){

    pi_00[idx_g] <- log(1 + gamma[idx_g]*s[idx_g])*(-1/gamma[idx_g])

  }

  if(any(idx_PVF)){

    pi_00[idx_PVF] <- -gamma[idx_PVF]/alpha[idx_PVF] *
      (1 - (gamma[idx_PVF]/(gamma[idx_PVF] + s[idx_PVF]))^alpha[idx_PVF])

  }

  return(pi_00)
}
