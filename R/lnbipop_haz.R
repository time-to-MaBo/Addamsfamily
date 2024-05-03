lnbipop_haz <- function(alpha = NULL, gamma = NULL, s, lambda_i, n,
                        idx_ag, idx_a, idx_g, idx_PVF){

  pop_haz <- rep(NA, n)

  if(any(idx_ag)){

    pop_haz[idx_ag] <- log(lambda_i[idx_ag]) -
      log(1 + gamma[idx_ag]/alpha[idx_ag] * (exp(alpha[idx_ag]*s[idx_ag]) - 1))

  }

  if(any(idx_a)){

    pop_haz[idx_a] <- log(lambda_i[idx_a]) - alpha[idx_a]*s[idx_a]

  }

  if(any(idx_g)){

    pop_haz[idx_g] <- log(lambda_i[idx_g]) - log(1 + gamma[idx_g]*s[idx_g])

  }

  if(any(idx_PVF)){

    pop_haz[idx_PVF] <- (alpha[idx_PVF] + 1) * log( gamma[idx_PVF]/(gamma[idx_PVF] + s[idx_PVF]) ) +
      log(lambda_i[idx_PVF])

  }

  return(pop_haz)
}
