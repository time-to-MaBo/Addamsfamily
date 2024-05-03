lnCRF <- function(alpha = NULL, gamma = NULL, s, n,
                  idx_ag, idx_a, idx_g, idx_PVF){

  crf <- rep(NA, n)

  if(any(idx_ag)){

    crf[idx_ag] <- log(1 + gamma[idx_ag] * exp(alpha[idx_ag]*s[idx_ag]))

  }

  if(any(idx_a)){

    crf[idx_a] <- log(1 + alpha[idx_a] * exp(alpha[idx_a]*s[idx_a]))

  }

  if(any(idx_g)){

    crf[idx_g] <- log(1 + gamma[idx_g])

  }

  if(any(idx_PVF)){

    crf[idx_PVF] <- log(
      (alpha[idx_PVF] + 1)/gamma[idx_PVF] * (1+s[idx_PVF]/gamma[idx_PVF])^alpha[idx_PVF] + 1
    )

  }

  return(crf)
}
