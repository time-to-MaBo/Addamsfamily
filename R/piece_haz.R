piece.haz <-function(stratum, partial.lambda, breaks, y, baseline){
  # Hazard Funktion insbesondere für out of sample estimates.
  # only Breslow or piecewise hazard individuals are allowed in here!
  # neu am 22.04
  # y: time
  # partial.lambda list with names of strata
  # breaks: time breaks of piecewise or partial hazard
  # stratum of individual wrt hazard
  # baseline specifies wether hazard is non-parametric or piecewise constant
  warning('I dont understand the function anymore. Carefully check!')
  Breslow <- 'Breslow' %in% baseline
  dt <- lapply(X = breaks, FUN = diff)

  Lam <- lapply(X = 1:length(y), FUN = function(i){
    strat_i <- stratum[i]
    b_i <- breaks[[strat_i]]
    dt_i <- dt[[strat_i]]
    lam <- partial.lambda[[strat_i]]
    if(Breslow){
      idx <- y[i] >= b_i[-1]
      dt_i[idx] <- 1
    }else{
      idx <- y[i] >= b_i[-length(b_i)] # is zero in b_i? I think so. Maybe delete fist b_i value and do not delete last entry in next row. No! does not work for peicewise
      t_k <- sum(idx)
      dt_i[t_k] <- y[i] - b_i[t_k]
    }
    dt_i[!idx] <- 0
    Lam_0 <- sum(dt_i*lam)
    return(Lam_0)
  })

  return(unlist(Lam))
}
