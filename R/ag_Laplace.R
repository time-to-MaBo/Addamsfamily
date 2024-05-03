ag.Laplace <- function(s,a,g,order,fac){
  # original version

  order.idx <- order>15
  L <- s*NA

  if(any(!order.idx)){
    L[!order.idx] <- mapply(FUN = function(alpha,gamma,s,n) eval(parse(text = ag.L[[n+1]])),
                            n = order[!order.idx], alpha = a[!order.idx], gamma = g[!order.idx], s = s[!order.idx],
                            SIMPLIFY = TRUE)
  }
  if(any(order.idx)){
    g <- g[order.idx]
    a <- a[order.idx]
    s <- s[order.idx]
    order <- order[order.idx]
    n_C <- sum(order.idx)
    del.Lap <- rep(NA,n_C)

    rfv <- g*exp(a*s)
    cumRFV1 <- 1 + (rfv-g)/a
    S <- ( (1-g/a)*exp(-a*s) + g/a )^(1/(a-g))

    for(i in 1:n_C){
      tmp <- fac[[order[i]]]
      del.Lap[i] <- sum(
        rfv[i]^tmp[,'idx2']/cumRFV1[i]^tmp[,'idx1'] * a[i]^tmp[,'exponent'] * tmp[,'factor']
      )
    }

    L[order.idx] <- del.Lap*S
  }

  return(L)

}
