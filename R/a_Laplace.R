a.Laplace <- function(a,s,idx,order){
  # Derivative of Laplace. You need to * (-^1)^del.n for loglihood(densities)
  L <- exp(lnS(alpha = a, gamma = NULL, s = s, n = sum(idx), idx_ag = FALSE, idx_a = idx[idx],
               idx_g = FALSE, idx_PVF = FALSE))
  if(any(order > 0)){ # new 25.02.2024 for predction of new observation where no event might have occured.
    s2 <- s[order > 0]
    order2 <- order[order>0]
    a <- a[order > 0]
    Bell <- rep(NA, length(s2))
    for(i in 1:length(s2)){
      del.n <- order2[i]
      a1 <- a[i]
      series <- exp(-a1*s2[i]) * (-a1)^(1:del.n)/a1
      M <- matrix(0, nrow = del.n, ncol = del.n)
      diag(M) <- series[1]
      if(del.n > 1){
        M[row(M)-1 == col(M)] <- -1
        for(r in 1:(del.n-1)){
          for(c in (r+1):del.n){
            M[r,c] <- choose(n = del.n-r, k = c-r) * series[c-r+1]
          }
        }
      }
      Bell[i] <- det(M)
    }
    L[order > 0] <- L[order > 0] * Bell
  }


  return(L)

}
