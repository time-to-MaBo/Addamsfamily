PVF.Laplace <- function(a,g,s,order,idx){
  # Derivative of Laplace. You need to * (-^1)^del.n for loglihood(densities)
  L <- exp(lnS(alpha = a, gamma = g, s = s, n = sum(idx), idx_ag = FALSE, idx_a = FALSE,
               idx_g = FALSE, idx_PVF = idx[idx]))

  if(any(order>0)){# new 25.02.2024. For prediction nof out of sample data where no event might have occured
    s2 <- s[order > 0]
    a <- a[order > 0]
    g <- g[order > 0]
    order2 <- order[order>0]
    Bell <- rep(NA, length(s2))
    for(i in 1:length(s2)){
      del.n <- order2[i]
      a1 <- a[i]
      g1 <- g[i]
      s1 <- s2[i]
      ph <- sapply(X = 1:del.n, FUN = function(x){
        orthopolynom::pochhammer(z=-a1-x+1, n=x)
      }, simplify = TRUE)
      series <- g1/a1  * ph/g1^(1:del.n) * ((g1+s1)/g1)^(-a1 - 1:del.n)
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
