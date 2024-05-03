g.Laplace <- function(g,s,order){
  # Derivative of Laplace. You need to * (-^1)^del.n for loglihood(densities)
  pochh <- mapply(FUN = function(g, o){
    orthopolynom::pochhammer(z = -1/g - o + 1, n = o)
  },
  g = g, o = order, SIMPLIFY = TRUE)

  L <- g^order * pochh * (g*s+1)^(-1/g - order)

  return(L)
}
