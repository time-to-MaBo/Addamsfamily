ln.d3Laplace <- function(alpha, gamma, s, n, idx_ag, idx_a, idx_g, idx_PVF){

  d3La <- rep(NA, n)

  if(any(idx_ag)){
    a <- alpha[idx_ag]
    g <- gamma[idx_ag]
    s_i <- s[idx_ag]
    innr_adptd <- -exp(-a*s_i)
    S <- exp(lnS(alpha=a, gamma=g, s=s_i, n=sum(idx_ag),
                 idx_ag=TRUE, idx_a=FALSE, idx_g=FALSE, idx_PVF=FALSE))
    d3La[idx_ag] <- (1-a+g)*innr_adptd*S^(1-a+g)*(
      (1-2*a+2*g)*S^(-2*a+2*g)*innr_adptd^2 -
        S^(-a+g)*2*a*innr_adptd - a*(
          S^(-a+g) * innr_adptd -
            a/(1-a+g)
        )
    )
  }

  if(any(idx_a)){
    a <- alpha[idx_a]
    s_i <- s[idx_a]
    innr_adptd <- -exp(-a*s_i)
    S <- exp(lnS(alpha=a, gamma=NULL, s=s_i, n=sum(idx_a),
                 idx_ag=FALSE, idx_a=TRUE, idx_g=FALSE, idx_PVF=FALSE))
    d3La[idx_a] <- S*innr_adptd*( innr_adptd^2 - 3*a*innr_adptd + a^2 )
  }

  if(any(idx_g)){
    g <- gamma[idx_g]
    s_i <- s[idx_g]
    d3La[idx_g] <- -(1+g)*(1+2*g)*(1+g*s_i)^(-(1+3*g)/g)
  }

  if(any(idx_PVF)){
    a <- alpha[idx_PVF]
    g <- gamma[idx_PVF]
    s_i <- s[idx_PVF]
    innr_adpt <- g/(g+s_i)
    S <- exp(lnS(alpha=a, gamma=g, s=s_i, n=sum(idx_PVF),
                 idx_ag=FALSE, idx_a=FALSE, idx_g=FALSE, idx_PVF=TRUE))
    d3La[idx_PVF] <- -S*innr_adpt^(a+3) * (
      innr_adpt^(2*a) + 3*(a+1)/g * innr_adpt^(a) +
        (a+1)*(a+2)/g^2
    )
  }

  ln.d3La <- log(-d3La)
  return(ln.d3La)

}
