CI.sig <- function(par, Sigma, level = 0.95, e = FALSE, p = TRUE){

  OUT <- matrix(NA, nrow = length(par), 2)
  colnames(OUT) <- c( 'CI_lower', 'CI_upper')

  quantile <- stats::qnorm(p = 1 - (1-level)/2)

  CI_lower<- par - quantile*Sigma
  CI_upper <- par + quantile*Sigma
  if(any(e)){
    CI_lower[e] <- exp(CI_lower[e])
    CI_upper[e] <- exp(CI_upper[e])
  }
  OUT[,1] <- CI_lower
  OUT[,2] <- CI_upper

  if(p){
    p <- 2*(1 - stats::pnorm(abs(par)/Sigma))
    OUT <- cbind(OUT,p)
  }

  return(OUT)
}
