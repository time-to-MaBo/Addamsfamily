risk.set_i <- function(z, X, beta, dt, tr.dt, m, d, n){
  # hier d?rfen keine GG-Elemente rein!

  psi_i <- as.vector(z * exp(X%*%beta))

  R <- t(psi_i*t(dt-tr.dt))

  R <- rowSums(R)
  # rows are sorted by date and stratum,
  ## columns by contribution to risk set of individual_i, i=1,...n, to specfic date

  R_i <- apply(X = m[,d==1,drop=FALSE], MARGIN = 2, FUN = function(x) t(x)%*%R)
  # rows are sorted by dates of uncesored individuals(could be simplified to unique date),

  return(R_i)

}
