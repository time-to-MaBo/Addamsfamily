#' Computes Survival function given survival status of cluster comrades.
#'
#' @param object Fitted object from ph.surv
#' @param newdata_t data frame with surivavl times for which survival function should be calculated (not the conditions!)
#' @param newdata_s data frame with conditions, i.e. known survival times up to s.
#' @param Lam_ij matrix with covariate weighted cumulative baseline hazard at time t in first column and at time f condition in second column (see predict.haz)
#' @param time name or position of time variable
#' @param status name or position of event indicator
#' @param frailty.strata name or position of frailty stratum
#' @param ID name or position of ID variable
#'
#' @return predicted conditional survival functions as matrix
#' @export
#'
#' @examples
predict.S_s <- function(object, newdata_t, newdata_s = 0, Lam_ij, time=1, status=2, frailty.strata, ID){
  # wrapper for derivatives of Laplace transform for predction of population survival curves given history up to s.
  # object: ph.surv estimation object
  # newdata_s: data frame, first column is time, second status, remaining columns X. Covariates shouldn't be necassary as already included in Lam_ij.?!
  # newdata_s: data frame of history up to s
  # Lam_ij: cumulative baseline hazards on individual level (Lambda_0 * exp(x_ij beta)), no frailty!
  ## first columns evaluated at time t, second column evaluated at time min(s,t)
  # time: gives position of event time column position in newdata_t, and newdata_s (character or integer)
  # status: gives position of event status column position in newdata_t, and newdata_s (character or integer)
  # frailty.strata: gives position (or name) of frailty stratum variable
  ## if frailty is not stratified NA or NULL has to be passed!

  if(is.null(dim(Lam_ij)) | dim(Lam_ij)[2]==1 ) Lam_ij <- cbind(Lam_ij,0)
  n <- dim(Lam_ij)[1]# war vorher 2, müsste aber falsche sein?!? ----

  if(is.null(dim(newdata_t))){
    RE.stratum <- rep(1,n)
    ID <- 1:n
    #!!!! .... !!!!! AUCH BEI HAZ <-
    stop('newdata is vector, i.e. no censoring no covariates Currently not implemented!')
  }else{
    if(is.null(frailty.strata) || is.na(frailty.strata)) {
      RE.stratum <- rep(1,n)
    }else{
      RE.stratum <- newdata_t[,frailty.strata]# individual level
    }
    ID <- newdata_t[,ID]
  }
  # wie mit no frailty model umgehen?

  ag.par_i <- object$model.par[RE.stratum,] # frailty distribution parameters on individual level
  # Indices of frailty distribution on individual level:
  type_i <- object$Z.dist$type[RE.stratum]
  ag_i <- type_i %in% c('alpha_gamma', 'aNB','aNB+p', 'NB+p', 'P+p')
  a_i <- type_i == 'alpha'
  g_i <- type_i == 'gamma'
  PVF_i <- type_i %in% c('PVF', 'IG', 'cP', 'Hougaard')
  rm('type_i')

  # compute sum of cumulative baseline hazards within a cluster:
  for(i in unique(ID)){
    Lam_si <- sum(Lam_ij[i == ID,2])
    Lam_ij[i==ID,1] <- Lam_ij[i==ID,1] + Lam_si - Lam_ij[i==ID,2]
    Lam_ij[i==ID,2] <- Lam_si
  }
  # determine the order opf the derivative in order to compute conditional survival function
  derivative_order <- unlist(tapply(X = newdata_s[,status]>0, INDEX = ID, FUN = function(d){
    rep(sum(d), length(d))})) # on individual level

  S.t_s <- rep(NA, n)

  if(any(ag_i)){
    # does not work when there is more than the 15th derivative necessary! (fac = NULL)
    S.t_s[ag_i] <- ag.Laplace(s = Lam_ij[ag_i,1], a = ag.par_i[ag_i,1],
                              g = ag.par_i[ag_i,2],order = derivative_order, fac = NULL) /
      ag.Laplace(s = Lam_ij[ag_i,2], a = ag.par_i[ag_i,1],
                 g = ag.par_i[ag_i,2],order = derivative_order, fac = NULL)
    # note that if individual experiences event up to s, result will be one.
    # These individuals will be sorted out later on.
    # Also clusters with incomplete cluster history will be sorted out later on
  }
  if(any(a_i)){
    S.t_s[a_i] <- a.Laplace(s = Lam_ij[a_i,1], a = ag.par_i[a_i,1], idx = a_i,
                            order = derivative_order) /
      a.Laplace(s = Lam_ij[a_i,2], a = ag.par_i[a_i,1], idx = a_i,
                order = derivative_order)
    # note that if individual experiences event up to s, result will be one.
    # These individuals will be sorted out later on.
    # Also clusters with incomplete cluster history will be sorted out later on
  }
  if(any(g_i)){
    S.t_s[g_i] <- g.Laplace(g = ag.par_i[g_i,2], order = derivative_order, s = Lam_ij[g_i,1] ) /
      g.Laplace(g = ag.par_i[g_i,2], order = derivative_order, s = Lam_ij[g_i,2] )
  }
  if(any(PVF_i)){
    S.t_s[PVF_i] <- PVF.Laplace(a = ag.par_i[PVF_i,1], g = ag.par_i[PVF_i,2], s = Lam_ij[PVF_i,1],
                                order = derivative_order, idx = PVF_i[PVF_i]) /
      PVF.Laplace(a = ag.par_i[PVF_i,1], g = ag.par_i[PVF_i,2], s = Lam_ij[PVF_i,2],
                  order = derivative_order, idx = PVF_i[PVF_i])
  }

  # handle individual who experienced event up to t seperately (surv prob = 0):
  S.t_s[newdata_s[,status]==1] <- 0
  # handle cluster with incomplete cluster history seperately:
  ###will be done out of function.
  return(S.t_s)
}
