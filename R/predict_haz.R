#' wrapper for piece.haz und GG.haz for predction of hazard rates (Lam_o*exp(xb))
#'
#' @param object ph.surv estimation object
#' @param newdata newdata: dataframe, first column is time, remaining colmns X
#' @param time gives position of time column or its name
#' @param haz.strata gives position of stratum variable of baseline hazard in the same way as time if baseline hazard is not stratified NA has to be passed!
#'
#' @return matrix of cumlative baseline hazards
#' @export
#'
#' @examples
predict.haz <- function(object, newdata, time = 1, haz.strata){
  # wrapper for piece.haz und GG.haz for predction of hazard rates (Lam_o*exp(xb))
  # ph.surv estimation object
  # newdata: dataframe, first column is time, remaining colmns X
  # time: gives position of time column or its name
  # haz.strata: gives position of stratum variable of baseline hazard in the same way as time
  ## if baseline hazard is not stratified NA has to be passed!

  if(is.null(dim(newdata))){
    y <- newdata
    X <- matrix(y*0,ncol = 1)
    beta <- matrix(0)
    n <- length(y)
    # ....
    stop('newdata is a vector, i.e. no censoring no coavriates. Currently not implemented!')
  }else{
    y <- newdata[,time]
    n <- length(y)
    if(all(colnames(object$plot.predict$X)=='zero')){
      X <- matrix(0)
    }else{
      X <- newdata[,colnames(object$plot.predict$X),drop = FALSE]
      for(k in dim(X)[2]){
        if(is.factor(X[,k])){
          X[,k] <- as.numeric(X[,k])-1
        }
      }
      mu <-attr(object$plot.predict$X,'mean')
      sd <- attr(object$plot.predict$X,'sd')
      X <- scale(x = X, center = mu, scale = sd)
    }
    if(is.null(haz.strata) || is.na(haz.strata)){
      stratum <- rep(1,n)
    }else{
      stratum <- newdata[,haz.strata]# individual level
    }
    strata <- object$plot.predict$strata # unique strata
    beta <- matrix(object$plot.predict$beta)
  }

  baseline <- object$baseline[[1]]

  partial_i <- baseline[stratum] %in% c('piecewise','Breslow')
  GG_i <- !partial_i
  if(any(partial_i) & any(GG_i)) stop('Mix of GG and piecewise hazard not implemented in predict.haz!')
  #GG.stratum <- as.numeric(stratum)

  breaks <- object$plot.predict$breaks
  partial.lambda <- object$plot.predict$partial.lambda
  GG.par <- object$T.dist$Cox_par

  Lambda <- rep(NA, n)

  if(any(partial_i)){
    Lambda[partial_i] <- piece.haz(stratum=stratum, partial.lambda=partial.lambda,
                                   breaks=breaks, y=y, baseline=baseline)
  }
  if(any(GG_i)){
    Lambda[GG_i] <- GG.haz(s=GG.par[,1],e=GG.par[,2],nu=GG.par[,3],y=y,str=stratum,n=n)[,1]
  }

  return(Lambda*exp(as.vector(X%*%beta)))
}
