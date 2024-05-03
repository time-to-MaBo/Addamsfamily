GG_tree <- function(sigma,eta,nu,Q,GG.idx,edge=1e-1,names.strata=1:Q){

  # Prepare Output
  temp <- replicate(list(), n = sum(GG.idx))
  OUT <- replicate(list(), n = Q)
  names(OUT) <- names.strata

  Cox_par <- matrix(NA, nrow = Q, ncol = 3)
  colnames(Cox_par) <- c('sigma', 'eta', 'nu')
  rownames(Cox_par) <- names.strata
  Cox_par[GG.idx,] <- cbind(sigma,eta,nu)

  if(any(GG.idx)){
    # get it to Cox et al. (2007) parameterization

    lambda <- nu
    #beta <- eta
    # sigma <- sigma

    ## which case? Cox et al. (2007), p.4356.
    ### Decision is made when parameters are very close to special case
    Weibull <- abs(lambda - 1) < edge
    Lognormal <- abs(lambda) < edge
    Half_normal <- (abs(lambda - sqrt(2)) < edge) & (abs(sigma - 1/sqrt(2)) < edge)
    Inv_Weibull <- (abs(lambda - -1)) < edge
    Inv_Ammag <-  lambda < 0 & ( (abs(lambda - -1/sigma) < edge) | (abs(1/lambda - -sigma) < edge) ) &
      !Inv_Weibull & !Lognormal
    Gamma <- lambda >0 & (abs(lambda-sigma) < edge) & !Weibull & !Lognormal
    Ammag <- lambda >0 & ( (abs(lambda - 1/sigma) < edge) | (abs(1/lambda - sigma) < edge) ) &
      (!Gamma & !Lognormal & !Half_normal)
    Inv_Gamma <- lambda < 0 & (abs(lambda - -sigma) < edge) & (!Lognormal & !Inv_Weibull & !Inv_Ammag)
    exponential <- abs(sigma-1) < edge & abs(nu-1) < edge
    Inv_exponential <- abs(sigma-1) < edge & abs(nu+1) < edge
    Weibull[exponential] <- FALSE
    GG3 <- !Gamma & !Ammag & !Lognormal & !Inv_Gamma & !Inv_Ammag & !Inv_Weibull & !Half_normal & !Weibull & !exponential

    temp[Gamma] <- 'T_0 ~ Gamma(shape=sigma^-2, rate=sigma^-2 exp{-eta}) (approx.)'
    temp[Ammag] <- 'T_0 ~ Ammag (approx.)'
    temp[Weibull] <- 'T_0 ~ Weibull(shape=1/sigma, scale=exp(eta)) (approx.)'
    temp[Lognormal] <- 'T_0 ~ Lognormal(mu=eta, sd=sigma) (approx.)'
    temp[Inv_Gamma] <- 'T_0 ~ Gamma^-1(shape=sigma^-2, scale=sigma^-2 / exp{-eta}) (approx.)'
    temp[Inv_Ammag] <- 'T_0 ~ Ammag^-1 (approx.)'
    temp[Inv_Weibull] <- 'T_0 ~ Weibull^-1(shape=1/sigma, scale=exp{eta}) (approx.)'
    temp[Half_normal] <- 'T_0 ~ HN(scale = exp(eta)) (approx.)'
    temp[exponential] <- 'T_0 ~ Exp(rate = exp(-eta)) (approx.)'
    temp[Inv_exponential] <- 'T ~ Exp^-1(scale = exp(eta)) (approx.)'
    temp[GG3] <- 'T_0 ~ GG(sigma,eta,nu)'

    OUT[GG.idx] <- temp
  }

  OUT[!GG.idx] <- list('Baseline hazard is (piecewise) constant/exponential')
  OUT$Cox_par <- Cox_par

  return(OUT)

}
