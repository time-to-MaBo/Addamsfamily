#' compares either a full (H1) and a nested (H0) model or, alternatively,
#' compares alpha_gamma vs gamma and alpha_gamma vs alpha.
#' list elements must be names in the nlm sense ($estimate, $minimum)
#'
#' @param alpha_gamma AF model (optional), alternatively H1, H0 might be specified
#' @param alpha Poisson restricted model (optional), alternatively H1, H0 might be specified
#' @param gamma gamma frailty model (optional), alternatively H1, H0 might be specified
#' @param H1 Some model (optional2), alternatively alpha_gamma, alpha, gamma might be specified
#' @param H0 Model nested in H1 model (optional2), alternatively alpha_gamma, alpha, gamma might be specified
#'
#' @return LRT, AIC, Deviance, Deviance-Test
#' @export
#'
#' @examples
LRT_AIC <- function(alpha_gamma = NULL, alpha = NULL, gamma = NULL,
                    H1 = NULL, H0 = NULL){

  if(!is.null(alpha_gamma) & !is.null(alpha) & !is.null(gamma)){

    l_alphagamma <- alpha_gamma$minimum * (-1)
    num_alphagamma <- length(alpha_gamma$estimate)
    l_alpha <- alpha$minimum * (-1)
    num_alpha <- length(alpha$estimate)
    l_gamma <- gamma$minimum * (-1)
    num_gamma <- length(gamma$estimate)

    lrt <- matrix(NA, ncol = 3, nrow = 2)
    colnames(lrt) <- c('LRT', 'df', 'p-value')
    rownames(lrt) <- c('ag vs a', 'ag vs g')
    lrt['ag vs a',] <- LRT(l_h1 = l_alphagamma, l_h0 = l_alpha,
                           df = num_alphagamma - num_alpha)
    lrt['ag vs g',] <- LRT(l_h1 = l_alphagamma, l_h0 = l_gamma,
                           df = num_alphagamma - num_gamma)

    AIC <- matrix(NA, nrow = 3)
    colnames(AIC) <- 'AIC'
    rownames(AIC) <- c('alpha-gamma', 'alpha', 'gamma')
    AIC['alpha-gamma',1] <- aic(loglihood = l_alphagamma, num_par = num_alphagamma)
    AIC['alpha',1] <- aic(loglihood = l_alpha, num_par = num_alpha)
    AIC['gamma',1] <- aic(loglihood = l_gamma, num_par = num_gamma)

  }
  else if(!is.null(H1) & !is.null(H0)){

    l_H1 <- H1$minimum * (-1)
    num_H1 <- length(H1$estimate)
    l_H0 <- H0$minimum * (-1)
    num_H0 <- length(H0$estimate)

    lrt <- LRT(l_h1 = H1, l_h0 = l_H0, df = num_H1 - num_H0)
    rownames(lrt) <- c('H1 vs H0')

    AIC <- matrix(NA, nrow = 2, ncol = 1)
    colnames(AIC) <- 'AIC'
    rownames(AIC) <- c('H1', 'H0')
    AIC['H1',1] <- aic(loglihood = l_H1, num_par = num_H1)
    AIC['H0', 1] <- aic(loglihood = l_H0, num_par = num_H0)

  }else{
    stop('Wrong input!')
  }

  print(lrt[,'p-value'])
  print(AIC)

  OUT <- list()
  OUT$AIC <- AIC
  OUT$LRT <- lrt
  return(OUT)

}
