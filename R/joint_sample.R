#' purpose of function: Jointly draw parameters from asymptotic distribution of ML estimates.
#' Samples can be used to empirically construct CI's for the RFV or th survival function.
#' Note that CIs wont be supplied by function, only samples will be drawn.
#' Consider to set seed beforehand!
#'
#' @param n_samp number of samples to be drawn from multivariate normal
#' @param marginal_par all parameters from marginal$par, i.e. from optimization, on optimzed scale
#' @param cond_par same for conditional model parameters (beta, hazard) (irrelevant for current status data)
#' @param H Hessian. Order of parameters has to be in line with c(marginal_par,cond_par).  no further parameters allowed!
#' @param exclude specify positions of parameters that should be excluded from being samples
#' @param H.singular Handling of hessian if it is singular.
#'
#' @return list of samples and input variance of multivariate normal form which samples were drawn.
#' @export
#'
#' @examples
joint_sample <- function(n_samp=1e5,marginal_par,cond_par=NULL,H,exclude = FALSE, H.singular = 'standard'){
  # purpose of function: Jointly draw parameters from asymptotic distribution of ML estimates.
  # Samples can be used to empirically construct CI's for the RFV or th survival function.
  # Note that CIs wont be supplied by function, only samples will be drawn.
  ## n_samp: number of samples to be drawn from multivariate normal
  ## maringal_par: all parameters from marginal$par, i.e. from optimization, on optimzed scale
  ## cond_par: same for conditional model parameters (beta, hazard)
  ### irrelaevant for current status data as conditional model not subject to optimization
  ### and all parameters are contained in marginal
  ## H: Hessian. Order of parameters has to be in line with c(marginal_par,cond_par)
  ### no further parameters allowed!
  ## H.singular: Handling of hessian if it is singular.
  ### copy paste from maxll

  par <- c(marginal_par,cond_par)
  var_names <- colnames(H)

  ### copy paste part from maxll:
  Sigma2 <- tryCatch(expr = solve(-H),
                     error = function(e) {warning('Hessian is singular. Handling can be specified by H.singular!');
                       diag(Inf, nrow = dim(H)[1])} )
  H.check <- all(is.infinite(diag(Sigma2)))
  H.check2 <- diag(Sigma2<0)
  if((H.check|any(H.check2,na.rm=TRUE)) & H.singular != 'standard'){
    if(any(H.check2)){
      tmp <- which(H.check2)
      warning(paste(c('Sigma^2 had negative Variances. This was the case for parmeters',
                      paste(tmp,collapse=','),'according to order in H.'),collapse=' '))
    }else if(H.check){
      tmp <- which(abs(diag(H)) < 1e-7)
      warning(paste(c('The following diagonal elements in H are very close to zero:',
                      paste(tmp,collapse=','),'. Might be the reason for H being singular.',
                      'Often happens if parameter is close to boundary!'),collapse=' '))
    }

    if(H.singular == 'gcholesky'){
      Sigma2 <- solve(bdsmatrix::gchol(-H))
    }else if(H.singular == 'ginverse'){
      Sigma2 <- MASS::ginv(-H)
    }else{
      stop(paste('H.singular can take the values standard, gcholesky and ginverse.',
                 'Your incorrect choice was',H.singular))
    }
  }
  ### end of copy paste part

  if(all(!is.na(exclude))){
    Sigma2 <- Sigma2[,-exclude]
    Sigma2 <- Sigma2[-exclude,]
    par <- par[-exclude]
    var_names <- var_names[-exclude]
  }

  idx <- rowall(t(Sigma2==0))
  if(any(idx)){
    warning('There was a woe/column which consists of zeros in covariance matrix.
            Usually occurs if a parameters i close to parameter boundary.
            The parameter was removed and in each sample a zero was imputed.')
    for(i in which(idx)){
      Sigma2 <- Sigma2[-i,]
      Sigma2 <- Sigma2[,-i]
      par <- par[-i]
    }
  }

  samples <- mnorm::rmnorm(n = n_samp, mean = par, sigma = Sigma2)

  if(any(idx)){
    for(i in which(idx)){
      if(i==1){
        samples <- cbind(0,samples)
      }else if(i<dim(Sigma2)[2]){
        samples <- cbind(samples[,1:(i-1)],0, samples[,i:dim(samples)[2]])
      }else{
        samples <- cbind(samples[,1:(i-1)],0)
      }
    }
  }

  colnames(samples) <- var_names
  rownames(samples) <- paste('draw no.', 1:n_samp)

  return(list(samples,Sigma2))

}
