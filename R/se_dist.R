#' Computes confidence intervals for for RC frailty values and conditional hazard ratios
#'
#' @param zstar
#' @param model.par
#' @param mu
#' @param Sigma
#' @param type
#' @param num.a
#' @param num.g
#' @param Kstar
#' @param mean.stratum
#' @param RE.stratum
#' @param CI
#' @param dist.par
#'
#' @return
#' @export
#'
#' @examples
se.dist <- function(zstar=0:10,model.par,mu=1,Sigma,type='alpha_gamma',
                    num.a=1,num.g=1,Kstar=0,mean.stratum=FALSE,
                    RE.stratum=row.names(model.par),CI=0.95,dist.par){
  # nur einzelner EW erlaubt, muss direkt nach alpha und gamma Eintr?gen in Sigma kommen!
  # FKT funktioniert nur mit Kstar = 1, kann man ?ber mean.stratum de facto auf 0 setzen

  J <- length(type)
  if(any(type!='alpha_gamma' | any(model.par[,'alpha']>=0))){
    stop('only aNB+p implemented yet (estimated via alpha_gamma)!')
  }

  if(is.character(mean.stratum)){
    mean.stratum <- which(RE.stratum == mean.stratum)
  }
  names(mu) <- RE.stratum
  Sigma <- Sigma[1:(num.a+num.g+Kstar),1:(num.a+num.g+Kstar)]

  z <- lnz <- cbind(zstar*NA,zstar*NA)
  colnames(z) <- colnames(lnz) <- RE.stratum
  row.names(z) <- row.names(lnz) <- paste('RC=',zstar+1,sep='')

  h <- matrix(NA,nrow=4,ncol=3)
  colnames(h) <- c('alpha','gamma','mu')
  rownames(h) <- c('phi','nu','p','pstar')

  hz <- matrix(NA,nrow=length(zstar),ncol=3)
  colnames(hz) <- c('alpha','gamma','mu')
  rownames(hz) <- paste('z*=',zstar,sep='')
  hz <- replicate(list(hz),n=J)
  hP <- hcdf <- hz
  SE.z <- matrix(NA, nrow = length(zstar), ncol = J)
  rownames(SE.z) <- paste('RC=',zstar+1,sep='')
  SE.P <- SE.dP <- SE.cdf <- SE.z
  SE.par <- matrix(NA, nrow = 4, ncol = J)
  rownames(SE.par) <- c('phi','nu','p','pstar')
  colnames(SE.par) <- colnames(SE.z) <- RE.stratum
  lnP <- P <- cdf <- lncdf <- matrix(NA,nrow=length(zstar),ncol=J)

  phistar <- abs(model.par[,'alpha'])*mu
  nu <- dist.par[,'nu']
  p <- dist.par[,'p']
  pstar <- p*mu

  CI.par <- CI.z <- CI.WHR <- CI.P <- CI.cdf <- list()

  for(i in 1:J){
    a <- model.par[i,'alpha']
    g <- model.par[i,'gamma']
    m <- mu[i]
    n <- nu[i]
    pi <- p[i]
    S <- Sigma[c(0,num.a,num.a+num.g+Kstar-i)+i,c(0,num.a,num.a+num.g+Kstar-i)+i]

    z[,i] <- abs(a)*m*zstar+pstar[i]
    lnz[,i] <- log(z[,i])

    # z CI
    hz[[i]][,1] <- (-zstar*m - m/(g-a) - a*m/(g-a)^2)/z[,i]
    hz[[i]][,2] <- (a/(g-a)^2 * m)*g/z[,i]
    hz[[i]][,3] <- (i%in%mean.stratum)+0

    SE.z[,i] <- sqrt(diag( hz[[i]] %*% S %*% t(hz[[i]]) ))

    CI.z[[i]] <- exp(cbind(
      lnz[,i] - stats::qnorm(p = 1-(1-CI)/2, mean = 0, sd = 1)*SE.z[,i],
      lnz[,i],
      lnz[,i] + stats::qnorm(p = 1-(1-CI)/2, mean = 0, sd = 1)*SE.z[,i]))

    # P(z*)-CI
    tmp1 <- digamma(zstar + n) - digamma(n) + log(pi)
    tmp2 <- n/pi- zstar/(1-pi)
    lnP[,i] <- stats::dnbinom(x = zstar, size = n, prob = pi, log = TRUE)
    hP[[i]][,1] <- tmp1/(g-a)^2 + tmp2*(a/(g-a) - 1)/(g-a)
    hP[[i]][,2] <- (-tmp1/(g-a)^2 + tmp2*a/(g-a)^2)*g
    hP[[i]][,3] <- 0
    hP[[i]] <- hP[[i]]/lnP[,i]

    SE.P[,i] <- sqrt(diag( hP[[i]] %*% S %*% t(hP[[i]]) ))
    CI.P[[i]] <- exp(-exp(cbind(
      log(-lnP[,i]) + stats::qnorm(p = 1-(1-CI)/2, mean = 0, sd = 1)*SE.P[,i],
      log(-lnP[,i]),
      log(-lnP[,i]) - stats::qnorm(p = 1-(1-CI)/2, mean = 0, sd = 1)*SE.P[,i])))

    P[,i] <- exp(lnP[,i])
    cdf[,i] <- cumsum(P[,i])
    lncdf[,i] <- log(cdf[,i])
    hcdf[[i]] <- apply(X = hP[[i]],MARGIN = 2, FUN = cumsum)*lnP[,i]*P[,i]/(lncdf[,i]*cdf[,i])
    SE.cdf[,i] <- sqrt(diag( hP[[i]] %*% S %*% t(hP[[i]]) ))
    CI.cdf[[i]] <- exp(-exp(cbind(
      log(-lncdf[,i]) + qnorm(p = 1-(1-CI)/2, mean = 0, sd = 1)*SE.cdf[,i],
      log(-lncdf[,i]),
      log(-lncdf[,i]) - qnorm(p = 1-(1-CI)/2, mean = 0, sd = 1)*SE.cdf[,i])))

    rownames(CI.z[[i]]) <- rownames(CI.P[[i]]) <- rownames(CI.cdf[[i]]) <- paste('RC=',zstar+1,sep='')

    # Within Stratum HR-CI
    lnHR <- diff(lnz[,i])
    dhz <- diff(hz[[i]])

    SE.WS <- sqrt(diag( dhz %*% S %*% t(dhz) ))

    CI.WHR[[i]] <- exp(cbind(
      lnHR - qnorm(p = 1-(1-CI)/2, mean = 0, sd = 1)*SE.WS,
      lnHR,
      lnHR + qnorm(p = 1-(1-CI)/2, mean = 0, sd = 1)*SE.WS))

    # Parameter CI
    h['phi',] <- c(1/a,0,(i%in%mean.stratum)+0)
    h['nu',] <- c(1/(g-a),g/(g-a),0)

    h['p',] <- c(1/a+1/(g-a),-g/(g-a),0)

    h['pstar',] <- h['p',]
    h['pstar',3] <- (i%in%mean.stratum)+0

    SE.par[,i] <- sqrt(diag( h %*% S %*% t(h) ))

    CI.par[[i]] <- exp(cbind(
      log(c(phistar[i],nu[i],p[i],pstar[i])) - qnorm(p = 1-(1-CI)/2, mean = 0, sd = 1)*SE.par[,i],
      log(c(phistar[i],nu[i],p[i],pstar[i])),
      log(c(phistar[i],nu[i],p[i],pstar[i])) + qnorm(p = 1-(1-CI)/2, mean = 0, sd = 1)*SE.par[,i]))
    row.names(CI.par[[i]]) <- c('phistar','nu','p','pstar')

    colnames(CI.par[[i]]) <- colnames(CI.z[[i]]) <-
      colnames(CI.WHR[[i]]) <- colnames(CI.P[[i]]) <- colnames(CI.cdf[[i]]) <-  c('CI-','Estimate','CI+')


  }

  # Across Stratum HR
  if(J==2){
    hz.AS <- cbind(-hz[[1]][,1], hz[[2]][,1],
                   -hz[[1]][,2], hz[[2]][,2],
                   hz[[2]][,3])

    SE.AHR <- sqrt(diag(hz.AS%*%Sigma%*%t(hz.AS)))
    CI.AHR <- exp(cbind(
      lnz[,2]-lnz[,1] - qnorm(p = 1-(1-CI)/2, mean = 0, sd = 1)*SE.AHR,
      lnz[,2]-lnz[,1],
      lnz[,2]-lnz[,1] + qnorm(p = 1-(1-CI)/2, mean = 0, sd = 1)*SE.AHR))
    colnames(CI.AHR) <- c('CI-','Estimate','CI+')


    hP.AS <- cbind(-hcdf[[1]][,1]*lncdf[,1], hcdf[[2]][,1]*lncdf[,2],
                   -hcdf[[1]][,2]*lncdf[,1], hcdf[[2]][,2]*lncdf[,2],
                   hcdf[[2]][,3])

    SE.dP <- sqrt(diag(hP.AS%*%Sigma%*%t(hP.AS)))
    CI.dP <- exp(cbind(
      lncdf[,2]-lncdf[,1] - qnorm(p = 1-(1-CI)/2, mean = 0, sd = 1)*SE.dP,
      lncdf[,2]-lncdf[,1],
      lncdf[,2]-lncdf[,1] + qnorm(p = 1-(1-CI)/2, mean = 0, sd = 1)*SE.dP))
    colnames(CI.AHR) <- colnames(CI.dP) <- c('CI-','Estimate','CI+')
  }else{
    CI.AHR <- 'only implemented for J==2'
  }

  names(CI.par) <- names(CI.z) <- names(CI.cdf) <- names(CI.WHR) <- names(CI.P) <- RE.stratum
  OUT <- list(CI.z=CI.z,CI.par=CI.par,CI.WHR=CI.WHR,CI.AHR=CI.AHR,CI.P=CI.P,
              CI.dP=CI.dP,CI.cdf=CI.cdf)
  attr(OUT,'CI-level') <- CI

  return(OUT)
}
