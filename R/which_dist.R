which_dist <- function(par, ag.idx, a.idx, g.idx,PVF.idx, IG.idx,
                       frailty_names=1:length(ag.idx),type,point_mass){

  PVF.idx <- PVF.idx & !IG.idx

  alpha <- par[,'alpha']
  gamma <- par[,'gamma']
  Q <- length(ag.idx)

  d_OUT <- rep(NA, Q)
  names(d_OUT) <- frailty_names
  dist_OUT <- replicate(NA, n = Q)
  names(dist_OUT) <- frailty_names
  par_OUT <- matrix(NA, nrow = Q, ncol = 9)
  colnames(par_OUT) <- c('nu', 'p', 'lambda', 'mu', 'shape', 'rate', '-alpha', 'theta', 'delta')
  rownames(par_OUT) <- frailty_names
  names(type) <- frailty_names

  idx_ag <- ag.idx
  case_N <- PVF.idx & (abs(alpha)/gamma <= 1e-4 & gamma^-1 <= 1e-4)
  case_a <- idx_ag & (gamma > alpha) & (alpha < 0)
  case_b <- g.idx | (PVF.idx & alpha >= 0 & alpha < 1e-4 & gamma^-1 > 1e-4 ) # alpha -> 0
  # gamma[type == 'PVF' & alpha == 0] <- gamma[type == 'PVF' & alpha == 0]^-1 # warum? gamma andersrum parameterisiert aus Gamma Sicht? m?sste falsch sein, da ich sp?ter eh ^-1 nehme!
  case_c <- idx_ag & (gamma > alpha) & (alpha > 0)
  case_d <- a.idx
  case_e <- idx_ag & (alpha > gamma)
  case_H <- PVF.idx & (alpha < 0 & alpha != -0.5) & !case_N
  case_cP <- PVF.idx & (alpha > 0 & (alpha^-1 > 1e-4|gamma^-1 > 1e-4 | abs(gamma/alpha-1) > 1e-4))  & !case_N
  case_IG <- IG.idx | (type == 'PVF' & alpha == -0.5) & !case_N
  case_P <- PVF.idx & (alpha > 0 & alpha^-1 <= 1e-4 & gamma^-1 <= 1e-4 & abs(gamma/alpha-1) <= 1e-4)

  if(any(case_a)){
    # case a
    nu <- 1/(gamma[case_a] - alpha[case_a])
    rnu <- round(nu,3)
    p <- -alpha[case_a]/(gamma[case_a] - alpha[case_a])
    rp <- round(p, 3)
    ra <- round(alpha[case_a], 3)
    d <- paste0('-alpha * NB(size = nu, prob = p) + p, with -alpha =', -ra,
                ', nu =', rnu, ' and p = ', rp, '.')
    new_par <- cbind(nu, p)

    par_OUT[case_a, c('nu', 'p')] <- new_par
    d_OUT[case_a] <- d
    dist_OUT[case_a] <- mapply(FUN = function(s,p){
      function(z_star){
        stats::dnbinom(x = z_star, size = s, prob = p)
      }
    },
    p = new_par[,2], s = new_par[,1])

  }

  if(any(case_c)){
    # case c
    nu <- 1/(gamma[case_c] - alpha[case_c])
    rnu <- round(nu,3)
    p <- alpha[case_c]/(gamma[case_c])
    rp <- round(p, 3)
    ra <- round(alpha[case_c], 3)
    d <- paste('alpha * NB(nu, p), with alpha =', ra,
               ', nu =', rnu, ' and p = ', rp, '.')
    new_par <- cbind(nu, p)

    par_OUT[case_c, c('nu', 'p')] <- new_par
    d_OUT[case_c] <- d
    for(i in 1:sum(case_c)){
      dist_OUT[case_c][i] <- list(function(z_star){
        stats::dnbinom(x = z_star, size = new_par[i,1],
                prob = new_par[i,2])
      })
    }

  }

  if(any(case_e)){
    # case e
    ra <- round(alpha[case_e], 3)
    p <- (alpha[case_e] - gamma[case_e])/alpha[case_e]
    rp <- round(p, 3)
    b <- (alpha[case_e] - gamma[case_e])^-1
    rb <- round(b, 0)
    if(any(b != rb)){
      warning('nu is not a postive integer, i.e.
              the estimates do not result in a valid density!
              Be aware that b is rounded towards a positive integer for this output!')
    }

    d <- paste('alpha * B(prob = p, size = b), with alpha = ', ra,
               'nu = ', rb, 'and p = ', rp, '.
               Note, that b has to be a postive integer! Unrounded Estimate of b is ', b, '.')
    new_par <- cbind(p, rb)
    d_dist <- stats::dbinom

    par_OUT[case_e, c('p', 'nu')] <- new_par
    d_OUT[case_e] <- d
    for(i in 1:sum(case_e)){
      dist_OUT[case_e][i] <- list(function(z_star){
        stats::dbinom(x = z_star, size = new_par[i,2],
               prob = new_par[i,1])
      })
    }

  }

  if(any(case_d)){
    # case d
    ra <- round(alpha[case_d], 3)
    alpha_inv <- alpha[case_d]^-1
    ra_i <- round(alpha_inv, 3)
    d <- paste('alpha * Pois(lambda = alpha^-1), with alpha =', ra,'and lambda = alpha^-1 = ', ra_i)
    new_par <- alpha_inv
    d_dist <- stats::dpois

    par_OUT[case_d, c('lambda')] <- new_par
    d_OUT[case_d] <- d
    for(i in 1:sum(case_d)){
      dist_OUT[case_d][i] <- list(function(z_star){
        stats::dpois(x = z_star, lambda = new_par[i])
      })
    }

  }

  if(any(case_b)){
    # case b
    gr <- round(gamma[case_b], 3)
    gr[PVF.idx[case_b]] <- gr[PVF.idx[case_b]]^-1
    ginv <- gamma[case_b]^-1
    ginv[PVF.idx[case_b]] <- ginv[PVF.idx[case_b]]^-1
    rginv <- round(ginv, 3)
    d <- paste('Gamma(shape, rate), with shape = rate = gamma^-1 =',
               rginv, 'and variance = gamma = ', gr)
    d[PVF.idx[case_b]] <- paste('Gamma(shape, rate), with shape = rate = gamma =',
                                rginv, 'and variance = gamma^-1 = ', gr)
    new_par <- cbind(ginv, ginv)

    par_OUT[case_b, c('shape', 'rate')] <- new_par
    d_OUT[case_b] <- d
    for(i in 1:sum(case_b)){
      dist_OUT[case_b][i] <- list(function(z){
        stats::dgamma(x = z, shape = new_par[i], rate = new_par[i])
      })
    }

  }

  if(any(case_H)){

    ra <- round(-alpha[case_H], 3)
    rg <- round(gamma[case_H], 3)
    delta <- gamma[case_H]^(alpha[case_H]+1)
    rd <- round(delta,3)
    d <- paste('Hougaard(alpha*, delta, theta), with alpha* = ',
               ra, ' and delta = ', rd,
               'and theta = ', rg)
    new_par <- cbind(-alpha[case_H], gamma[case_H], delta)

    par_OUT[case_H, c('-alpha', 'theta','delta')] <- new_par
    d_OUT[case_H] <- d
    for(i in 1:sum(case_H)){
      dist_OUT[case_H][i] <- list(function(z){
        #Hougaard(z = z, alpha = -new_par[i,1], delta = new_par[i,3], theta = new_par[i,2])
        'Hougaard density is currently not implemented'
      })
    }

  }

  if(any(case_IG)){

    rg <- round(gamma[case_IG]/2,3)
    d <- paste('Inverse Gaussian(mu, shape), with  mu = ', 1,
               'and shape = ', rg)
    new_par <- cbind(1, gamma[case_IG]/2)
    par_OUT[case_IG, c('mu', 'shape')] <- new_par
    d_OUT[case_IG] <- d
    for(i in 1:sum(case_IG)){
      dist_OUT[case_IG][i] <- list(function(z_star){
        statmod::dinvgauss(y = z, m = new_par[i,1], s = new_par[i,2]^-1)
      })
    }

  }

  if(any(case_cP)){

    g <- gamma[case_cP]
    a <- alpha[case_cP]
    ra <- round(a,3)
    rg <- round(g,3)
    lam <- gamma[case_cP]/a
    rlam <- round(lam,3)
    d <- paste('Compound Poisson distribution:
                X_i ~ Gamma(shape = ',ra,'rate = ',rg,'); i=1,..,N;
                and N ~ Pois(lambda = ',rlam , ')')
    new_par <- cbind(a, g,lam)
    par_OUT[case_cP, c('shape', 'rate', 'lambda')] <- new_par
    d_OUT[case_cP] <- d
    for(i in 1:sum(case_cP)){
      dist_OUT[case_cP][i] <- list('compund Poisson Gamma. Count variable is random')
    }

  }

  if(any(case_P)){

    lam <- 1
    rlam <- round(lam, 3)
    d <- paste('Poisson(lambda = ',rlam,')')
    new_par <- cbind(lam)

    par_OUT[case_P, c('lambda')] <- new_par
    d_OUT[case_P] <- d
    for(i in 1:sum(case_P)){
      dist_OUT[case_P][i] <- list(function(z){
        stats::dpois(x = z, lambda = new_par[i])
      })
    }

  }

  if(any(case_N)){

    d <- paste('(Z-1)/(sqrt{(alpha + 1)/gamma}) ~ N(0,1)')

    for(i in 1:sum(case_N)){
      dist_OUT[case_N][i] <- list(function(z){
        x <- (z-1)/sqrt((alpha + 1)/gamma)
        stats::dnorm(x = x, mean = 0, sd = 1)
      })
    }
  }

  par_OUT <- cbind(par_OUT,'point mass' = 0)
  if(!is.null(point_mass)){
    #d <- paste(d,'and support shifted by',round(point_mass,3)) # d ist falsch hier, m?sste an d_OUT.
    par_OUT[,'point mass'] <- point_mass
  }
  par_OUT[case_a,'point mass'] <- par_OUT[case_a,'p']

  case <- cbind(case_a,case_b,case_c,case_d,case_e,
                case_H,case_cP,case_IG,case_P, case_N)
  new_case <- rep(NA, nrow = dim(case)[1])
  new_case[case[,1]] <- '-aNB+p'
  new_case[case[,2]] <- 'G'
  new_case[case[,3]] <- 'aNB'
  new_case[case[,4]] <- 'aP'
  new_case[case[,5]] <- 'aB'
  new_case[case[,6]] <- 'H'
  new_case[case[,7]] <- 'cP'
  new_case[case[,8]] <- 'IG'
  new_case[case[,9]] <- 'P'
  new_case[case[,10]] <- 'N'
  names(new_case) <- frailty_names

  OUT <- list(distribution = d_OUT, par = par_OUT, case = new_case, type = type,
              dist_fun = dist_OUT)
  return(OUT)
}
