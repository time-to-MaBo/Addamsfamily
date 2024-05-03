z_hat.slim <- function(dist_par, alpha, gamma, lambda_i1, lambda_i2,
                       d_1, d_2, Lambda_i1, Lambda_i2, s_i, n, EZ = TRUE, z_0 = 0:3,
                       max = TRUE, distribution = FALSE, idx, frailty_names,
                       case, type, ag_i, a_i, g_i, PVF_i, ag.idx){

  # Output Preparation:
  Z <- list(z = NULL, g = NULL)
  z_hat <- rep(NA, n)
  para_cont <- matrix(NA, nrow = n, ncol = 7) # discrete cases remain NA.
  colnames(para_cont) <- c('shape', 'rate', 'lambda', 'mu', '-alpha', 'delta', 'theta')
  pdf_cont <- replicate(NA, n = n, simplify = FALSE)
  OUT_cont <- list()

  # Loop Conditions
  discrete <- any(c('-aNB+p','aNB','aB','aP','P') %in% case)
  continuous <- any(c('G','IG','H','cP') %in% case)
  log_dist <- distribution
  log_max <- max

  dist.par_i <- apply(X = dist_par[idx,], MARGIN = 1, FUN = function(x) list(x))
  dist.par_i <- lapply(X = dist.par_i, FUN = function(x) unlist(x))

  if(discrete){
    discrete_case <- case
    discrete_case[case %in% c('H', 'G', 'cP', 'IG')] <- NA
    idx_i <- !is.na(discrete_case)[idx]
    while(log_dist | log_max){

      Z_temp <- condfrailtydist(z_star = z_0, n = n, alpha = alpha, gamma = gamma,
                                d_1 = d_1[idx_i], d_2 = d_2[idx_i],
                                lambda_i1 = lambda_i1[idx_i], lambda_i2 = lambda_i2[idx_i], s_i = s_i[idx_i],
                                idx = idx[idx_i], case = discrete_case,
                                dist.par_i = dist.par_i[idx_i],
                                frailty_names = frailty_names, idx_i = idx_i,
                                ag_i = ag_i[idx_i], a_i = a_i[idx_i], g_i = FALSE,
                                PVF_i = FALSE, ag_idx = which(ag.idx))

      Z$z <- cbind(Z$z, Z_temp$discrete$z)
      Z$g <- cbind(Z$g, Z_temp$discrete$g)
      G <- 1 - rowSums(Z$g)

      z_0 <- (max(z_0)+1):(2*max(z_0))

      log_max <- !all(rowany(G < Z$g), na.rm = TRUE) & max # na.rm f?r continuous distribution
      log_dist <- !all(G < 1e-14, na.rm = TRUE) & distribution
    }

    idx_i <- which(idx %in% which(case %in% c('-aNB+p', 'aNB', 'aB', 'P', 'aP')))
    if(EZ){ #z_hat = E[Z|data]

      z_hat[idx_i] <- rowSums(Z$g*Z$z)

    }else{ # z_hat = conditional mode

      idx_max <- apply(X = Z$g[idx_i,,drop = FALSE], MARGIN = 1, FUN = function(x){
        which(x == max(x))
      })

      z_hat[idx_i] <- mapply(FUN = function(i,j){
        Z$z[i,j]
      },
      i = idx_i, j = idx_max)

    }

  }

  if(continuous){
    cont_case <- case
    cont_case[case %in% c('P', '-aNB+p', 'aNB', 'aB', 'aP', 'G', 'cP')] <- NA
    if(any(c('G', 'cP') %in% cont_case)){

      idx_i <- which( idx %in% which(case %in% c('G','cP')) )
      N <- rep(1, length(idx_i))
      idx_i2 <- c()

      if(any('cP' %in% case)){

        idx_i2 <- which( idx %in% which(case %in% c('cP')) )

        N[idx_i2] <- sapply(X = dist.par_i[idx_i2], FUN = function(x) x['lambda'])
        if(!EZ){
          N[idx_i2] <- sapply(X = dist.par_i[idx_i2], FUN = function(x) floor(x['lambda']))
        }
        dist.par_i[idx_i2] <- mapply(FUN = function(d.p_i,k){
          d.p_i['shape'] <- d.p_i['shape']*k
          d.p_i
        },
        d.p_i = dist.par_i[idx_i2], k = N, SIMPLIFY = FALSE)

      }

      para_cont[idx_i, 'lambda'] <- mapply(FUN = function(x,c){
        if(c){
          x['lambda']
        }else{
          Inf
        }
      },
      x = dist.par_i[idx_i],
      c = (idx_i %in% idx_i2), SIMPLIFY = TRUE
      )

      para_cont[idx_i, 'shape'] <- mapply(FUN = function(x,y,k){
        k*(x['shape'] + y)
      },
      x = dist.par_i[idx_i], y = d_1[idx_i] + d_2[idx_i],
      k = N>0, SIMPLIFY = TRUE
      )

      para_cont[idx_i, 'rate'] <- mapply(FUN = function(x,y){
        x['rate'] + y
      },
      x = dist.par_i[idx_i], y = s_i[idx_i],
      SIMPLIFY = TRUE
      )

      pdf_cont[idx_i] <- mapply(FUN = function(s,r,l){
        function(z){
          if(z == 0 & is.finite(l)){
            # exp(-l) + floor(
            #   pgamma(.Machine$double.xmin, shape = s, rate = r))*(1-exp(-l))
            exp(-l)
          }else{
            dgamma(z, shape = s, rate = r)*(1-exp(-l))
          }
        }
      },
      s = para_cont[idx_i, 'shape'],
      r = para_cont[idx_i, 'rate'], l = para_cont[idx_i, 'lambda']
      )

      if(EZ){
        z_hat[idx_i] <- para_cont[idx_i, 'shape'] / para_cont[idx_i,'rate']
      }else{
        z_hat[idx_i] <- (para_cont[idx_i, 'shape']-1) / para_cont[idx_i,'rate']
        z_hat[z_hat < 0  & !is.na(z_hat)] <- 0 + .Machine$double.xmin*(N>0)
      }

    }

    if(any(c('IG') %in% case)){

      idx_i <- which(idx %in% which(case %in% 'IG'))

      para_cont[idx_i, 'mu'] <- mapply(FUN = function(x){
        x['mu']
      },
      x = dist.par_i[idx_i], SIMPLIFY = TRUE
      )

      para_cont[idx_i, 'shape'] <- mapply(FUN = function(x){
        x['shape']
      },
      x = dist.par_i[idx_i], SIMPLIFY = TRUE
      )

      pdf_cont[idx_i] <- mapply(FUN = function(s,m){
        function(z){
          statmod::dinvgauss(x = z, shape = s, mean = m)
        }
      },
      s = para_cont[idx_i, 'shape'],
      m = para_cont[idx_i, 'mu']
      )

    }

    if(any('H' %in% case)){

      idx_i <- which(idx %in% which(case %in% 'H'))

      para_cont[idx_i, '-alpha'] <- mapply(FUN = function(x){
        x['-alpha']
      },
      x = dist.par_i[idx_i], SIMPLIFY = TRUE
      )

      para_cont[idx_i, 'theta'] <- mapply(FUN = function(x){
        x['theta']
      },
      x = dist.par_i[idx_i], SIMPLIFY = TRUE
      )

      para_cont[idx_i, 'delta'] <- mapply(FUN = function(x){
        x['delta']
      },
      x = dist.par_i[idx_i], SIMPLIFY = TRUE
      )

      pdf_cont[idx_i] <- mapply(FUN = function(a,d,t){
        function(z){
          #Hougaard(z = z, alpha = a, delta = d, theta = t, k = 200)
          'Hougaard densityis currently not implemented'
        }
      },
      a = para_cont[idx_i, '-alpha'],
      d = para_cont[idx_i, 'delta'],
      t = para_cont[idx_i, 'theta']
      )

    }

    idx_i <- which(idx %in% which(case %in% c('IG','H')))
    if(length(idx_i)>0){
      pdf_cont[idx_i] <- condfrailtydist(d_1 = d_1[idx_i], d_2 = d_2[idx_i], n = n,
                                         lambda_i1 = lambda_i1[idx_i], lambda_i2 = lambda_i2[idx_i], s_i = s_i[idx_i],
                                         alpha = alpha, gamma = gamma,
                                         idx = idx[idx_i], case = cont_case,
                                         dist.par_i = dist.par_i[idx_i], pdf_cont = pdf_cont[idx_i],
                                         frailty_names = frailty_names,
                                         ag_i = FALSE, a_i = FALSE, g_i = FALSE,
                                         PVF_i = PVF_i[idx_i])$continuous$pdf
      # hier war vorher noch ein [idx_i] dran

      if(EZ){
        z_hat[idx_i] <- sapply(X = pdf[idx_i],
                               FUN = function(pdf){
                                 stats::integrate(f = function(x) x*pdf(x),
                                           lower = 0, upper = Inf)
                               })
        # klann ich analytisch berechnen mindestens wenn beide ?berlebt haben
      }else{
        pseudohood <- function(z, pdf = pdf_cont[idx_i],scale = -1e-3){
          ph <- sum(log(
            mapply(FUN = function(dist,zi){
              dist(zi)
            },
            dist = pdf, zi = z, SIMPLIFY = TRUE) )) * scale
          if(is.infinite(ph)){
            print('Mode:Infinite value during of pseudohood')
            ph <- sqrt(.Machine$double.xmax) * sign(ph)
          }
          return(ph)
        }
        z_hat[idx_i] <- optimax(par = rep(0.1, length(idx_i)), fn = pseudohood,
                                method = 'L-BFGS-B',
                                lower = .Machine$double.xmin,upper = .Machine$double.xmax)$par
      }

    }

    OUT_cont <- list(pdf = pdf_cont)

  }

  return(z_hat)

}
