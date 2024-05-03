breaks_piece <- function(y, d, strata, stratum, Qstar, style = 'Breslow', by = NULL){
  # style = c('Breslow', 'Intervall'). by is parameter for intervall

  breaks <- list()

  if(length(strata) == 0){

  }else{

    maxy <- stats::aggregate(y, list(strata), max)[,2]

    y <- y[d == 1]
    strata <- strata[d == 1]

    for(i in 1:Qstar){
      y_temp <- y[strata == stratum[i]]
      y_temp <- y_temp[-which(y_temp == max(y_temp))]
      if(style == 'Breslow'){
        breaks[[i]] <- sort( unique(c(0, y_temp, maxy[i])) )
      }
      if(style == 'Intervall'){
        if(by == 'max'){
          breaks[[i]] <- c(0, maxy[i])
        }else{
          temp <- unique(
            c(seq(from = 0, to = max(y_temp), by = by), maxy[i])
          )
          if(min(y_temp)>temp[2]) temp <- c(0,min(y_temp),temp[min(y_temp)<temp])
          breaks[[i]] <- c(0, temp[temp>=by])
        }
      }
    }
  }

  idx <- sapply(X = breaks, FUN = function(x) {
    ifelse(length(x) == 1, T, F)
  })
  if(any(idx)){
    breaks[idx] <- replicate(list(c(0,1)), n = sum(idx))
  }

  return(breaks)
}
