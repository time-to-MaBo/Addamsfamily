combos <- function(order){

  idx <- list()

  M <- matrix(c(1,0,-1,0), nrow = 1)
  colnames(M) <- c('idx1', 'idx2', 'factor', 'exponent')
  idx[[1]] <- M
  if(order < 2){
    return(idx)
  }
  for(i in 2:order){
    # step1: multiply with -A_{1,0}
    M0 <- M
    M0[,'idx1'] <- M0[,'idx1'] + 1
    M0[,'factor'] <- M0[,'factor']*-1
    # step2: take the derivative of D^{(i-1)}
    # not sorted by idx in A!
    ## step2a: take the "derivative of denominator"
    M1 <- M
    M1[,'factor'] <- M1[,'factor'] * (-M1[,'idx1'])
    M1[,c('idx1','idx2')] <- M1[,c('idx1','idx2')] + 1

    ## step2b: take "derivative of nominator"
    M2 <- M
    M2[,'factor'] <- M2[,'factor'] * M2[,'idx2']
    M2[,'exponent'] <- M2[,'exponent'] + 1
    M2 <- M2[!(M2[,'factor']==0),,drop=FALSE]

    M <- rbind(M0,M1,M2)

    # step 3: find unique terms by idx1,idx2 and add factors
    uni <- unique(M[,c('idx1','idx2','exponent')])
    M.uni <- matrix(nrow = 0, ncol= 4)
    for(j in 1:(dim(uni)[1])){
      tmp <- M[,'idx1'] == uni[j,1] &
        M[,'idx2'] == uni[j,2] &
        M[,'exponent'] == uni[j,3]
      tmp <- which(tmp)
      tmp2 <- M[tmp[1],,drop=FALSE]
      tmp2[,'factor'] <- sum(M[tmp,'factor'])
      M.uni <- rbind(M.uni, tmp2)
    }

    idx[[i]] <- M.uni
    M <- M.uni
  }

  return(idx)

}
